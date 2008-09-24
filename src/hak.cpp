/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
    This file is part of libDAI.

    libDAI is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    libDAI is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libDAI; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <map>
#include <dai/hak.h>
#include <dai/util.h>
#include <dai/diffs.h>
#include <dai/exceptions.h>


namespace dai {


using namespace std;


const char *HAK::Name = "HAK";


void HAK::setProperties( const PropertySet &opts ) {
    assert( opts.hasKey("tol") );
    assert( opts.hasKey("maxiter") );
    assert( opts.hasKey("verbose") );
    assert( opts.hasKey("doubleloop") );
    assert( opts.hasKey("clusters") );
    
    props.tol = opts.getStringAs<double>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.verbose = opts.getStringAs<size_t>("verbose");
    props.doubleloop = opts.getStringAs<bool>("doubleloop");
    props.clusters = opts.getStringAs<Properties::ClustersType>("clusters");

    if( opts.hasKey("loopdepth") )
        props.loopdepth = opts.getStringAs<size_t>("loopdepth");
    else
        assert( props.clusters != Properties::ClustersType::LOOP );
}


PropertySet HAK::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "doubleloop", props.doubleloop );
    opts.Set( "clusters", props.clusters );
    opts.Set( "loopdepth", props.loopdepth );
    return opts;
}


string HAK::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "doubleloop=" << props.doubleloop << ",";
    s << "clusters=" << props.clusters << ",";
    s << "loopdepth=" << props.loopdepth << "]";
    return s.str();
}


void HAK::constructMessages() {
    // Create outer beliefs
    _Qa.clear();
    _Qa.reserve(nrORs());
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        _Qa.push_back( Factor( OR(alpha).vars() ) );

    // Create inner beliefs
    _Qb.clear();
    _Qb.reserve(nrIRs());
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        _Qb.push_back( Factor( IR(beta) ) );
    
    // Create messages
    _muab.clear();
    _muab.reserve( nrORs() );
    _muba.clear();
    _muba.reserve( nrORs() );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        _muab.push_back( vector<Factor>() );
        _muba.push_back( vector<Factor>() );
        _muab[alpha].reserve( nbOR(alpha).size() );
        _muba[alpha].reserve( nbOR(alpha).size() );
        foreach( const Neighbor &beta, nbOR(alpha) ) {
            _muab[alpha].push_back( Factor( IR(beta) ) );
            _muba[alpha].push_back( Factor( IR(beta) ) );
        }
    }
}


HAK::HAK(const RegionGraph & rg, const PropertySet &opts ) : DAIAlgRG(rg) {
    setProperties( opts );

    constructMessages();
}


void HAK::findLoopClusters( const FactorGraph & fg, std::set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars ) {
    for( VarSet::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        VarSet ind = fg.delta( fg.findVar( *in ) );
        if( (newcl.size()) >= 2 && (ind >> root) ) {
            allcl.insert( newcl | *in );
        }
        else if( length > 1 )
            findLoopClusters( fg, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


HAK::HAK(const FactorGraph & fg, const PropertySet &opts) : DAIAlgRG(), props(), maxdiff(0.0) {
    setProperties( opts );

    vector<VarSet> cl;
    if( props.clusters == Properties::ClustersType::MIN ) {
        cl = fg.Cliques();
    } else if( props.clusters == Properties::ClustersType::DELTA ) {
        for( size_t i = 0; i < fg.nrVars(); i++ )
            cl.push_back(fg.Delta(i)); 
    } else if( props.clusters == Properties::ClustersType::LOOP ) {
        cl = fg.Cliques();
        set<VarSet> scl;
        for( size_t i0 = 0; i0 < fg.nrVars(); i0++ ) {
            VarSet i0d = fg.delta(i0);
            if( props.loopdepth > 1 )
                findLoopClusters( fg, scl, fg.var(i0), fg.var(i0), props.loopdepth - 1, fg.delta(i0) );
        }
        for( set<VarSet>::const_iterator c = scl.begin(); c != scl.end(); c++ )
            cl.push_back(*c);
        if( props.verbose >= 3 ) {
            cout << "HAK uses the following clusters: " << endl;
            for( vector<VarSet>::const_iterator cli = cl.begin(); cli != cl.end(); cli++ )
                cout << *cli << endl;
        }
    } else
        DAI_THROW(INTERNAL_ERROR);

    RegionGraph rg(fg,cl);
    RegionGraph::operator=(rg);
    constructMessages();

    if( props.verbose >= 3 )
        cout << "HAK regiongraph: " << *this << endl;
}


string HAK::identify() const { 
    return string(Name) + printProperties();
}


void HAK::init( const VarSet &ns ) {
    for( vector<Factor>::iterator alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
        if( alpha->vars().intersects( ns ) )
            alpha->fill( 1.0 / alpha->states() );

    for( size_t beta = 0; beta < nrIRs(); beta++ )
        if( IR(beta).intersects( ns ) ) {
            _Qb[beta].fill( 1.0 );
            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                muab( alpha, _beta ).fill( 1.0 / IR(beta).states() );
                muba( alpha, _beta ).fill( 1.0 / IR(beta).states() );
            }
        }
}


void HAK::init() {
    for( vector<Factor>::iterator alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
        alpha->fill( 1.0 / alpha->states() );

    for( vector<Factor>::iterator beta = _Qb.begin(); beta != _Qb.end(); beta++ )
        beta->fill( 1.0 / beta->states() );

    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        foreach( const Neighbor &beta, nbOR(alpha) ) {
            size_t _beta = beta.iter;
            muab( alpha, _beta ).fill( 1.0 / muab( alpha, _beta ).states() );
            muba( alpha, _beta ).fill( 1.0 / muab( alpha, _beta ).states() );
        }
}


double HAK::doGBP() {
    if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";
    if( props.verbose >= 3)
        cout << endl;

    double tic = toc();

    // Check whether counting numbers won't lead to problems
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        assert( nbIR(beta).size() + IR(beta).c() != 0.0 );

    // Keep old beliefs to check convergence
    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        old_beliefs.push_back( belief( var(i) ) );

    // Differences in single node beliefs
    Diffs diffs(nrVars(), 1.0);

    size_t iter = 0;
    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( iter = 0; iter < props.maxiter && diffs.maxDiff() > props.tol; iter++ ) {
        for( size_t beta = 0; beta < nrIRs(); beta++ ) {
            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                muab( alpha, _beta ) = _Qa[alpha].marginal(IR(beta)).divided_by( muba(alpha,_beta) );
            }

            Factor Qb_new;
            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                Qb_new *= muab(alpha,_beta) ^ (1 / (nbIR(beta).size() + IR(beta).c()));
            }

            Qb_new.normalize( Prob::NORMPROB );
            if( Qb_new.hasNaNs() ) {
                cout << "HAK::doGBP:  Qb_new has NaNs!" << endl;
                return 1.0;
            }
//          _Qb[beta] = Qb_new.makeZero(1e-100);    // damping?
            _Qb[beta] = Qb_new;

            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;

                muba(alpha,_beta) = _Qb[beta].divided_by( muab(alpha,_beta) );

                Factor Qa_new = OR(alpha);
                foreach( const Neighbor &gamma, nbOR(alpha) )
                    Qa_new *= muba(alpha,gamma.iter);
                Qa_new ^= (1.0 / OR(alpha).c());
                Qa_new.normalize( Prob::NORMPROB );
                if( Qa_new.hasNaNs() ) {
                    cout << "HAK::doGBP:  Qa_new has NaNs!" << endl;
                    return 1.0;
                }
//              _Qa[alpha] = Qa_new.makeZero(1e-100); // damping?
                _Qa[alpha] = Qa_new;
            }
        }

        // Calculate new single variable beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            Factor new_belief = belief( var( i ) );
            diffs.push( dist( new_belief, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = new_belief;
        }

        if( props.verbose >= 3 )
            cout << "HAK::doGBP:  maxdiff " << diffs.maxDiff() << " after " << iter+1 << " passes" << endl;
    }

    if( diffs.maxDiff() > maxdiff )
        maxdiff = diffs.maxDiff();

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
            cout << "HAK::doGBP:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( props.verbose >= 2 )
                cout << "HAK::doGBP:  ";
            cout << "converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.maxDiff();
}


double HAK::doDoubleLoop() {
    if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";
    if( props.verbose >= 3)
        cout << endl;

    double tic = toc();

    // Save original outer regions
    vector<FRegion> org_ORs = ORs;

    // Save original inner counting numbers and set negative counting numbers to zero
    vector<double> org_IR_cs( nrIRs(), 0.0 );
    for( size_t beta = 0; beta < nrIRs(); beta++ ) {
        org_IR_cs[beta] = IR(beta).c();
        if( IR(beta).c() < 0.0 )
            IR(beta).c() = 0.0;
    }

    // Keep old beliefs to check convergence
    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        old_beliefs.push_back( belief( var(i) ) );

    // Differences in single node beliefs
    Diffs diffs(nrVars(), 1.0);

    size_t  outer_maxiter   = props.maxiter;
    double  outer_tol       = props.tol;
    size_t  outer_verbose   = props.verbose;
    double  org_maxdiff     = maxdiff;

    // Set parameters for inner loop
    props.maxiter = 5;
    props.verbose = outer_verbose ? outer_verbose - 1 : 0;

    size_t outer_iter = 0;
    for( outer_iter = 0; outer_iter < outer_maxiter && diffs.maxDiff() > outer_tol; outer_iter++ ) {
        // Calculate new outer regions
        for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
            OR(alpha) = org_ORs[alpha];
            foreach( const Neighbor &beta, nbOR(alpha) )
                OR(alpha) *= _Qb[beta] ^ ((IR(beta).c() - org_IR_cs[beta]) / nbIR(beta).size());
        }

        // Inner loop
        if( isnan( doGBP() ) )
            return 1.0;

        // Calculate new single variable beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor new_belief = belief( var( i ) );
            diffs.push( dist( new_belief, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = new_belief;
        }

        if( props.verbose >= 3 )
            cout << "HAK::doDoubleLoop:  maxdiff " << diffs.maxDiff() << " after " << outer_iter+1 << " passes" << endl;
    }

    // restore _maxiter, _verbose and _maxdiff
    props.maxiter = outer_maxiter;
    props.verbose = outer_verbose;
    maxdiff = org_maxdiff;

    if( diffs.maxDiff() > maxdiff )
        maxdiff = diffs.maxDiff();

    // Restore original outer regions
    ORs = org_ORs;

    // Restore original inner counting numbers
    for( size_t beta = 0; beta < nrIRs(); ++beta )
        IR(beta).c() = org_IR_cs[beta];

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
                cout << "HAK::doDoubleLoop:  WARNING: not converged within " << outer_maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.maxDiff() << endl;
            } else {
                if( props.verbose >= 3 )
                    cout << "HAK::doDoubleLoop:  ";
                cout << "converged in " << outer_iter << " passes (" << toc() - tic << " clocks)." << endl;
            }
        }

    return diffs.maxDiff();
}


double HAK::run() {
    if( props.doubleloop )
        return doDoubleLoop();
    else
        return doGBP();
}


Factor HAK::belief( const VarSet &ns ) const {
    vector<Factor>::const_iterator beta;
    for( beta = _Qb.begin(); beta != _Qb.end(); beta++ )
        if( beta->vars() >> ns )
            break;
    if( beta != _Qb.end() )
        return( beta->marginal(ns) );
    else {
        vector<Factor>::const_iterator alpha;
        for( alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
            if( alpha->vars() >> ns )
                break;
        assert( alpha != _Qa.end() );
        return( alpha->marginal(ns) );
    }
}


Factor HAK::belief( const Var &n ) const {
    return belief( (VarSet)n );
}


vector<Factor> HAK::beliefs() const {
    vector<Factor> result;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        result.push_back( Qb(beta) );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        result.push_back( Qa(alpha) );
    return result;
}


Real HAK::logZ() const {
    Real sum = 0.0;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        sum += IR(beta).c() * Qb(beta).entropy();
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        sum += OR(alpha).c() * Qa(alpha).entropy();
        sum += (OR(alpha).log0() * Qa(alpha)).totalSum();
    }
    return sum;
}


} // end of namespace dai
