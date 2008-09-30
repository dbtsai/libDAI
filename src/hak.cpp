/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

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
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<double>("damping");
    else
        props.damping = 0.0;
}


PropertySet HAK::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "doubleloop", props.doubleloop );
    opts.Set( "clusters", props.clusters );
    opts.Set( "loopdepth", props.loopdepth );
    opts.Set( "damping", props.damping );
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
    s << "loopdepth=" << props.loopdepth << ",";
    s << "damping=" << props.damping << "]";
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


HAK::HAK( const RegionGraph &rg, const PropertySet &opts ) : DAIAlgRG(rg), _Qa(), _Qb(), _muab(), _muba(), _maxdiff(0.0), _iters(0U), props() {
    setProperties( opts );

    constructMessages();
}


void HAK::findLoopClusters( const FactorGraph & fg, std::set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars ) {
    for( VarSet::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        VarSet ind = fg.delta( fg.findVar( *in ) );
        if( (newcl.size()) >= 2 && ind.contains( root ) ) {
            allcl.insert( newcl | *in );
        }
        else if( length > 1 )
            findLoopClusters( fg, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


HAK::HAK(const FactorGraph & fg, const PropertySet &opts) : DAIAlgRG(), _Qa(), _Qb(), _muab(), _muba(), _maxdiff(0.0), _iters(0U), props() {
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
            cout << Name << " uses the following clusters: " << endl;
            for( vector<VarSet>::const_iterator cli = cl.begin(); cli != cl.end(); cli++ )
                cout << *cli << endl;
        }
    } else
        DAI_THROW(INTERNAL_ERROR);

    RegionGraph rg(fg,cl);
    RegionGraph::operator=(rg);
    constructMessages();

    if( props.verbose >= 3 )
        cout << Name << " regiongraph: " << *this << endl;
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
                muab( alpha, _beta ).fill( 1.0 );
                muba( alpha, _beta ).fill( 1.0 );
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

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( _iters = 0; _iters < props.maxiter && diffs.maxDiff() > props.tol; _iters++ ) {
        for( size_t beta = 0; beta < nrIRs(); beta++ ) {
            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                muab( alpha, _beta ) = _Qa[alpha].marginal(IR(beta)).divided_by( muba(alpha,_beta) );
                /* TODO: INVESTIGATE THIS PROBLEM
                 *
                 * In some cases, the muab's can have very large entries because the muba's have very
                 * small entries. This may cause NANs later on (e.g., multiplying large quantities may
                 * result in +inf; normalization then tries to calculate inf / inf which is NAN). 
                 * A fix of this problem would consist in normalizing the messages muab.
                 * However, it is not obvious whether this is a real solution, because it has a
                 * negative performance impact and the NAN's seem to be a symptom of a fundamental
                 * numerical unstability.
                 */
                 muab(alpha,_beta).normalize(); 
            }

            Factor Qb_new;
            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                Qb_new *= muab(alpha,_beta) ^ (1 / (nbIR(beta).size() + IR(beta).c()));
            }

            Qb_new.normalize();
            if( Qb_new.hasNaNs() ) {
                // TODO: WHAT TO DO IN THIS CASE?
                cout << Name << "::doGBP:  Qb_new has NaNs!" << endl;
                return 1.0;
            }
            /* TODO: WHAT IS THE PURPOSE OF THE FOLLOWING CODE?
             *
             *   _Qb[beta] = Qb_new.makeZero(1e-100);
             */

            if( props.doubleloop || props.damping == 0.0 )
                _Qb[beta] = Qb_new; // no damping for double loop
            else
                _Qb[beta] = (Qb_new^(1.0 - props.damping)) * (_Qb[beta]^props.damping);

            foreach( const Neighbor &alpha, nbIR(beta) ) {
                size_t _beta = alpha.dual;
                muba(alpha,_beta) = _Qb[beta].divided_by( muab(alpha,_beta) );

                /* TODO: INVESTIGATE WHETHER THIS HACK (INVENTED BY KEES) TO PREVENT NANS MAKES SENSE 
                 *
                 *   muba(beta,*alpha).makePositive(1e-100);
                 *
                 */

                Factor Qa_new = OR(alpha);
                foreach( const Neighbor &gamma, nbOR(alpha) )
                    Qa_new *= muba(alpha,gamma.iter);
                Qa_new ^= (1.0 / OR(alpha).c());
                Qa_new.normalize();
                if( Qa_new.hasNaNs() ) {
                    cout << Name << "::doGBP:  Qa_new has NaNs!" << endl;
                    return 1.0;
                }
                /* TODO: WHAT IS THE PURPOSE OF THE FOLLOWING CODE?
                 *
                 *   _Qb[beta] = Qb_new.makeZero(1e-100);
                 */

                if( props.doubleloop || props.damping == 0.0 )
                    _Qa[alpha] = Qa_new; // no damping for double loop
                else
                    // FIXME: GEOMETRIC DAMPING IS SLOW!
                _Qa[alpha] = (Qa_new^(1.0 - props.damping)) * (_Qa[alpha]^props.damping);
            }
        }

        // Calculate new single variable beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            Factor new_belief = belief( var( i ) );
            diffs.push( dist( new_belief, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = new_belief;
        }

        if( props.verbose >= 3 )
            cout << Name << "::doGBP:  maxdiff " << diffs.maxDiff() << " after " << _iters+1 << " passes" << endl;
    }

    if( diffs.maxDiff() > _maxdiff )
        _maxdiff = diffs.maxDiff();

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
            cout << Name << "::doGBP:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( props.verbose >= 2 )
                cout << Name << "::doGBP:  ";
            cout << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
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

    size_t outer_maxiter   = props.maxiter;
    double outer_tol       = props.tol;
    size_t outer_verbose   = props.verbose;
    double org_maxdiff     = _maxdiff;

    // Set parameters for inner loop
    props.maxiter = 5;
    props.verbose = outer_verbose ? outer_verbose - 1 : 0;

    size_t outer_iter = 0;
    size_t total_iter = 0;
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

        total_iter += Iterations();

        if( props.verbose >= 3 )
            cout << Name << "::doDoubleLoop:  maxdiff " << diffs.maxDiff() << " after " << total_iter << " passes" << endl;
    }

    // restore _maxiter, _verbose and _maxdiff
    props.maxiter = outer_maxiter;
    props.verbose = outer_verbose;
    _maxdiff = org_maxdiff;

    _iters = total_iter;
    if( diffs.maxDiff() > _maxdiff )
        _maxdiff = diffs.maxDiff();

    // Restore original outer regions
    ORs = org_ORs;

    // Restore original inner counting numbers
    for( size_t beta = 0; beta < nrIRs(); ++beta )
        IR(beta).c() = org_IR_cs[beta];

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
                cout << Name << "::doDoubleLoop:  WARNING: not converged within " << outer_maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << diffs.maxDiff() << endl;
            } else {
                if( props.verbose >= 3 )
                    cout << Name << "::doDoubleLoop:  ";
                cout << "converged in " << total_iter << " passes (" << toc() - tic << " seconds)." << endl;
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
