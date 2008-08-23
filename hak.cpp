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
#include "hak.h"
#include "util.h"
#include "diffs.h"


const char *HAK::Name = "HAK";


bool HAK::checkProperties() {
    if( !HasProperty("tol") )
        return false;
    if (!HasProperty("maxiter") )
        return false;
    if (!HasProperty("verbose") )
        return false;
    if( !HasProperty("doubleloop") )
        return false;
    if( !HasProperty("clusters") )
        return false;
    
    ConvertPropertyTo<double>("tol");
    ConvertPropertyTo<size_t>("maxiter");
    ConvertPropertyTo<size_t>("verbose");
    ConvertPropertyTo<bool>("doubleloop");
    ConvertPropertyTo<ClustersType>("clusters");

    if( HasProperty("loopdepth") )
        ConvertPropertyTo<size_t>("loopdepth");
    else if( Clusters() == ClustersType::LOOP )
        return false;

    return true;
}


void HAK::constructMessages() {
    // Create outer beliefs
    _Qa.clear();
    _Qa.reserve(nr_ORs());
    for( size_t alpha = 0; alpha < nr_ORs(); alpha++ )
        _Qa.push_back( Factor( OR(alpha).vars() ) );

    // Create inner beliefs
    _Qb.clear();
    _Qb.reserve(nr_IRs());
    for( size_t beta = 0; beta < nr_IRs(); beta++ )
        _Qb.push_back( Factor( IR(beta) ) );
    
    // Create messages
    _muab.clear();
    _muab.reserve(nr_Redges());
    _muba.clear();
    _muba.reserve(nr_Redges());
    for( vector<R_edge_t>::const_iterator ab = Redges().begin(); ab != Redges().end(); ab++ ) {
        _muab.push_back( Factor( IR(ab->second) ) );
        _muba.push_back( Factor( IR(ab->second) ) );
    }
}


HAK::HAK(const RegionGraph & rg, const Properties &opts) : DAIAlgRG(rg, opts) {
    assert( checkProperties() );

    constructMessages();
}


void HAK::findLoopClusters( const FactorGraph & fg, set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars ) {
    for( VarSet::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        VarSet ind = fg.delta( *in );
        if( (newcl.size()) >= 2 && (ind >> root) ) {
            allcl.insert( newcl | *in );
        }
        else if( length > 1 )
            findLoopClusters( fg, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


HAK::HAK(const FactorGraph & fg, const Properties &opts) : DAIAlgRG(opts) {
    assert( checkProperties() );

    vector<VarSet> cl;
    if( Clusters() == ClustersType::MIN ) {
        cl = fg.Cliques();
    } else if( Clusters() == ClustersType::DELTA ) {
        for( size_t i = 0; i < fg.nrVars(); i++ )
            cl.push_back(fg.Delta(fg.var(i))); 
    } else if( Clusters() == ClustersType::LOOP ) {
        cl = fg.Cliques();
        set<VarSet> scl;
        for( vector<Var>::const_iterator i0 = fg.vars().begin(); i0 != fg.vars().end(); i0++ ) {
            VarSet i0d = fg.delta(*i0);
            if( LoopDepth() > 1 )
                findLoopClusters( fg, scl, *i0, *i0, LoopDepth() - 1, fg.delta(*i0) );
        }
        for( set<VarSet>::const_iterator c = scl.begin(); c != scl.end(); c++ )
            cl.push_back(*c);
        if( Verbose() >= 3 ) {
            cout << "HAK uses the following clusters: " << endl;
            for( vector<VarSet>::const_iterator cli = cl.begin(); cli != cl.end(); cli++ )
                cout << *cli << endl;
        }
    } else
        throw "Invalid Clusters type";

    RegionGraph rg(fg,cl);
    RegionGraph::operator=(rg);
    constructMessages();

    if( Verbose() >= 3 )
        cout << "HAK regiongraph: " << *this << endl;
}


string HAK::identify() const { 
    stringstream result (stringstream::out);
    result << Name << GetProperties();
    return result.str();
}


void HAK::init( const VarSet &ns ) {
    for( vector<Factor>::iterator alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
        if( alpha->vars() && ns )
            alpha->fill( 1.0 / alpha->stateSpace() );

    for( size_t beta = 0; beta < nr_IRs(); beta++ )
        if( IR(beta) && ns ) {
            _Qb[beta].fill( 1.0 / IR(beta).stateSpace() );
            for( R_nb_cit alpha = nbIR(beta).begin(); alpha != nbIR(beta).end(); alpha++ ) {
                muab(*alpha,beta).fill( 1.0 / IR(beta).stateSpace() );
                muba(beta,*alpha).fill( 1.0 / IR(beta).stateSpace() );
            }
        }
}


void HAK::init() {
    assert( checkProperties() );

    for( vector<Factor>::iterator alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
        alpha->fill( 1.0 / alpha->stateSpace() );

    for( vector<Factor>::iterator beta = _Qb.begin(); beta != _Qb.end(); beta++ )
        beta->fill( 1.0 / beta->stateSpace() );

    for( size_t ab = 0; ab < nr_Redges(); ab++ ) {
        _muab[ab].fill( 1.0 / _muab[ab].stateSpace() );
        _muba[ab].fill( 1.0 / _muba[ab].stateSpace() );
    }
}


double HAK::doGBP() {
    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";
    if( Verbose() >= 3)
        cout << endl;

    clock_t tic = toc();

    // Check whether counting numbers won't lead to problems
    for( size_t beta = 0; beta < nr_IRs(); beta++ )
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
    for( iter = 0; iter < MaxIter() && diffs.max() > Tol(); iter++ ) {
        for( size_t beta = 0; beta < nr_IRs(); beta++ ) {
            for( R_nb_cit alpha = nbIR(beta).begin(); alpha != nbIR(beta).end(); alpha++ )
                muab(*alpha,beta) = _Qa[*alpha].marginal(IR(beta)).divided_by( muba(beta,*alpha) );

            Factor Qb_new;
            for( R_nb_cit alpha = nbIR(beta).begin(); alpha != nbIR(beta).end(); alpha++ )
                Qb_new *= muab(*alpha,beta) ^ (1 / (nbIR(beta).size() + IR(beta).c()));
            Qb_new.normalize( _normtype );
            if( Qb_new.hasNaNs() ) {
                cout << "HAK::doGBP:  Qb_new has NaNs!" << endl;
                return NAN;
            }
//          _Qb[beta] = Qb_new.makeZero(1e-100);    // damping?
            _Qb[beta] = Qb_new;

            for( R_nb_cit alpha = nbIR(beta).begin(); alpha != nbIR(beta).end(); alpha++ ) {
                muba(beta,*alpha) = _Qb[beta].divided_by( muab(*alpha,beta) );

                Factor Qa_new = OR(*alpha);
                for( R_nb_cit gamma = nbOR(*alpha).begin(); gamma != nbOR(*alpha).end(); gamma++ )
                    Qa_new *= muba(*gamma,*alpha);
                Qa_new ^= (1.0 / OR(*alpha).c());
                Qa_new.normalize( _normtype );
                if( Qa_new.hasNaNs() ) {
                    cout << "HAK::doGBP:  Qa_new has NaNs!" << endl;
                    return NAN;
                }
//              _Qa[*alpha] = Qa_new.makeZero(1e-100); // damping?
                _Qa[*alpha] = Qa_new;
            }
        }

        // Calculate new single variable beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            Factor new_belief = belief( var( i ) );
            diffs.push( dist( new_belief, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = new_belief;
        }

        if( Verbose() >= 3 )
            cout << "HAK::doGBP:  maxdiff " << diffs.max() << " after " << iter+1 << " passes" << endl;
    }

    updateMaxDiff( diffs.max() );

    if( Verbose() >= 1 ) {
        if( diffs.max() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
            cout << "HAK::doGBP:  WARNING: not converged within " << MaxIter() << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.max() << endl;
        } else {
            if( Verbose() >= 2 )
                cout << "HAK::doGBP:  ";
            cout << "converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.max();
}


double HAK::doDoubleLoop() {
    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";
    if( Verbose() >= 3)
        cout << endl;

    clock_t tic = toc();

    // Save original outer regions
    vector<FRegion> org_ORs = ORs();

    // Save original inner counting numbers and set negative counting numbers to zero
    vector<double> org_IR_cs( nr_IRs(), 0.0 );
    for( size_t beta = 0; beta < nr_IRs(); beta++ ) {
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

    size_t   outer_maxiter   = MaxIter();
    double  outer_tol       = Tol();
    size_t   outer_verbose   = Verbose();
    double  org_maxdiff     = MaxDiff();

    // Set parameters for inner loop
    MaxIter( 5 );
    Verbose( outer_verbose ? outer_verbose - 1 : 0 );

    size_t outer_iter = 0;
    for( outer_iter = 0; outer_iter < outer_maxiter && diffs.max() > outer_tol; outer_iter++ ) {
        // Calculate new outer regions
        for( size_t alpha = 0; alpha < nr_ORs(); alpha++ ) {
            OR(alpha) = org_ORs[alpha];
            for( R_nb_cit beta = nbOR(alpha).begin(); beta != nbOR(alpha).end(); beta++ )
                OR(alpha) *= _Qb[*beta] ^ ((IR(*beta).c() - org_IR_cs[*beta]) / nbIR(*beta).size());
        }

        // Inner loop
        if( isnan( doGBP() ) )
            return NAN;

        // Calculate new single variable beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            Factor new_belief = belief( var( i ) );
            diffs.push( dist( new_belief, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = new_belief;
        }

        if( Verbose() >= 3 )
            cout << "HAK::doDoubleLoop:  maxdiff " << diffs.max() << " after " << outer_iter+1 << " passes" << endl;
    }

    // restore _maxiter, _verbose and _maxdiff
    MaxIter( outer_maxiter );
    Verbose( outer_verbose );
    MaxDiff( org_maxdiff );

    updateMaxDiff( diffs.max() );

    // Restore original outer regions
    ORs() = org_ORs;

    // Restore original inner counting numbers
    for( size_t beta = 0; beta < nr_IRs(); beta++ )
        IR(beta).c() = org_IR_cs[beta];

    if( Verbose() >= 1 ) {
        if( diffs.max() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
                cout << "HAK::doDoubleLoop:  WARNING: not converged within " << outer_maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.max() << endl;
            } else {
                if( Verbose() >= 3 )
                    cout << "HAK::doDoubleLoop:  ";
                cout << "converged in " << outer_iter << " passes (" << toc() - tic << " clocks)." << endl;
            }
        }

    return diffs.max();
}


double HAK::run() {
    if( DoubleLoop() )
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
    for( size_t beta = 0; beta < nr_IRs(); beta++ )
        result.push_back( Qb(beta) );
    for( size_t alpha = 0; alpha < nr_ORs(); alpha++ )
        result.push_back( Qa(alpha) );
    return result;
}


Complex HAK::logZ() const {
    Complex sum = 0.0;
    for( size_t beta = 0; beta < nr_IRs(); beta++ )
        sum += Complex(IR(beta).c()) * Qb(beta).entropy();
    for( size_t alpha = 0; alpha < nr_ORs(); alpha++ ) {
        sum += Complex(OR(alpha).c()) * Qa(alpha).entropy();
        sum += (OR(alpha).log0() * Qa(alpha)).totalSum();
    }
    return sum;
}
