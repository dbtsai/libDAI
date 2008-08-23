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


#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include "lc.h"
#include "diffs.h"
#include "util.h"
#include "alldai.h"
#include "x2x.h"


using namespace std;


const char *LC::Name = "LC";


bool LC::checkProperties() {
    if( !HasProperty("cavity") )
        return false;
    if( !HasProperty("updates") )
        return false;
    if( !HasProperty("tol") )
        return false;
    if (!HasProperty("maxiter") )
        return false;
    if (!HasProperty("verbose") )
        return false;

    ConvertPropertyTo<CavityType>("cavity");
    ConvertPropertyTo<UpdateType>("updates");
    ConvertPropertyTo<double>("tol");
    ConvertPropertyTo<size_t>("maxiter");
    ConvertPropertyTo<size_t>("verbose");

    if (HasProperty("cavainame") )
        ConvertPropertyTo<string>("cavainame");
    if (HasProperty("cavaiopts") )
        ConvertPropertyTo<Properties>("cavaiopts");
    if( HasProperty("reinit") )
        ConvertPropertyTo<bool>("reinit");
    
    return true;
}


LC::LC(const FactorGraph & fg, const Properties &opts) : DAIAlgFG(fg, opts) {
    assert( checkProperties() );

    // calc iI
    for( size_t i=0; i < nrVars(); i++ ) {
        for( _nb_cit I = nb1(i).begin(); I != nb1(i).end(); I++ ) {
                _iI_type _iI_entry;
                _iI_entry.i = i;
                _iI_entry.I = *I;

                _iI.push_back(_iI_entry);
            }
    }

    // create pancakes
    _pancakes.resize(nrVars());
   
    // create cavitydists
    for( size_t i=0; i < nrVars(); i++ )
        _cavitydists.push_back(Factor(delta(var(i))));

    // create phis
    _phis.reserve(nr_edges());
    for( size_t iI = 0; iI < nr_edges(); iI++ ) {
        size_t i = edge(iI).first;
        size_t I = edge(iI).second;
        _phis.push_back( Factor( factor(I).vars() / var(i) ) );
    }

    // create beliefs
    for( size_t i=0; i < nrVars(); i++ )
        _beliefs.push_back(Factor(var(i)));
}


string LC::identify() const { 
    stringstream result (stringstream::out);
    result << Name << GetProperties();
    return result.str();
}


void LC::CalcBelief (size_t i) {
    _beliefs[i] = _pancakes[i].marginal(var(i));
}


double LC::CalcCavityDist (size_t i, const string &name, const Properties &opts) {
    Factor Bi;
    double maxdiff = 0;

    if( Verbose() >= 2 )
        cout << "Initing cavity " << var(i) << "(" << delta(var(i)).size() << " vars, " << delta(var(i)).stateSpace() << " states)" << endl;

    if( Cavity() == CavityType::UNIFORM )
        Bi = Factor(delta(var(i)));
    else {
        InfAlg *cav = newInfAlg( name, *this, opts );
        cav->makeCavity( var(i) );

        if( Cavity() == CavityType::FULL )
            Bi = calcMarginal( *cav, cav->delta(var(i)), reInit() );
        else if( Cavity() == CavityType::PAIR )
            Bi = calcMarginal2ndO( *cav, cav->delta(var(i)), reInit() );
        else if( Cavity() == CavityType::PAIR2 ) {
            vector<Factor> pairbeliefs = calcPairBeliefsNew( *cav, cav->delta(var(i)), reInit() );
            for( size_t ij = 0; ij < pairbeliefs.size(); ij++ )
                Bi *= pairbeliefs[ij];
        } else if( Cavity() == CavityType::PAIRINT ) {
            Bi = calcMarginal( *cav, cav->delta(var(i)), reInit() );
            
            // Set interactions of order > 2 to zero
            size_t N = delta(var(i)).size();
            double *p = &(*Bi.p().p().begin());
            x2x::p2logp (N, p);
            x2x::logp2w (N, p);
            x2x::fill (N, p, 2, 0.0);
            x2x::w2logp (N, p);
//            x2x::logpnorm (N, p);
            x2x::logp2p (N, p);
        } else if( Cavity() == CavityType::PAIRCUM ) {
            Bi = calcMarginal( *cav, cav->delta(var(i)), reInit() );
            
            // Set cumulants of order > 2 to zero
            size_t N = delta(var(i)).size();
            double *p = &(*Bi.p().p().begin());
            x2x::p2m (N, p);
            x2x::m2c (N, p, N);
            x2x::fill (N, p, 2, 0.0);
            x2x::c2m (N, p, N);
            x2x::m2p (N, p);
        }
        maxdiff = cav->MaxDiff();
        delete cav;
    }
    Bi.normalize( _normtype );
    _cavitydists[i] = Bi;

    return maxdiff;
}


double LC::InitCavityDists (const string &name, const Properties &opts) {
    clock_t tic = toc();

    if( Verbose() >= 1 ) {
        cout << "LC::InitCavityDists:  ";
        if( Cavity() == CavityType::UNIFORM )
            cout << "Using uniform initial cavity distributions" << endl;
        else if( Cavity() == CavityType::FULL )
            cout << "Using full " << name << opts << "...";
        else if( Cavity() == CavityType::PAIR )
            cout << "Using pairwise " << name << opts << "...";
        else if( Cavity() == CavityType::PAIR2 )
            cout << "Using pairwise(new) " << name << opts << "...";
    }

    double maxdiff = 0.0;
    for( size_t i = 0; i < nrVars(); i++ ) {
        double md = CalcCavityDist(i, name, opts);
        if( md > maxdiff )
            maxdiff = md;
    }
    init();

    if( Verbose() >= 1 ) {
        cout << "used " << toc() - tic << " clocks." << endl;
    }

    return maxdiff;
}


long LC::SetCavityDists (vector<Factor> &Q) {
    if( Verbose() >= 1 ) 
        cout << "LC::SetCavityDists:  Setting initial cavity distributions" << endl;
    if( Q.size() != nrVars() )
        return -1;
    for( size_t i = 0; i < nrVars(); i++ ) {
        if( _cavitydists[i].vars() != Q[i].vars() ) {
            return i+1;
        } else
            _cavitydists[i] = Q[i];
    }
    init();
    return 0;
}


void LC::init() {
    for( size_t iI = 0; iI < nr_edges(); iI++ ) {
        if( Updates() == UpdateType::SEQRND )
            _phis[iI].randomize();
        else
            _phis[iI].fill(1.0);
    }
    for( size_t i = 0; i < nrVars(); i++ ) {
        _pancakes[i] = _cavitydists[i];
        
        for( _nb_cit I = nb1(i).begin(); I != nb1(i).end(); I++ ) {
            _pancakes[i] *= factor(*I);
            if( Updates() == UpdateType::SEQRND )
              _pancakes[i] *= _phis[VV2E(i,*I)];
        }
        
        _pancakes[i].normalize( _normtype );

        CalcBelief(i);
    }
}


Factor LC::NewPancake (size_t iI, bool & hasNaNs) {
    size_t i = _iI[iI].i;
    size_t I = _iI[iI].I;
    iI = VV2E(i, I);

    Factor piet = _pancakes[i];

    // recalculate _pancake[i]
    VarSet Ivars = factor(I).vars();
    Factor A_I;
    for( VarSet::const_iterator k = Ivars.begin(); k != Ivars.end(); k++ )
        if( var(i) != *k )
            A_I *= (_pancakes[findVar(*k)] * factor(I).inverse()).part_sum( Ivars / var(i) );
    if( Ivars.size() > 1 )
        A_I ^= (1.0 / (Ivars.size() - 1));
    Factor A_Ii = (_pancakes[i] * factor(I).inverse() * _phis[iI].inverse()).part_sum( Ivars / var(i) );
    Factor quot = A_I.divided_by(A_Ii);

    piet *= quot.divided_by( _phis[iI] ).normalized( _normtype );
    _phis[iI] = quot.normalized( _normtype );

    piet.normalize( _normtype );

    if( piet.hasNaNs() ) {
        cout << "LC::NewPancake(" << iI << "):  has NaNs!" << endl;
        hasNaNs = true;
    }

    return piet;
}


double LC::run() {
    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";
    if( Verbose() >= 2 )
        cout << endl;

    clock_t tic = toc();
    Diffs diffs(nrVars(), 1.0);

    updateMaxDiff( InitCavityDists(GetPropertyAs<string>("cavainame"), GetPropertyAs<Properties>("cavaiopts")) );

    vector<Factor> old_beliefs;
    for(size_t i=0; i < nrVars(); i++ )
        old_beliefs.push_back(belief(i));

    bool hasNaNs = false;
    for( size_t i=0; i < nrVars(); i++ )
        if( _pancakes[i].hasNaNs() ) {
            hasNaNs = true;
            break;
        }
    if( hasNaNs ) {
        cout << "LC::run:  initial _pancakes has NaNs!" << endl;
        return NAN;
    }

    vector<long> update_seq(nr_iI(),0);
    for( size_t k=0; k < nr_iI(); k++ )
        update_seq[k] = k;

    size_t iter=0;

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( iter=0; iter < MaxIter() && diffs.max() > Tol(); iter++ ) {
        // Sequential updates
        if( Updates() == UpdateType::SEQRND )
            random_shuffle( update_seq.begin(), update_seq.end() );
        
        for( size_t t=0; t < nr_iI(); t++ ) {
            long iI = update_seq[t];
            long i = _iI[iI].i;
            _pancakes[i] = NewPancake(iI, hasNaNs);
            if( hasNaNs )
                return NAN;
            CalcBelief(i);
        }

        // compare new beliefs with old ones
        for(size_t i=0; i < nrVars(); i++ ) {
            diffs.push( dist( belief(i), old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = belief(i);
        }

        if( Verbose() >= 3 )
            cout << "LC::run:  maxdiff " << diffs.max() << " after " << iter+1 << " passes" << endl;
    }

    updateMaxDiff( diffs.max() );

    if( Verbose() >= 1 ) {
        if( diffs.max() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
                cout << "LC::run:  WARNING: not converged within " << MaxIter() << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.max() << endl;
        } else {
            if( Verbose() >= 2 )
                cout << "LC::run:  ";
                cout << "converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.max();
}
