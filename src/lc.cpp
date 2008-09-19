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
#include <dai/lc.h>
#include <dai/diffs.h>
#include <dai/util.h>
#include <dai/alldai.h>
#include <dai/x2x.h>


namespace dai {


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

    // create pancakes
    _pancakes.resize(nrVars());
   
    // create cavitydists
    for( size_t i=0; i < nrVars(); i++ )
        _cavitydists.push_back(Factor(delta(i)));

    // create phis
    _phis.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ ) {
        _phis.push_back( vector<Factor>() );
        _phis[i].reserve( nbV(i).size() );
        foreach( const Neighbor &I, nbV(i) )
            _phis[i].push_back( Factor( factor(I).vars() / var(i) ) );
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


double LC::CalcCavityDist (size_t i, const std::string &name, const Properties &opts) {
    Factor Bi;
    double maxdiff = 0;

    if( Verbose() >= 2 )
        cout << "Initing cavity " << var(i) << "(" << delta(i).size() << " vars, " << delta(i).states() << " states)" << endl;

    if( Cavity() == CavityType::UNIFORM )
        Bi = Factor(delta(i));
    else {
        InfAlg *cav = newInfAlg( name, *this, opts );
        cav->makeCavity( i );

        if( Cavity() == CavityType::FULL )
            Bi = calcMarginal( *cav, cav->fg().delta(i), reInit() );
        else if( Cavity() == CavityType::PAIR )
            Bi = calcMarginal2ndO( *cav, cav->fg().delta(i), reInit() );
        else if( Cavity() == CavityType::PAIR2 ) {
            vector<Factor> pairbeliefs = calcPairBeliefsNew( *cav, cav->fg().delta(i), reInit() );
            for( size_t ij = 0; ij < pairbeliefs.size(); ij++ )
                Bi *= pairbeliefs[ij];
        } else if( Cavity() == CavityType::PAIRINT ) {
            Bi = calcMarginal( *cav, cav->fg().delta(i), reInit() );
            
            // Set interactions of order > 2 to zero
            size_t N = delta(i).size();
            Real *p = &(*Bi.p().p().begin());
            x2x::p2logp (N, p);
            x2x::logp2w (N, p);
            x2x::fill (N, p, 2, 0.0);
            x2x::w2logp (N, p);
//            x2x::logpnorm (N, p);
            x2x::logp2p (N, p);
        } else if( Cavity() == CavityType::PAIRCUM ) {
            Bi = calcMarginal( *cav, cav->fg().delta(i), reInit() );
            
            // Set cumulants of order > 2 to zero
            size_t N = delta(i).size();
            Real *p = &(*Bi.p().p().begin());
            x2x::p2m (N, p);
            x2x::m2c (N, p, N);
            x2x::fill (N, p, 2, 0.0);
            x2x::c2m (N, p, N);
            x2x::m2p (N, p);
        }
        maxdiff = cav->MaxDiff();
        delete cav;
    }
    Bi.normalize( Prob::NORMPROB );
    _cavitydists[i] = Bi;

    return maxdiff;
}


double LC::InitCavityDists (const std::string &name, const Properties &opts) {
    double tic = toc();

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


long LC::SetCavityDists( std::vector<Factor> &Q ) {
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
    for( size_t i = 0; i < nrVars(); ++i )
        foreach( const Neighbor &I, nbV(i) )
            if( Updates() == UpdateType::SEQRND )
                _phis[i][I.iter].randomize();
            else
                _phis[i][I.iter].fill(1.0);
    for( size_t i = 0; i < nrVars(); i++ ) {
        _pancakes[i] = _cavitydists[i];
        
        foreach( const Neighbor &I, nbV(i) ) {
            _pancakes[i] *= factor(I);
            if( Updates() == UpdateType::SEQRND )
              _pancakes[i] *= _phis[i][I.iter];
        }
        
        _pancakes[i].normalize( Prob::NORMPROB );

        CalcBelief(i);
    }
}


Factor LC::NewPancake (size_t i, size_t _I, bool & hasNaNs) {
    size_t I = nbV(i)[_I];
    Factor piet = _pancakes[i];

    // recalculate _pancake[i]
    VarSet Ivars = factor(I).vars();
    Factor A_I;
    for( VarSet::const_iterator k = Ivars.begin(); k != Ivars.end(); k++ )
        if( var(i) != *k )
            A_I *= (_pancakes[findVar(*k)] * factor(I).inverse()).part_sum( Ivars / var(i) );
    if( Ivars.size() > 1 )
        A_I ^= (1.0 / (Ivars.size() - 1));
    Factor A_Ii = (_pancakes[i] * factor(I).inverse() * _phis[i][_I].inverse()).part_sum( Ivars / var(i) );
    Factor quot = A_I.divided_by(A_Ii);

    piet *= quot.divided_by( _phis[i][_I] ).normalized( Prob::NORMPROB );
    _phis[i][_I] = quot.normalized( Prob::NORMPROB );

    piet.normalize( Prob::NORMPROB );

    if( piet.hasNaNs() ) {
        cout << "LC::NewPancake(" << i << ", " << _I << "):  has NaNs!" << endl;
        hasNaNs = true;
    }

    return piet;
}


double LC::run() {
    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";
    if( Verbose() >= 2 )
        cout << endl;

    double tic = toc();
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
        return -1.0;
    }

    size_t nredges = nrEdges();
    vector<Edge> update_seq;
    update_seq.reserve( nredges );
    for( size_t i = 0; i < nrVars(); ++i )
        foreach( const Neighbor &I, nbV(i) )
            update_seq.push_back( Edge( i, I.iter ) );

    size_t iter = 0;

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( iter=0; iter < MaxIter() && diffs.maxDiff() > Tol(); iter++ ) {
        // Sequential updates
        if( Updates() == UpdateType::SEQRND )
            random_shuffle( update_seq.begin(), update_seq.end() );
        
        for( size_t t=0; t < nredges; t++ ) {
            size_t i = update_seq[t].first;
            size_t _I = update_seq[t].second;
            _pancakes[i] = NewPancake( i, _I, hasNaNs);
            if( hasNaNs )
                return -1.0;
            CalcBelief( i );
        }

        // compare new beliefs with old ones
        for(size_t i=0; i < nrVars(); i++ ) {
            diffs.push( dist( belief(i), old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = belief(i);
        }

        if( Verbose() >= 3 )
            cout << "LC::run:  maxdiff " << diffs.maxDiff() << " after " << iter+1 << " passes" << endl;
    }

    updateMaxDiff( diffs.maxDiff() );

    if( Verbose() >= 1 ) {
        if( diffs.maxDiff() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
                cout << "LC::run:  WARNING: not converged within " << MaxIter() << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( Verbose() >= 2 )
                cout << "LC::run:  ";
                cout << "converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.maxDiff();
}


} // end of namespace dai
