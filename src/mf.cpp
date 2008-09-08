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
#include <sstream>
#include <map>
#include <set>
#include <dai/mf.h>
#include <dai/diffs.h>
#include <dai/util.h>


namespace dai {


using namespace std;


const char *MF::Name = "MF";


bool MF::checkProperties() {
    if( !HasProperty("tol") )
        return false;
    if (!HasProperty("maxiter") )
        return false;
    if (!HasProperty("verbose") )
        return false;
    
    ConvertPropertyTo<double>("tol");
    ConvertPropertyTo<size_t>("maxiter");
    ConvertPropertyTo<size_t>("verbose");

    return true;
}


void MF::Regenerate() {
//    DAIAlgFG::Regenerate();

    // clear beliefs
    _beliefs.clear();
    _beliefs.reserve( nrVars() );

    // create beliefs
    for( size_t i = 0; i < nrVars(); ++i )
        _beliefs.push_back(Factor(var(i)));
}


string MF::identify() const { 
    stringstream result (stringstream::out);
    result << Name << GetProperties();
    return result.str();
}


void MF::init() {
    assert( checkProperties() );

    for( vector<Factor>::iterator qi = _beliefs.begin(); qi != _beliefs.end(); qi++ )
        qi->fill(1.0);
}


double MF::run() {
    clock_t tic = toc();

    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";

    size_t pass_size = _beliefs.size();
    Diffs diffs(pass_size * 3, 1.0);

    size_t t=0;
    for( t=0; t < (MaxIter()*pass_size) && diffs.max() > Tol(); t++ ) {
        // choose random Var i
        size_t i = (size_t) (nrVars() * rnd_uniform());

        Factor jan;
        Factor piet;
        foreach( const Neighbor &I, nbV(i) ) {
            Factor henk;
            foreach( const Neighbor &j, nbF(I) ) // for all j in I \ i
                if( j != i )
                    henk *= _beliefs[j];
            piet = factor(I).log0();
            piet *= henk;
            piet = piet.part_sum(var(i));
            piet = piet.exp();
            jan *= piet; 
        }

        jan.normalize( _normtype );

        if( jan.hasNaNs() ) {
            cout << "MF::run():  ERROR: jan has NaNs!" << endl;
            return NAN;
        }

        diffs.push( dist( jan, _beliefs[i], Prob::DISTLINF ) );

        _beliefs[i] = jan;
    }

    updateMaxDiff( diffs.max() );

    if( Verbose() >= 1 ) {
        if( diffs.max() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
            cout << "MF::run:  WARNING: not converged within " << MaxIter() << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.max() << endl;
        } else {
            if( Verbose() >= 2 )
                cout << "MF::run:  ";
            cout << "converged in " << t / pass_size << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.max();
}


Factor MF::beliefV (size_t i) const {
    Factor piet;
    piet = _beliefs[i];
    piet.normalize( Prob::NORMPROB );
    return(piet);
}


Factor MF::belief (const VarSet &ns) const {
    if( ns.size() == 1 )
        return belief( *(ns.begin()) );
    else {
        assert( ns.size() == 1 );
        return Factor();
    }
}


Factor MF::belief (const Var &n) const {
    return( beliefV( findVar( n ) ) );
}


vector<Factor> MF::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( beliefV(i) );
    return result;
}


Complex MF::logZ() const {
    Complex sum = 0.0;
    
    for(size_t i=0; i < nrVars(); i++ )
        sum -= beliefV(i).entropy();
    for(size_t I=0; I < nrFactors(); I++ ) {
        Factor henk;
        foreach( const Neighbor &j, nbF(I) )  // for all j in I
            henk *= _beliefs[j];
        henk.normalize( Prob::NORMPROB );
        Factor piet;
        piet = factor(I).log0();
        piet *= henk;
        sum -= Complex( piet.totalSum() );
    }

    return -sum;
}


void MF::init( const VarSet &ns ) {
    for( size_t i = 0; i < nrVars(); i++ ) {
        if( ns && var(i) )
            _beliefs[i].fill( 1.0 );
    }
}


} // end of namespace dai
