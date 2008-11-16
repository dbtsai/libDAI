/*  Copyright (C) 2008  Frederik Eaton [frederik at ofb dot net]

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
#include <algorithm>
#include <dai/gibbs.h>
#include <dai/util.h>
#include <dai/properties.h>


namespace dai {


using namespace std;


const char *Gibbs::Name = "GIBBS";


void Gibbs::setProperties( const PropertySet &opts ) {
    assert( opts.hasKey("iters") );
    props.iters = opts.getStringAs<size_t>("iters");

    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
}


PropertySet Gibbs::getProperties() const {
    PropertySet opts;
    opts.Set( "iters", props.iters );
    opts.Set( "verbose", props.verbose );
    return opts;
}


string Gibbs::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "iters=" << props.iters << ",";
    s << "verbose=" << props.verbose << "]";
    return s.str();
}


void Gibbs::construct() {
    _var_counts.clear();
    _var_counts.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        _var_counts.push_back( _count_t( var(i).states(), 0 ) );
    
    _factor_counts.clear();
    _factor_counts.reserve( nrFactors() );
    for( size_t I = 0; I < nrFactors(); I++ )
        _factor_counts.push_back( _count_t( factor(I).states(), 0 ) );

    _sample_count = 0;

    _factor_entries.clear();
    _factor_entries.resize( nrFactors(), 0 );

    _state.clear();
    _state.resize( nrVars(), 0 );
}


void Gibbs::calc_factor_entries() {
    for( size_t I = 0; I < nrFactors(); I++ )
        _factor_entries[I] = get_factor_entry( I );
}

void Gibbs::update_factor_entries( size_t i ) {
    foreach( const Neighbor &I, nbV(i) )
        _factor_entries[I] = get_factor_entry( I );
}


void Gibbs::update_counts() {
    for( size_t i = 0; i < nrVars(); i++ )
        _var_counts[i][_state[i]]++;
    for( size_t I = 0; I < nrFactors(); I++ )
        _factor_counts[I][_factor_entries[I]]++;
//        _factor_counts[I][get_factor_entry(I)]++;
    _sample_count++;
}


inline size_t Gibbs::get_factor_entry( size_t I ) {
    size_t f_entry = 0;
    VarSet::const_reverse_iterator check = factor(I).vars().rbegin();
    for( int _j = nbF(I).size() - 1; _j >= 0; _j-- ) {
        size_t j = nbF(I)[_j];     // FIXME
        assert( var(j) == *check );
        f_entry *= var(j).states();
        f_entry += _state[j];
        check++;
    }
    return f_entry;
}


inline size_t Gibbs::get_factor_entry_interval( size_t I, size_t i ) {
    size_t skip = 1;
    VarSet::const_iterator check = factor(I).vars().begin();
    for( size_t _j = 0; _j < nbF(I).size(); _j++ ) {
        size_t j = nbF(I)[_j];     // FIXME
        assert( var(j) == *check );
        if( i == j )
            break;
        else
            skip *= var(j).states();
        check++;
    }
    return skip;
}


Prob Gibbs::get_var_dist( size_t i ) {
    assert( i < nrVars() );
    size_t i_states = var(i).states();
    Prob i_given_MB( i_states, 1.0 );

    // use markov blanket of var(i) to calculate distribution
    foreach( const Neighbor &I, nbV(i) ) {
        const Factor &f_I = factor(I);
        size_t I_skip = get_factor_entry_interval( I, i );
//        size_t I_entry = get_factor_entry(I) - (_state[i] * I_skip);
        size_t I_entry = _factor_entries[I] - (_state[i] * I_skip);
        for( size_t st_i = 0; st_i < i_states; st_i++ ) {
            i_given_MB[st_i] *= f_I[I_entry];
            I_entry += I_skip;
        }
    }

    return i_given_MB.normalized();
}


inline void Gibbs::resample_var( size_t i ) {
    // draw randomly from conditional distribution and update _state
    size_t new_state = get_var_dist(i).draw();
    if( new_state != _state[i] ) {
        _state[i] = new_state;
        update_factor_entries( i );
    }
}


void Gibbs::randomize_state() {
    for( size_t i = 0; i < nrVars(); i++ )
        _state[i] = rnd_int( 0, var(i).states() - 1 );
}


void Gibbs::init() {
    for( size_t i = 0; i < nrVars(); i++ )
        fill( _var_counts[i].begin(), _var_counts[i].end(), 0 );
    for( size_t I = 0; I < nrFactors(); I++ )
        fill( _factor_counts[I].begin(), _factor_counts[I].end(), 0 );
    _sample_count = 0;
}


double Gibbs::run() {
    if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";
    if( props.verbose >= 3 )
        cout << endl;

    double tic = toc();
    
    randomize_state();

    calc_factor_entries();
    for( size_t iter = 0; iter < props.iters; iter++ ) {
        for( size_t j = 0; j < nrVars(); j++ )
            resample_var( j );
        update_counts();
    }

    if( props.verbose >= 3 ) {
        for( size_t i = 0; i < nrVars(); i++ ) {
            cerr << "belief for variable " << var(i) << ": " << beliefV(i) << endl;
            cerr << "counts for variable " << var(i) << ": " << Prob( _var_counts[i].begin(), _var_counts[i].end() ) << endl;
        }
    }
    
    if( props.verbose >= 3 )
        cout << "Gibbs::run:  ran " << props.iters << " passes (" << toc() - tic << " clocks)." << endl;

    return 0.0;
}


inline Factor Gibbs::beliefV( size_t i ) const {
    return Factor( var(i), _var_counts[i].begin() ).normalized();
}


inline Factor Gibbs::beliefF( size_t I ) const {
    return Factor( factor(I).vars(), _factor_counts[I].begin() ).normalized();
}


Factor Gibbs::belief( const Var &n ) const {
    return( beliefV( findVar( n ) ) );
}


vector<Factor> Gibbs::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( beliefV(i) );
    for( size_t I = 0; I < nrFactors(); I++ )
        result.push_back( beliefF(I) );
    return result;
}


Factor Gibbs::belief( const VarSet &ns ) const {
    if( ns.size() == 1 )
        return belief( *(ns.begin()) );
    else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        assert( I != nrFactors() );
        return beliefF(I).marginal(ns);
    }
}


} // end of namespace dai
