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
}


void Gibbs::update_counts( _state_t &st ) {
    for( size_t i = 0; i < nrVars(); i++ )
        _var_counts[i][st[i]]++;
    for( size_t I = 0; I < nrFactors(); I++ ) {
        if( 0 ) {
/*            multind mi( factor(I).vars() );
            _state_t f_st( factor(I).vars().size() );
            int k = 0;
            foreach( size_t j, nbF(I) )
                f_st[k++] = st[j];
            _factor_counts[I][mi.li(f_st)]++;*/
        } else {
            size_t ent = get_factor_entry(st, I);
            _factor_counts[I][ent]++;
        }
    }
    _sample_count++;
}


inline
size_t Gibbs::get_factor_entry(const _state_t &st, int factor) {
  size_t f_entry=0;
  int rank = nbF(factor).size();
  for(int j=rank-1; j>=0; j--) {
      int jn = nbF(factor)[j];
      f_entry *= var(jn).states();
      f_entry += st[jn];
  }
  return f_entry;
}


Prob Gibbs::get_var_dist( _state_t &st, size_t i ) {
    assert( st.size() == vars().size() );
    assert( i < nrVars() );
    if( 1 ) {
        // use markov blanket of n to calculate distribution
        size_t dim = var(i).states();
        Neighbors &facts = nbV(i);

        Prob values( dim, 1.0 );

        for( size_t I = 0; I < facts.size(); I++ ) {
            size_t fa = facts[I];
            const Factor &f = factor(fa);
            int save_ind = st[i];
            for( size_t k = 0; k < dim; k++ ) {
                st[i] = k;
                int f_entry = get_factor_entry(st, fa);
                values[k] *= f[f_entry];
            }
            st[i] = save_ind;
        }

        return values.normalized();
    } else {
/*        Var vi = var(i);
        Factor d(vi);
        assert(vi.states()>0);
        assert(vi.label()>=0);
        // loop over factors containing i (nbV(i)):
        foreach(size_t I, nbV(i)) {
            // use multind to find linear state for variables != i in factor
            assert(I<nrFactors());
            assert(factor(I).vars().size() > 0);
            VarSet vs (factor(I).vars() / vi);
            multind mi(vs);
            _state_t I_st(vs.size());
            int k=0;
            foreach(size_t l, nbF(I)) {
                if(l!=i) I_st[k++] = st[l];
            }
            // use slice(ns,ns_state) to get beliefs for variable i
            // multiply all these beliefs together
            d *= factor(I).slice(vs, mi.li(I_st));
        }
        d.p().normalize();
        return d.p();*/
    }
}


void Gibbs::resample_var( _state_t &st, size_t i ) {
    // draw randomly from conditional distribution and update 'st'
    st[i] = get_var_dist( st, i ).draw();
}


void Gibbs::randomize_state( _state_t &st ) {
    assert( st.size() == nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        st[i] = rnd_int( 0, var(i).states() - 1 );
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
    
    vector<size_t> state( nrVars() );
    randomize_state( state );

    for( size_t iter = 0; iter < props.iters; iter++ ) {
        for( size_t j = 0; j < nrVars(); j++ )
            resample_var( state, j );
        update_counts( state );
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


Factor Gibbs::beliefV( size_t i ) const {
    Prob p( _var_counts[i].begin(), _var_counts[i].end() );
    p.normalize();
    return( Factor( var(i), p ) );
}


Factor Gibbs::beliefF( size_t I ) const {
    Prob p( _factor_counts[I].begin(), _factor_counts[I].end() );
    p.normalize();
    return( Factor( factor(I).vars(), p ) );
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
