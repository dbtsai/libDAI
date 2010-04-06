/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <dai/mf.h>
#include <dai/util.h>


namespace dai {


using namespace std;


const char *MF::Name = "MF";


void MF::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );
    DAI_ASSERT( opts.hasKey("maxiter") );

    props.tol = opts.getStringAs<Real>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0U;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<Real>("damping");
    else
        props.damping = 0.0;
}


PropertySet MF::getProperties() const {
    PropertySet opts;
    opts.set( "tol", props.tol );
    opts.set( "maxiter", props.maxiter );
    opts.set( "verbose", props.verbose );
    opts.set( "damping", props.damping );
    return opts;
}


string MF::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "damping=" << props.damping << "]";
    return s.str();
}


void MF::construct() {
    // create beliefs
    _beliefs.clear();
    _beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        _beliefs.push_back( Factor( var(i) ) );
}


string MF::identify() const {
    return string(Name) + printProperties();
}


void MF::init() {
    for( vector<Factor>::iterator qi = _beliefs.begin(); qi != _beliefs.end(); qi++ )
        qi->fill(1.0);
}


Factor MF::calcNewBelief( size_t i ) {
    Factor result;
    foreach( const Neighbor &I, nbV(i) ) {
        Factor henk;
        foreach( const Neighbor &j, nbF(I) ) // for all j in I \ i
            if( j != i )
                henk *= _beliefs[j];
        Factor piet = factor(I).log(true);
        piet *= henk;
        piet = piet.marginal(var(i), false);
        piet = piet.exp();
        result *= piet;
    }
    result.normalize();
    return result;
}


Real MF::run() {
    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";

    double tic = toc();

    vector<size_t> update_seq;
    update_seq.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        update_seq.push_back( i );

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    Real maxDiff = INFINITY;
    for( _iters = 0; _iters < props.maxiter && maxDiff > props.tol; _iters++ ) {
        random_shuffle( update_seq.begin(), update_seq.end() );

        maxDiff = -INFINITY;
        foreach( const size_t &i, update_seq ) {
            Factor nb = calcNewBelief( i );

            if( nb.hasNaNs() ) {
                cerr << Name << "::run():  ERROR: new belief of variable " << var(i) << " has NaNs!" << endl;
                return 1.0;
            }

            if( props.damping != 0.0 )
                nb = (nb^(1.0 - props.damping)) * (_beliefs[i]^props.damping);

            maxDiff = std::max( maxDiff, dist( nb, _beliefs[i], DISTLINF ) );
            _beliefs[i] = nb;
        }

        if( props.verbose >= 3 )
            cerr << Name << "::run:  maxdiff " << maxDiff << " after " << _iters+1 << " passes" << endl;
    }

    if( maxDiff > _maxdiff )
        _maxdiff = maxDiff;

    if( props.verbose >= 1 ) {
        if( maxDiff > props.tol ) {
            if( props.verbose == 1 )
                cerr << endl;
            cerr << Name << "::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << maxDiff << endl;
        } else {
            if( props.verbose >= 3 )
                cerr << Name << "::run:  ";
            cerr << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    return maxDiff;
}


Factor MF::beliefV( size_t i ) const {
    return _beliefs[i].normalized();
}


Factor MF::belief (const VarSet &ns) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 )
        return beliefV( findVar( *(ns.begin()) ) );
    else {
        DAI_THROW(BELIEF_NOT_AVAILABLE);
        return Factor();
    }
}


vector<Factor> MF::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( beliefV(i) );
    return result;
}


Real MF::logZ() const {
    Real s = 0.0;

    for( size_t i = 0; i < nrVars(); i++ )
        s -= beliefV(i).entropy();
    for( size_t I = 0; I < nrFactors(); I++ ) {
        Factor henk;
        foreach( const Neighbor &j, nbF(I) )  // for all j in I
            henk *= _beliefs[j];
        henk.normalize();
        Factor piet;
        piet = factor(I).log(true);
        piet *= henk;
        s -= piet.sum();
    }

    return -s;
}


void MF::init( const VarSet &ns ) {
    for( size_t i = 0; i < nrVars(); i++ ) {
        if( ns.contains(var(i) ) )
            _beliefs[i].fill( 1.0 );
    }
}


} // end of namespace dai
