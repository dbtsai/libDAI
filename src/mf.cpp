/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
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

    props.tol = opts.getStringAs<double>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0U;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<double>("damping");
    else
        props.damping = 0.0;
}


PropertySet MF::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "damping", props.damping );
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


double MF::run() {
    double tic = toc();

    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";

    size_t pass_size = _beliefs.size();
    Diffs diffs(pass_size * 3, 1.0);

    size_t t=0;
    for( t=0; t < (props.maxiter*pass_size) && diffs.maxDiff() > props.tol; t++ ) {
        // choose random Var i
        size_t i = (size_t) (nrVars() * rnd_uniform());

        Factor jan;
        Factor piet;
        foreach( const Neighbor &I, nbV(i) ) {
            Factor henk;
            foreach( const Neighbor &j, nbF(I) ) // for all j in I \ i
                if( j != i )
                    henk *= _beliefs[j];
            piet = factor(I).log(true);
            piet *= henk;
            piet = piet.marginal(var(i), false);
            piet = piet.exp();
            jan *= piet;
        }

        jan.normalize();

        if( jan.hasNaNs() ) {
            cerr << Name << "::run():  ERROR: jan has NaNs!" << endl;
            return 1.0;
        }

        if( props.damping != 0.0 )
            jan = (jan^(1.0 - props.damping)) * (_beliefs[i]^props.damping);
        diffs.push( dist( jan, _beliefs[i], Prob::DISTLINF ) );

        _beliefs[i] = jan;
    }

    _iters = t / pass_size;
    if( diffs.maxDiff() > _maxdiff )
        _maxdiff = diffs.maxDiff();

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cerr << endl;
            cerr << Name << "::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( props.verbose >= 2 )
                cerr << Name << "::run:  ";
            cerr << "converged in " << t / pass_size << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    return diffs.maxDiff();
}


Factor MF::beliefV( size_t i ) const {
    Factor piet;
    piet = _beliefs[i];
    piet.normalize();
    return(piet);
}


Factor MF::belief (const VarSet &ns) const {
    if( ns.size() == 1 )
        return belief( *(ns.begin()) );
    else {
        DAI_ASSERT( ns.size() == 1 );
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


Real MF::logZ() const {
    Real s = 0.0;

    for(size_t i=0; i < nrVars(); i++ )
        s -= beliefV(i).entropy();
    for(size_t I=0; I < nrFactors(); I++ ) {
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
