/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <dai/exactinf.h>
#include <sstream>


namespace dai {


using namespace std;


const char *ExactInf::Name = "EXACT";


void ExactInf::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("verbose") );

    props.verbose = opts.getStringAs<size_t>("verbose");
}


PropertySet ExactInf::getProperties() const {
    PropertySet opts;
    opts.Set( "verbose", props.verbose );
    return opts;
}


string ExactInf::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "verbose=" << props.verbose << "]";
    return s.str();
}


void ExactInf::construct() {
    // clear variable beliefs and reserve space
    _beliefsV.clear();
    _beliefsV.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefsV.push_back( Factor( var(i) ) );

    // clear factor beliefs and reserve space
    _beliefsF.clear();
    _beliefsF.reserve( nrFactors() );
    for( size_t I = 0; I < nrFactors(); I++ )
        _beliefsF.push_back( Factor( factor(I).vars() ) );
}


void ExactInf::init() {
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefsV[i].fill( 1.0 );
    for( size_t I = 0; I < nrFactors(); I++ )
        _beliefsF[I].fill( 1.0 );
}


Real ExactInf::run() {
    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";

    Factor P;
    for( size_t I = 0; I < nrFactors(); I++ )
        P *= factor(I);

    Real Z = P.sum();
    _logZ = std::log(Z);
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefsV[i] = P.marginal(var(i));
    for( size_t I = 0; I < nrFactors(); I++ )
        _beliefsF[I] = P.marginal(factor(I).vars());

    if( props.verbose >= 1 )
        cerr << "finished" << endl;

    return 0.0;
}


Factor ExactInf::calcMarginal( const VarSet &vs ) const {
    Factor P;
    for( size_t I = 0; I < nrFactors(); I++ )
        P *= factor(I);
    return P.marginal( vs, true );
}


vector<Factor> ExactInf::beliefs() const {
    vector<Factor> result = _beliefsV;
    result.insert( result.end(), _beliefsF.begin(), _beliefsF.end() );
    return result;
}


Factor ExactInf::belief( const VarSet &ns ) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 ) {
        return beliefV( findVar( *(ns.begin()) ) );
    } else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        if( I == nrFactors() )
            DAI_THROW(BELIEF_NOT_AVAILABLE);
        return beliefF(I).marginal(ns);
    }
}


string ExactInf::identify() const {
    return string(Name) + printProperties();
}


} // end of namespace dai
