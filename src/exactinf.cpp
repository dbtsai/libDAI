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


#include <dai/exactinf.h>
#include <sstream>


namespace dai {


using namespace std;


const char *ExactInf::Name = "EXACT";


void ExactInf::setProperties( const PropertySet &opts ) {
    assert( opts.hasKey("verbose") );
    
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


double ExactInf::run() {
    if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";

    Factor P;
    for( size_t I = 0; I < nrFactors(); I++ )
        P *= factor(I);

    Real Z = P.totalSum();
    _logZ = std::log(Z);
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefsV[i] = P.marginal(var(i));
    for( size_t I = 0; I < nrFactors(); I++ )
        _beliefsF[I] = P.marginal(factor(I).vars());

    if( props.verbose >= 1 )
        cout << "finished" << endl;

    return 0.0;
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
        return belief( *(ns.begin()) );
    } else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        assert( I != nrFactors() );
        return beliefF(I).marginal(ns);
    }
}


string ExactInf::identify() const { 
    return string(Name) + printProperties();
}


} // end of namespace dai
