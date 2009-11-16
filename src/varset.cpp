/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij      [joris dot mooij at libdai dot org]
 *  Copyright (C) 2002-2007  Radboud University Nijmegen, The Netherlands
 */


#include <dai/varset.h>


namespace dai {


using namespace std;


size_t calcLinearState( const VarSet &vs, const std::map<Var, size_t> &state ) {
    size_t prod = 1;
    size_t st = 0;
    for( VarSet::const_iterator v = vs.begin(); v != vs.end(); v++ ) {
        std::map<Var, size_t>::const_iterator m = state.find( *v );
        if( m != state.end() )
            st += prod * m->second;
        prod *= v->states();
    }
    return st;
}


std::map<Var, size_t> calcState( const VarSet &vs, size_t linearState ) {
    std::map<Var, size_t> state;
    for( VarSet::const_iterator v = vs.begin(); v != vs.end(); v++ ) {
        state[*v] = linearState % v->states();
        linearState /= v->states();
    }
    DAI_ASSERT( linearState == 0 );
    return state;
}


} // end of namespace dai
