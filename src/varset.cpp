/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Copyright (C) 2002  Martijn Leisink  [martijn@mbfys.kun.nl]
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


/*#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <dai/var.h>
#include <dai/util.h>
*/
#include <dai/varset.h>


namespace dai {


using namespace std;


/// Calculates the product of number of states of all variables in vars
size_t nrStates( const VarSet &vars ) {
    size_t states = 1;
    for( VarSet::const_iterator n = vars.begin(); n != vars.end(); n++ )
        states *= n->states();
    return states;
}


/// calcState calculates the linear index of vars that corresponds
/// to the states of the variables given in states, implicitly assuming
/// states[m] = 0 for all m in this VarSet which are not in states.
size_t calcState( const VarSet &vars, const std::map<Var, size_t> &states ) {
    size_t prod = 1;
    size_t state = 0;
    for( VarSet::const_iterator n = vars.begin(); n != vars.end(); n++ ) {
        map<Var, size_t>::const_iterator m = states.find( *n );
        if( m != states.end() )
            state += prod * m->second;
        prod *= n->states();
    }
    return state;
}


/// Sends a VarSet to an output stream
std::ostream& operator<< (std::ostream &os, const VarSet& ns) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
        os << *n;
    return( os );
}


} // end of namespace dai
