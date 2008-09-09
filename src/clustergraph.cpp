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


#include <set>
#include <vector>
#include <iostream>
#include <dai/varset.h>
#include <dai/clustergraph.h>


namespace dai {


using namespace std;


ClusterGraph ClusterGraph::VarElim( const std::vector<Var> & ElimSeq ) const {
    const long verbose = 0;

    // Make a copy
    ClusterGraph _Cl(*this);

    ClusterGraph result;
    _Cl.eraseNonMaximal();
    
    // Do variable elimination
    for( vector<Var>::const_iterator n = ElimSeq.begin(); n != ElimSeq.end(); n++ ) {
        assert( _Cl.vars() && *n );

        if( verbose >= 1 )
            cout << "Cost of eliminating " << *n << ": " << _Cl.eliminationCost( *n ) << " new edges" << endl;
        
        result.insert( _Cl.Delta(*n) );

        if( verbose >= 1 )
            cout << "_Cl = " << _Cl << endl;

        if( verbose >= 1 )
            cout << "After inserting " << _Cl.delta(*n) << ", _Cl = ";
        _Cl.insert( _Cl.delta(*n) );
        if( verbose >= 1 )
            cout << _Cl << endl;

        if( verbose >= 1 )
            cout << "After erasing clusters that contain " << *n <<  ", _Cl = ";
        _Cl.eraseSubsuming( *n );
        if( verbose >= 1 )
            cout << _Cl << endl;

        if( verbose >= 1 )
            cout << "After erasing nonmaximal clusters, _Cl = ";
        _Cl.eraseNonMaximal();
        if( verbose >= 1 )
            cout << _Cl << endl;
    }

    return result;
}


ClusterGraph ClusterGraph::VarElim_MinFill() const {
    const long verbose = 0;

    // Make a copy
    ClusterGraph _Cl(*this);
    VarSet _vars( vars() );

    ClusterGraph result;
    _Cl.eraseNonMaximal();
    
    // Do variable elimination
    while( !_vars.empty() ) {
        if( verbose >= 1 )
            cout << "Var  Eliminiation cost" << endl;
        VarSet::const_iterator lowest = _vars.end();
        size_t lowest_cost = -1UL;
        for( VarSet::const_iterator n = _vars.begin(); n != _vars.end(); n++ ) {
            size_t cost = _Cl.eliminationCost( *n );
            if( verbose >= 1 )
                cout << *n << "  " << cost << endl;
            if( lowest == _vars.end() || lowest_cost > cost ) {
                lowest = n;
                lowest_cost = cost;
            }
        }
        Var n = *lowest;

        if( verbose >= 1 )
            cout << "Lowest: " << n << " (" << lowest_cost << ")" << endl;

        result.insert( _Cl.Delta(n) );

        _Cl.insert( _Cl.delta(n) );
        _Cl.eraseSubsuming( n );
        _Cl.eraseNonMaximal();
        _vars /= n;

    }

    return result;
}


} // end of namespace dai
