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


#ifndef __BIPGRAPH_H__
#define __BIPGRAPH_H__


#include <vector>
#include <algorithm>


/// A BipartiteGraph represents a graph with two types of nodes
template <class T1, class T2> class BipartiteGraph {
    public:
        typedef std::pair<size_t,size_t>    _edge_t;
        typedef std::vector<size_t>         _nb_t;
        typedef _nb_t::const_iterator       _nb_cit;
        

    private:
        /// Vertices of first type
        std::vector<T1>                     _V1;

        /// Vertices of second type
        std::vector<T2>                     _V2;
        
        /// Edges, which are pairs (v1,v2) with v1 in _V1 and v2 in _V2
        std::vector<_edge_t>                _E12;

        /// Conversion matrix that computes which index of _E12 corresponds to a vertex index pair (v1,v2)
        std::vector<std::vector<size_t> >   _E12ind;

        /// Neighbour indices of vertices of first type
        std::vector<_nb_t>                  _nb1;

        /// Neighbour indices of vertices of second type
        std::vector<_nb_t>                  _nb2;


    public:
        /// Default constructor
        BipartiteGraph<T1,T2> () {};

        /// Copy constructor
        BipartiteGraph<T1,T2> ( const BipartiteGraph<T1,T2> & x ) : _V1(x._V1), _V2(x._V2), _E12(x._E12), _E12ind(x._E12ind), _nb1(x._nb1), _nb2(x._nb2) {};

        /// Assignment operator
        BipartiteGraph<T1,T2> & operator=(const BipartiteGraph<T1,T2> & x) {
            if( this != &x ) {
                _V1 =       x._V1;
                _V2 =       x._V2;
                _E12 =      x._E12;
                _E12ind =   x._E12ind;
                _nb1 =      x._nb1;
                _nb2 =      x._nb2;
            }
            return *this;
        }
        
        /// Provides read access to node of first type
        const T1 & V1( size_t i1 ) const { return _V1[i1]; }
        /// Provides full access to node of first type
        T1 & V1( size_t i1 ) { return _V1[i1]; }
        /// Provides read access to all nodes of first type
        const std::vector<T1> & V1s() const { return _V1; }
        /// Provides full access to all nodes of first type
        std::vector<T1> & V1s() { return _V1; }

        /// Provides read access to node of second type
        const T2 & V2( size_t i2 ) const { return _V2[i2]; }
        /// Provides full access to node of second type
        T2 & V2( size_t i2 ) { return _V2[i2]; }
        /// Provides read access to all nodes of second type
        const std::vector<T2> & V2s() const { return _V2; }
        /// Provides full access to all nodes of second type
        std::vector<T2> & V2s() { return _V2; }

        /// Provides read access to edge
        const _edge_t & edge(size_t ind) const { return _E12[ind]; }
        /// Provides full access to edge
        _edge_t & edge(size_t ind) { return _E12[ind]; }
        /// Provides read access to all edges
        const std::vector<_edge_t> & edges() const { return _E12; }
        /// Provides full access to all edges
        std::vector<_edge_t> & edges() { return _E12; }
        /// Returns number of edges
        size_t nr_edges() const { return _E12.size(); }

        /// Provides read access to neighbours of node of first type
        const _nb_t & nb1( size_t i1 ) const { return _nb1[i1]; }
        /// Provides full access to neighbours of node of first type
        _nb_t & nb1( size_t i1 ) { return _nb1[i1]; }

        /// Provides read access to neighbours of node of second type
        const _nb_t & nb2( size_t i2 ) const { return _nb2[i2]; }
        /// Provides full access to neighbours of node of second type
        _nb_t & nb2( size_t i2 ) { return _nb2[i2]; }

        /// Converts the pair of indices (i1,i2) to the corresponding edge index
        size_t VV2E( const size_t i1, const size_t i2 ) const { return _E12ind[i1][i2]; }


        /// Regenerates internal structures (_E12ind, _nb1 and _nb2) based upon _V1, _V2 and _E12
        void Regenerate() {
            // Calculate _nb1 and _nb2

            // Start with empty vectors
            _nb1.clear();
            _nb1.resize(_V1.size());
            // Start with empty vectors
            _nb2.clear();
            _nb2.resize(_V2.size());
            // Each edge yields a neighbour pair
            for( std::vector<_edge_t>::const_iterator e = _E12.begin(); e != _E12.end(); e++ ) {
                _nb1[e->first].push_back(e->second);
                _nb2[e->second].push_back(e->first);
            }
            // Remove duplicates from _nb1
            for( size_t i1 = 0; i1 < _V1.size(); i1++ ) {
                _nb_t::iterator new_end = unique(_nb1[i1].begin(), _nb1[i1].end());
                _nb1[i1].erase( new_end, _nb1[i1].end() );
            }
            // Remove duplicates from _nb2
            for( size_t i2 = 0; i2 < _V2.size(); i2++ ) {
                _nb_t::iterator new_end = unique(_nb2[i2].begin(), _nb2[i2].end());
                _nb2[i2].erase( new_end, _nb2[i2].end() );
            }

            // Calculate _E12ind
            
            // Allocate data structures
            _E12ind.clear();
            std::vector<size_t> col(_V2.size(),-1);
            _E12ind.assign(_V1.size(), col);
            // Assign elements
            for( size_t k = 0; k < _E12.size(); k++ )
                _E12ind[_E12[k].first][_E12[k].second] = k;
        }
};


#endif
