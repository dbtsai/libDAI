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


#ifndef __defined_libdai_bipgraph_h
#define __defined_libdai_bipgraph_h


#include <vector>
#include <algorithm>
#include <dai/util.h>


namespace dai {


/// A BipartiteGraph represents a bipartite graph, with two types of nodes (both are numbered
/// as 0,1,2,...), with edges only between nodes of different type. The edges are stored as 
/// lists of adjacent nodes for each node.
class BipartiteGraph {
    public:
        /// A Neighbor describes a neighboring node of some other node.
        /** Iterating over all neighbors of some node i can be done in the following way:
         *  \code
         *      foreach( const BipartiteGraph::Neighbor &I, nb1(i) ) {
         *          size_t _I = I.iter;
         *          size_t _i = I.dual;
         *          // I == I.node;
         *          // The _I'th neighbor of i is I, and the _i'th neighbor of I is i:
         *          // nb1(i)[_I] == I, nb2(I)[_i] == i
         *      }
         *  \endcode
         */
        struct Neighbor {
            /// iter corresponds to the index of this Neighbor entry in the list of neighbors
            unsigned iter;
            /// node contains the number of the neighboring node
            unsigned node;
            /// dual contains the "dual" iter
            unsigned dual;
            /// cast to unsigned returns node
            operator unsigned () const { return node; };
        };

        /// Neighbors is a vector of Neigbor entries
        typedef std::vector<Neighbor> Neighbors;

    private:
        /// _nb1 contains for each node of the first kind a list of its neighbors
        std::vector<Neighbors> _nb1;
        /// _nb2 contains for each node of the second kind a list of its neighbors
        std::vector<Neighbors> _nb2;

    public:
        /// Default constructor
        BipartiteGraph() : _nb1(), _nb2() {}

        /// Copy constructor
        BipartiteGraph( const BipartiteGraph & x ) : _nb1(x._nb1), _nb2(x._nb2) {}

        /// Assignment operator
        BipartiteGraph & operator=( const BipartiteGraph & x ) {
            if( this != &x ) {
                _nb1 = x._nb1;
                _nb2 = x._nb2;
            }
            return *this;
        }

        /// Create bipartite graph from a range of edges, encoded as pairs of node numbers
        /// (more precisely, a std::pair<unsigned, unsigned> where the first integer corresponds
        /// to the node of the first type, and the second integer corresponds to the node of the
        /// second type). nr1 is the number of nodes of the first type, nr2 the number of nodes
        /// of the second type.
        template<typename EdgeInputIterator>
        void create( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end );

        /// Construct bipartite graph from a range of edges, encoded as pairs of node numbers
        /// (more precisely, a std::pair<unsigned, unsigned> where the first integer corresponds
        /// to the node of the first type, and the second integer corresponds to the node of the
        /// second type). nr1 is the number of nodes of the first type, nr2 the number of nodes
        /// of the second type.
        template<typename EdgeInputIterator>
        BipartiteGraph( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) : _nb1( nr1 ), _nb2( nr2 ) {
            create( nr1, nr2, begin, end );
        }

        /// Returns constant reference to the _i2'th neighbor of node i1 of first type
        const Neighbor & nb1( size_t i1, size_t _i2 ) const { return _nb1[i1][_i2]; }
        /// Returns reference to the _i2'th neighbor of node i1 of first type
        Neighbor & nb1( size_t i1, size_t _i2 ) { return _nb1[i1][_i2]; }

        /// Returns constant reference to the _i1'th neighbor of node i2 of second type
        const Neighbor & nb2( size_t i2, size_t _i1 ) const { return _nb2[i2][_i1]; }
        /// Returns reference to the _i1'th neighbor of node i2 of second type
        Neighbor & nb2( size_t i2, size_t _i1 ) { return _nb2[i2][_i1]; }

        /// Returns constant reference to all neighbors of node of first type
        const Neighbors & nb1( size_t i1 ) const { return _nb1[i1]; }
        /// Returns reference to all neighbors of node of first type
        Neighbors & nb1( size_t i1 ) { return _nb1[i1]; }

        /// Returns constant reference to all neighbors of node of second type
        const Neighbors & nb2( size_t i2 ) const { return _nb2[i2]; }
        /// Returns reference to all neighbors of node of second type
        Neighbors & nb2( size_t i2 ) { return _nb2[i2]; }

        /// Returns number of nodes of first type
        size_t nr1() const { return _nb1.size(); }
        /// Returns number of nodes of second type
        size_t nr2() const { return _nb2.size(); }
        
        /// Calculates the number of edges
        size_t nrEdges() const {
            size_t sum = 0;
            for( size_t i1 = 0; i1 < nr1(); i1++ )
                sum += nb1(i1).size();
            return sum;
        }

        /// Returns true if the graph is connected
        /// FIXME: this should be optimized
        bool isConnected() const {
            if( nr1() == 0 ) {
                return true;
            } else {
                std::vector<bool> incomponent1( nr1(), false );
                std::vector<bool> incomponent2( nr2(), false );

                incomponent1[0] = true;
                bool found_new_nodes;
                do {
                    found_new_nodes = false;

                    // For all nodes of second type, check if they are connected with the (growing) component
                    for( size_t n2 = 0; n2 < nr2(); n2++ )
                        if( !incomponent2[n2] ) {
                            foreach( const Neighbor &n1, nb2(n2) ) {
                                if( incomponent1[n1] ) {
                                    found_new_nodes = true;
                                    incomponent2[n2] = true;
                                    break;
                                }
                            }
                        }

                    // For all nodes of first type, check if they are connected with the (growing) component
                    for( size_t n1 = 0; n1 < nr1(); n1++ )
                        if( !incomponent1[n1] ) {
                            foreach( const Neighbor &n2, nb1(n1) ) {
                                if( incomponent2[n2] ) {
                                    found_new_nodes = true;
                                    incomponent1[n1] = true;
                                    break;
                                }
                            }
                        }
                } while( found_new_nodes );

                // Check if there are remaining nodes (not in the component)
                bool all_connected = true;
                for( size_t n1 = 0; (n1 < nr1()) && all_connected; n1++ )
                    if( !incomponent1[n1] )
                        all_connected = false;
                for( size_t n2 = 0; (n2 < nr2()) && all_connected; n2++ )
                    if( !incomponent2[n2] )
                        all_connected = false;

                return all_connected;
            }
        }

};


template<typename EdgeInputIterator>
void BipartiteGraph::create( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) {
    _nb1.clear();
    _nb1.resize( nr1 );
    _nb2.clear();
    _nb2.resize( nr2 );

    for( EdgeInputIterator e = begin; e != end; e++ ) {
        // Each edge yields a neighbor pair
        Neighbor nb_1;
        nb_1.iter = _nb1[e->first].size();
        nb_1.node = e->second;
        nb_1.dual = _nb2[e->second].size();

        Neighbor nb_2;
        nb_2.iter = nb_1.dual;
        nb_2.node = e->first;
        nb_2.dual = nb_1.iter;

        _nb1[e->first].push_back( nb_1 );
        _nb2[e->second].push_back( nb_2 );
    }
}


} // end of namespace dai


#endif
