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


#include <ostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <dai/util.h>


namespace dai {


/// A BipartiteGraph represents the neighborhood structure of nodes in a bipartite graph.
/** A bipartite graph has two types of nodes: type 1 and type 2. Edges can occur only between
 *  nodes of different type. Nodes are indexed by an unsigned integer, edges are indexed as
 *  a pair of unsigned integers, where the pair (a,b) means the b'th neighbor of the a'th node.
 *  The BipartiteGraph stores neighborhood structures as vectors of vectors of Neighbor entries:
 *  each node has a vector of Neighbor entries, which is also known as a Neighbors type.
 */
class BipartiteGraph {
    public:
        /// A Neighbor describes a neighboring node of some other node in a BipartiteGraph.
        /** Iterating over all neighbors of node n1 of type 1 can be done in the following way:
         *  \code
         *      foreach( const BipartiteGraph::Neighbor &n2, nb1(n1) ) {
         *          size_t _n2 = n2.iter;
         *          size_t _n1 = n2.dual;
         *          // n2 == n2.node;
         *          // The _n2'th neighbor of n1 is n2, and the _n1'th neighbor of n2 is n1:
         *          // nb1(n1)[_n2] == n2, nb2(n2)[_n1] == n1 
         *      }
         *  \endcode
         */
        struct Neighbor {
            /// Corresponds to the index of this Neighbor entry in the vector of neighbors
            unsigned iter;
            /// Contains the number of the neighboring node
            unsigned node;
            /// Contains the "dual" iter
            unsigned dual;
            /// Cast to unsigned returns node member
            operator unsigned () const { return node; }
            /// Default constructor
            Neighbor() {}
            /// Constructor
            Neighbor( size_t iter, size_t node, size_t dual ) : iter(iter), node(node), dual(dual) {}
        };

        /// Neighbors is a vector of Neighbor entries; each node has an associated Neighbors variable, which describes its neighbors.
        typedef std::vector<Neighbor> Neighbors;

        /// Edge is used as index of an edge: an Edge(a,b) corresponds to the edge between the a'th node and its b'th neighbor. It depends on the context whether the first node (with index a) is of type 1 or of type 2.
        typedef std::pair<size_t,size_t> Edge;

    private:
        /// _nb1 contains for each node of type 1 a vector of its neighbors
        std::vector<Neighbors> _nb1;
        /// _nb2 contains for each node of type 2 a vector of its neighbors
        std::vector<Neighbors> _nb2;

        /// Used internally by isTree()
        struct levelType {
            std::vector<size_t> ind1;       // indices of vertices of type 1
            std::vector<size_t> ind2;       // indices of vertices of type 2
        };

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

        /// Create bipartite graph from a range of edges. 
        /** nr1 is the number of nodes of type 1, nr2 the number of nodes of type 2. 
         *  The value_type of an EdgeInputIterator should be Edge.
         */
        template<typename EdgeInputIterator>
        void create( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end );

        /// Construct bipartite graph from a range of edges. 
        /** nr1 is the number of nodes of type 1, nr2 the number of nodes of type 2. 
         *  The value_type of an EdgeInputIterator should be Edge.
         */
        template<typename EdgeInputIterator>
        BipartiteGraph( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) : _nb1( nr1 ), _nb2( nr2 ) {
            create( nr1, nr2, begin, end );
        }

        /// Returns constant reference to the _i2'th neighbor of node i1 of type 1
        const Neighbor & nb1( size_t i1, size_t _i2 ) const { 
#ifdef DAI_DEBUG
            assert( i1 < _nb1.size() );
            assert( _i2 < _nb1[i1].size() );
#endif
            return _nb1[i1][_i2]; 
        }
        /// Returns reference to the _i2'th neighbor of node i1 of type 1
        Neighbor & nb1( size_t i1, size_t _i2 ) {
#ifdef DAI_DEBUG
            assert( i1 < _nb1.size() );
            assert( _i2 < _nb1[i1].size() );
#endif
            return _nb1[i1][_i2]; 
        }

        /// Returns constant reference to the _i1'th neighbor of node i2 of type 2
        const Neighbor & nb2( size_t i2, size_t _i1 ) const { 
#ifdef DAI_DEBUG
            assert( i2 < _nb2.size() );
            assert( _i1 < _nb2[i2].size() );
#endif
            return _nb2[i2][_i1]; 
        }
        /// Returns reference to the _i1'th neighbor of node i2 of type 2
        Neighbor & nb2( size_t i2, size_t _i1 ) { 
#ifdef DAI_DEBUG
            assert( i2 < _nb2.size() );
            assert( _i1 < _nb2[i2].size() );
#endif
            return _nb2[i2][_i1]; 
        }

        /// Returns constant reference to all neighbors of node i1 of type 1
        const Neighbors & nb1( size_t i1 ) const { 
#ifdef DAI_DEBUG
            assert( i1 < _nb1.size() );
#endif
            return _nb1[i1]; 
        }
        /// Returns reference to all neighbors of node of i1 type 1
        Neighbors & nb1( size_t i1 ) { 
#ifdef DAI_DEBUG
            assert( i1 < _nb1.size() );
#endif
            return _nb1[i1]; 
        }

        /// Returns constant reference to all neighbors of node i2 of type 2
        const Neighbors & nb2( size_t i2 ) const { 
#ifdef DAI_DEBUG
            assert( i2 < _nb2.size() );
#endif
            return _nb2[i2]; 
        }
        /// Returns reference to all neighbors of node i2 of type 2
        Neighbors & nb2( size_t i2 ) { 
#ifdef DAI_DEBUG
            assert( i2 < _nb2.size() );
#endif
            return _nb2[i2]; 
        }

        /// Returns number of nodes of type 1
        size_t nr1() const { return _nb1.size(); }
        /// Returns number of nodes of type 2
        size_t nr2() const { return _nb2.size(); }
        
        /// Calculates the number of edges, using O(nr1()) time
        size_t nrEdges() const {
            size_t sum = 0;
            for( size_t i1 = 0; i1 < nr1(); i1++ )
                sum += nb1(i1).size();
            return sum;
        }
        
        /// Add node of type 1 without neighbors.
        void add1() { _nb1.push_back( Neighbors() ); }
        
        /// Add node of type 2 without neighbors.
        void add2() { _nb2.push_back( Neighbors() ); }

        /// Add node of type 1 with neighbors specified by a range.
        /** The value_type of an NodeInputIterator should be a size_t, corresponding to
         *  the indices of nodes of type 2 that should become neighbors of the added node.
         *  For improved efficiency, the size of the range may be specified by sizeHint.
         */
        template <typename NodeInputIterator>
        void add1( NodeInputIterator begin, NodeInputIterator end, size_t sizeHint = 0 ) {
            Neighbors nbs1new;
            nbs1new.reserve( sizeHint );
            size_t iter = 0;
            for( NodeInputIterator it = begin; it != end; ++it ) {
                assert( *it < nr2() );
                Neighbor nb1new( iter, *it, nb2(*it).size() );
                Neighbor nb2new( nb2(*it).size(), nr1(), iter++ );
                nbs1new.push_back( nb1new );
                nb2( *it ).push_back( nb2new );
            }
            _nb1.push_back( nbs1new );
        }

        /// Add node of type 2 with neighbors specified by a range.
        /** The value_type of an NodeInputIterator should be a size_t, corresponding to
         *  the indices of nodes of type 1 that should become neighbors of the added node.
         *  For improved efficiency, the size of the range may be specified by sizeHint.
         */
        template <typename NodeInputIterator>
        void add2( NodeInputIterator begin, NodeInputIterator end, size_t sizeHint = 0 ) {
            Neighbors nbs2new;
            nbs2new.reserve( sizeHint );
            size_t iter = 0;
            for( NodeInputIterator it = begin; it != end; ++it ) {
                assert( *it < nr1() );
                Neighbor nb2new( iter, *it, nb1(*it).size() );
                Neighbor nb1new( nb1(*it).size(), nr2(), iter++ );
                nbs2new.push_back( nb2new );
                nb1( *it ).push_back( nb1new );
            }
            _nb2.push_back( nbs2new );
        }

        /// Remove node n1 of type 1 and all incident edges.
        void erase1( size_t n1 );

        /// Remove node n2 of type 2 and all incident edges.
        void erase2( size_t n2 );

        /// Add edge between node n1 of type 1 and node n2 of type 2.
        /** If check == true, only adds the edge if it does not exist already.
         */
        void addEdge( size_t n1, size_t n2, bool check = true ) {
            assert( n1 < nr1() );
            assert( n2 < nr2() );
            bool exists = false;
            if( check ) {
                // Check whether the edge already exists
                foreach( const Neighbor &nb2, nb1(n1) )
                    if( nb2 == n2 ) {
                        exists = true;
                        break;
                    }
            }
            if( !exists ) { // Add edge
                Neighbor nb_1( _nb1[n1].size(), n2, _nb2[n2].size() );
                Neighbor nb_2( nb_1.dual, n1, nb_1.iter );
                _nb1[n1].push_back( nb_1 );
                _nb2[n2].push_back( nb_2 );
            }
        }

        /// Calculate second-order neighbors (i.e., neighbors of neighbors) of node n1 of type 1.
        /** If include == true, include n1 itself, otherwise exclude n1.
         */
        std::vector<size_t> delta1( size_t n1, bool include = false ) const;

        /// Calculate second-order neighbors (i.e., neighbors of neighbors) of node n2 of type 2.
        /** If include == true, include n2 itself, otherwise exclude n2.
         */
        std::vector<size_t> delta2( size_t n2, bool include = false ) const;

        /// Returns true if the graph is connected
        bool isConnected() const;

        /// Returns true if the graph is a tree, i.e., if it is singly connected and connected.
        /** This is equivalent to whether for each pair of vertices in the graph, there exists
         *  a unique path in the graph that starts at the first and ends at the second vertex.
         */
        bool isTree() const;

        /// Stream to output stream os in graphviz .dot syntax
        void display( std::ostream& os ) const;

        /// Checks internal consistency
        void check() const;
};


template<typename EdgeInputIterator>
void BipartiteGraph::create( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) {
    _nb1.clear();
    _nb1.resize( nr1 );
    _nb2.clear();
    _nb2.resize( nr2 );

    for( EdgeInputIterator e = begin; e != end; e++ ) {
#ifdef DAI_DEBUG
        addEdge( e->first, e->second, true );
#else
        addEdge( e->first, e->second, false );
#endif
    }
}


} // end of namespace dai


#endif
