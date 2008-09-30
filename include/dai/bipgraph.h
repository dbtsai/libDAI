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


/// \file
/// \brief Defines BipartiteGraph class


#ifndef __defined_libdai_bipgraph_h
#define __defined_libdai_bipgraph_h


#include <ostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <dai/util.h>


namespace dai {


/// Represents the neighborhood structure of nodes in a bipartite graph.
/** A bipartite graph has two types of nodes: type 1 and type 2. Edges can occur only between
 *  nodes of different type. Nodes are indexed by an unsigned integer, edges are indexed as
 *  a pair of unsigned integers, where the pair (a,b) means the b'th neighbor of the a'th node.
 *
 *  The BipartiteGraph stores for each node of type 1 a vector of Neighbor structures, where
 *  each Neighbor corresponds to a neighboring node of type 2. In addition, each node of type 2
 *  stores a vector of Neighbor structures describing its neighboring nodes of type 1.
 */
class BipartiteGraph {
    public:
        /// Describes a neighboring node of some other node in a BipartiteGraph.
        /** Iterating over all neighbors of the n1'th node of type 1 can be done in the following way:
         *  \code
         *      size_t n1 = ...;
         *      foreach( const BipartiteGraph::Neighbor &n2, nb1(n1) ) {
         *          size_t _n2 = n2.iter;
         *          size_t _n1 = n2.dual;
         *          std::cout << "The " << _n2 << "'th neighbor of the " << n1 << "'th node of type 1 is: the " << n2 << "'th node of type 2" << endl;
         *
         *          // The _n2'th neighbor of n1 is n2:
         *          assert( nb1(n1)[_n2] == n2 );
         *          // The _n1'th neighbor of n2 is n1:
         *          assert( nb2(n2)[_n1] == n1 );
         *          // n2 can be used as an abbreviation of n2.node:
         *          assert( static_cast<size_t>(n2) == n2.node );
         *      }
         *  \endcode
         */
        struct Neighbor {
            /// Corresponds to the index of this Neighbor entry in the vector of neighbors
            size_t iter;
            /// Contains the number of the neighboring node
            size_t node;
            /// Contains the "dual" iter
            size_t dual;

            /// Default constructor
            Neighbor() {}
            /// Constructor that sets the Neighbor members accordingly to the parameters
            Neighbor( size_t iter, size_t node, size_t dual ) : iter(iter), node(node), dual(dual) {}

            /// Cast to size_t returns member node
            operator size_t () const { return node; }
        };

        /// Describes the neighbors of some node.
        typedef std::vector<Neighbor> Neighbors;

        /// Used as index of an edge: an Edge(a,b) corresponds to the edge between the a'th node and its b'th neighbor (it depends on the context whether the first node (with index a) is of type 1 or of type 2).
        typedef std::pair<size_t,size_t> Edge;

    private:
        /// Contains for each node of type 1 a vector of its neighbors
        std::vector<Neighbors> _nb1;

        /// Contains for each node of type 2 a vector of its neighbors
        std::vector<Neighbors> _nb2;

        /// Used internally by isTree()
        struct levelType {
            std::vector<size_t> ind1;       // indices of nodes of type 1
            std::vector<size_t> ind2;       // indices of nodes of type 2
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

        /// Constructs BipartiteGraph from a range of edges. 
        /** \tparam EdgeInputIterator Iterator with value_type Edge.
         *  \param nr1 The number of nodes of type 1.
         *  \param nr2 The number of nodes of type 2. 
         *  \param begin Points to the first Edge.
         *  \param end Points just beyond the last Edge.
         */
        template<typename EdgeInputIterator>
        BipartiteGraph( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) : _nb1( nr1 ), _nb2( nr2 ) {
            construct( nr1, nr2, begin, end );
        }

        /// (Re)constructs BipartiteGraph from a range of edges. 
        /** \tparam EdgeInputIterator Iterator with value_type Edge.
         *  \param nr1 The number of nodes of type 1.
         *  \param nr2 The number of nodes of type 2. 
         *  \param begin Points to the first Edge.
         *  \param end Points just beyond the last Edge.
         */
        template<typename EdgeInputIterator>
        void construct( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end );

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
        
        /// Adds a node of type 1 without neighbors.
        void add1() { _nb1.push_back( Neighbors() ); }
        
        /// Adds a node of type 2 without neighbors.
        void add2() { _nb2.push_back( Neighbors() ); }

        /// Adds a node of type 1, with neighbors specified by a range of indices of nodes of type 2.
        /** \tparam NodeInputIterator Iterator with value_type size_t, corresponding to
         *  the indices of nodes of type 2 that should become neighbors of the added node.
         *  \param begin Points to the index of the first neighbor.
         *  \param end Points just beyond the index of the last neighbor.
         *  \param sizeHint For improved efficiency, the size of the range may be specified by sizeHint.
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

        /// Adds a node of type 2, with neighbors specified by a range of indices of nodes of type 1.
        /** \tparam NodeInputIterator Iterator with value_type size_t, corresponding to
         *  the indices of nodes of type 1 that should become neighbors of the added node.
         *  \param begin Points to the index of the first neighbor.
         *  \param end Points just beyond the index of the last neighbor.
         *  \param sizeHint For improved efficiency, the size of the range may be specified by sizeHint.
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

        /// Removes node n1 of type 1 and all incident edges.
        void erase1( size_t n1 );

        /// Removes node n2 of type 2 and all incident edges.
        void erase2( size_t n2 );

        /// Adds an edge between node n1 of type 1 and node n2 of type 2.
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

        /// Calculates second-order neighbors (i.e., neighbors of neighbors) of node n1 of type 1.
        /** If include == true, includes n1 itself, otherwise excludes n1.
         */
        std::vector<size_t> delta1( size_t n1, bool include = false ) const;

        /// Calculates second-order neighbors (i.e., neighbors of neighbors) of node n2 of type 2.
        /** If include == true, includes n2 itself, otherwise excludes n2.
         */
        std::vector<size_t> delta2( size_t n2, bool include = false ) const;

        /// Returns true if the graph is connected
        bool isConnected() const;

        /// Returns true if the graph is a tree, i.e., if it is singly connected and connected.
        /** This is equivalent to whether for each pair of nodes in the graph, there exists
         *  a unique path in the graph that starts at the first and ends at the second node.
         */
        bool isTree() const;

        /// Writes this BipartiteGraph to an output stream in GraphViz .dot syntax
        void printDot( std::ostream& os ) const;

    private:
        /// Checks internal consistency
        void check() const;
};


template<typename EdgeInputIterator>
void BipartiteGraph::construct( size_t nr1, size_t nr2, EdgeInputIterator begin, EdgeInputIterator end ) {
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
