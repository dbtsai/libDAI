/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines the GraphAL class, which represents an undirected graph as an adjacency list


#ifndef __defined_libdai_graph_h
#define __defined_libdai_graph_h


#include <ostream>
#include <vector>
#include <algorithm>
#include <dai/util.h>
#include <dai/exceptions.h>


namespace dai {


/// Represents the neighborhood structure of nodes in an undirected graph.
/** A graph has nodes connected by edges. Nodes are indexed by an unsigned integer. 
 *  If there are nrNodes() nodes, they are numbered 0,1,2,...,nrNodes()-1. An edge
 *  between node \a n1 and node \a n2 is represented by a GraphAL::Edge(\a n1,\a n2).
 *
 *  GraphAL is implemented as a sparse adjacency list, i.e., it stores for each node a list of
 *  its neighboring nodes. The list of neighboring nodes is implemented as a vector of Neighbor
 *  structures (accessible by the nb() method). Thus, each node has an associated variable of 
 *  type GraphAL::Neighbors, which is a vector of Neighbor structures, describing its 
 *  neighboring nodes.
 */
class GraphAL {
    public:
        /// Describes the neighbor relationship of two nodes in a GraphAL.
        /** Sometimes we want to do an action, such as sending a
         *  message, for all edges in a graph. However, most graphs
         *  will be sparse, so we need some way of storing a set of
         *  the neighbors of a node, which is both fast and
         *  memory-efficient. We also need to be able to go between
         *  viewing node \a a as a neighbor of node \a b, and node \a b
         *  as a neighbor of node \a a. The Neighbor struct solves
         *  both of these problems. Each node has a list of neighbors,
         *  stored as a std::vector<\link Neighbor \endlink>, and 
         *  extra information is included in the Neighbor struct which 
         *  allows us to access a node as a neighbor of its neighbor 
         *  (the \c dual member).
         *
         *  By convention, variable identifiers naming indices into a
         *  vector of neighbors are prefixed with an underscore ("_").
         *  The neighbor list which they point into is then understood
         *  from context. For example:
         *
         *  \code
         *  Real MR::T( size_t i, size_t _j );
         *  \endcode
         *
         *  Here, \a i is the "absolute" index of node i, but \a _j is
         *  understood as a "relative" index, giving node \a j 's entry in
         *  <tt>nb(i)</tt>. The corresponding Neighbor structure can be
         *  accessed as <tt>nb(i,_j)</tt> or <tt>nb(i)[_j]</tt>. The 
         *  absolute index of \a _j, which would be called \a j, can be 
         *  recovered from the \c node member: <tt>nb(i,_j).node</tt>. 
         *  The \c iter member gives the relative index \a _j, and the 
         *  \c dual member gives the "dual" relative index, i.e., the 
         *  index of \a i in \a j 's neighbor list.
         *
         *  \code
         *  Neighbor n = nb(i,_j);
         *  n.node == j &&
         *  n.iter == _j &&
         *  nb(n.node,n.dual).node == i
         *  \endcode
         *
         *  There is no easy way to transform a pair of absolute node
         *  indices \a i and \a j into a Neighbor structure relative
         *  to one of the nodes. Such a feature has never yet been
         *  found to be necessary. Iteration over edges can always be
         *  accomplished using the Neighbor lists, and by writing
         *  functions that accept relative indices:
         *  \code
         *  for( size_t i = 0; i < nrNodes(); ++i )
         *      foreach( const Neighbor &j, nb(i) )
         *          T( i, j.iter );
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
            /// Constructor that sets the Neighbor members according to the parameters
            Neighbor( size_t iter, size_t node, size_t dual ) : iter(iter), node(node), dual(dual) {}

            /// Cast to \c size_t returns \c node member
            operator size_t () const { return node; }
        };

        /// Describes the neighbors of some node.
        typedef std::vector<Neighbor> Neighbors;

        /// Represents an edge: an Edge(\a n1,\a n2) corresponds to the edge between node \a n1 and node \a n2.
        typedef std::pair<size_t,size_t> Edge;

    private:
        /// Contains for each node a vector of its neighbors
        std::vector<Neighbors> _nb;

    public:
    /// \name Constructors and destructors
    //@{
        /// Default constructor (creates an empty graph).
        GraphAL() : _nb() {}

        /// Constructs GraphAL with \a nr nodes and no edges.
        GraphAL( size_t nr ) : _nb( nr ) {}

        /// Constructs GraphAL from a range of edges.
        /** \tparam EdgeInputIterator Iterator that iterates over instances of GraphAL::Edge.
         *  \param nr The number of nodes.
         *  \param begin Points to the first edge.
         *  \param end Points just beyond the last edge.
         */
        template<typename EdgeInputIterator>
        GraphAL( size_t nr, EdgeInputIterator begin, EdgeInputIterator end ) : _nb() {
            construct( nr, begin, end );
        }
    //@}

    /// \name Accessors and mutators
    //@{
        /// Returns constant reference to the \a _n2 'th neighbor of node \a n1
        const Neighbor & nb( size_t n1, size_t _n2 ) const {
            DAI_DEBASSERT( n1 < _nb.size() );
            DAI_DEBASSERT( _n2 < _nb[n1].size() );
            return _nb[n1][_n2];
        }
        /// Returns reference to the \a _n2 'th neighbor of node \a n1
        Neighbor & nb( size_t n1, size_t _n2 ) {
            DAI_DEBASSERT( n1 < _nb.size() );
            DAI_DEBASSERT( _n2 < _nb[n1].size() );
            return _nb[n1][_n2];
        }

        /// Returns constant reference to all neighbors of node \a n
        const Neighbors & nb( size_t n ) const {
            DAI_DEBASSERT( n < _nb.size() );
            return _nb[n];
        }
        /// Returns reference to all neighbors of node \a n
        Neighbors & nb( size_t n ) {
            DAI_DEBASSERT( n < _nb.size() );
            return _nb[n];
        }
    //@}

    /// \name Adding nodes and edges
    //@{
        /// (Re)constructs GraphAL from a range of edges.
        /** \tparam EdgeInputIterator Iterator that iterates over instances of GraphAL::Edge.
         *  \param nr The number of nodes.
         *  \param begin Points to the first edge.
         *  \param end Points just beyond the last edge.
         */
        template<typename EdgeInputIterator>
        void construct( size_t nr, EdgeInputIterator begin, EdgeInputIterator end );

        /// Adds a node without neighbors and returns the index of the added node.
        size_t addNode() { _nb.push_back( Neighbors() ); return _nb.size() - 1; }

        /// Adds a node, with neighbors specified by a range of nodes.
        /** \tparam NodeInputIterator Iterator that iterates over instances of \c size_t.
         *  \param begin Points to the first index of the nodes that should become neighbors of the added node.
         *  \param end Points just beyond the last index of the nodes that should become neighbors of the added node.
         *  \param sizeHint For improved efficiency, the size of the range may be specified by \a sizeHint.
         *  \returns Index of the added node.
         */
        template <typename NodeInputIterator>
        size_t addNode( NodeInputIterator begin, NodeInputIterator end, size_t sizeHint = 0 ) {
            Neighbors nbsnew;
            nbsnew.reserve( sizeHint );
            size_t iter = 0;
            for( NodeInputIterator it = begin; it != end; ++it ) {
                DAI_ASSERT( *it < nrNodes() );
                Neighbor nb1new( iter, *it, nb(*it).size() );
                Neighbor nb2new( nb(*it).size(), nrNodes(), iter++ );
                nbsnew.push_back( nb1new );
                nb( *it ).push_back( nb2new );
            }
            _nb.push_back( nbsnew );
            return _nb.size() - 1;
        }

        /// Adds an edge between node \a n1 and node \a n2.
        /** If \a check == \c true, only adds the edge if it does not exist already.
         */
        void addEdge( size_t n1, size_t n2, bool check = true );
    //@}

    /// \name Erasing nodes and edges
    //@{
        /// Removes node \a n and all incident edges; indices of other nodes are changed accordingly.
        void eraseNode( size_t n );

        /// Removes edge between node \a n1 and node \a n2.
        void eraseEdge( size_t n1, size_t n2 );
    //@}

    /// \name Queries
    //@{
        /// Returns number of nodes
        size_t nrNodes() const { return _nb.size(); }

        /// Calculates the number of edges, time complexity: O(nrNodes())
        size_t nrEdges() const {
            size_t sum = 0;
            for( size_t i = 0; i < nrNodes(); i++ )
                sum += nb(i).size();
            return sum / 2;
        }

        /// Returns true if the graph contains an edge between nodes \a n1 and \a n2
        /** \note The time complexity is linear in the number of neighbors of \a n1 or \a n2
         */
        bool hasEdge( size_t n1, size_t n2 ) {
            if( nb(n1).size() < nb(n2).size() ) {
                for( size_t _n2 = 0; _n2 < nb(n1).size(); _n2++ )
                    if( nb( n1, _n2 ) == n2 )
                        return true;
            } else {
                for( size_t _n1 = 0; _n1 < nb(n2).size(); _n1++ )
                    if( nb( n2, _n1 ) == n1 )
                        return true;
            }
            return false;
        }

        /// Returns the index of a given node \a n2 amongst the neighbors of \a n1
        /** \note The time complexity is linear in the number of neighbors of \a n1
         *  \throw OBJECT_NOT_FOUND if \a n2 is not a neighbor of \a n1
         */
        size_t findNb( size_t n1, size_t n2 ) {
            for( size_t _n2 = 0; _n2 < nb(n1).size(); _n2++ )
                if( nb( n1, _n2 ) == n2 )
                    return _n2;
            DAI_THROW(OBJECT_NOT_FOUND);
            return nb(n1).size();
        }

        /// Returns true if the graph is connected
        bool isConnected() const;

        /// Returns true if the graph is a tree, i.e., if it is singly connected and connected.
        bool isTree() const;

        /// Asserts internal consistency
        void checkConsistency() const;
    //@}

    /// \name Input and output
    //@{
        /// Writes this GraphAL to an output stream in GraphViz .dot syntax
        void printDot( std::ostream& os ) const;
    //@}
};


template<typename EdgeInputIterator>
void GraphAL::construct( size_t nr, EdgeInputIterator begin, EdgeInputIterator end ) {
    _nb.clear();
    _nb.resize( nr );

    for( EdgeInputIterator e = begin; e != end; e++ ) {
#ifdef DAI_DEBUG
        addEdge( e->first, e->second, true );
#else
        addEdge( e->first, e->second, false );
#endif
    }
}


/// Creates a fully-connected graph with \a N nodes
GraphAL createGraphFull( size_t N );
/// Creates a two-dimensional rectangular grid of \a N1 by \a N2 nodes, which can be \a periodic
GraphAL createGraphGrid( size_t N1, size_t N2, bool periodic );
/// Creates a three-dimensional rectangular grid of \a N1 by \a N2 by \a N3 nodes, which can be \a periodic
GraphAL createGraphGrid3D( size_t N1, size_t N2, size_t N3, bool periodic );
/// Creates a graph consisting of a single loop of \a N nodes
GraphAL createGraphLoop( size_t N );
/// Creates a random tree-structured graph of \a N nodes
GraphAL createGraphTree( size_t N );
/// Creates a random regular graph of \a N nodes with uniform connectivity \a d
/** Algorithm 1 in [\ref StW99].
 *  Draws a random graph of size \a N and uniform degree \a d
 *  from an almost uniform probability distribution over these graphs
 *  (which becomes uniform in the limit that \a d is small and \a N goes
 *  to infinity).
 */
GraphAL createGraphRegular( size_t N, size_t d );


} // end of namespace dai


#endif
