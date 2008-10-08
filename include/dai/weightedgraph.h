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
/// \brief Defines some utility functions for weighted graphs
/// \todo Improve documentation


#ifndef __defined_libdai_weightedgraph_h
#define __defined_libdai_weightedgraph_h


#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <cassert>
#include <limits>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>


namespace dai {


/// Represents a directed edge pointing from n1 to n2
class DEdge {
    public:
        size_t n1;  ///< First node index
        size_t n2;  ///< Second node index
    
        /// Default constructor
        DEdge() {}

        /// Constructor
        DEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}

        /// Tests for equality
        bool operator==( const DEdge &x ) const { return ((n1 == x.n1) && (n2 == x.n2)); }

        /// Tests for inequality
        bool operator!=( const DEdge &x ) const { return !(*this == x); }

        /// Smaller-than operator (performs lexicographical comparison)
        bool operator<( const DEdge &x ) const {
            return( (n1 < x.n1) || ((n1 == x.n1) && (n2 < x.n2)) );
        }

        /// Writes a DEdge to an output stream
        friend std::ostream & operator << (std::ostream & os, const DEdge & e) {
            os << "(" << e.n1 << "," << e.n2 << ")";
            return os;
        }
};


/// Undirected edge between nodes n1 and n2
class UEdge {
    public:
        size_t  n1;  ///< First node index
        size_t  n2;  ///< Second node index
    
        /// Default constructor
        UEdge() {}

        /// Constructor
        UEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}

        /// Construct from DEdge
        UEdge( const DEdge & e ) : n1(e.n1), n2(e.n2) {}

        /// Tests for inequality (disregarding the ordering of n1 and n2)
        bool operator==( const UEdge &x ) {
            return ((n1 == x.n1) && (n2 == x.n2)) || ((n1 == x.n2) && (n2 == x.n1));
        }

        /// Smaller-than operator
        bool operator<( const UEdge &x ) const {
            size_t s = n1, l = n2;
            if( s > l )
                std::swap( s, l );
            size_t xs = x.n1, xl = x.n2;
            if( xs > xl )
                std::swap( xs, xl );
            return( (s < xs) || ((s == xs) && (l < xl)) );
        }

        /// Writes a UEdge to an output stream
        friend std::ostream & operator << (std::ostream & os, const UEdge & e) {
            if( e.n1 < e.n2 )
                os << "{" << e.n1 << "," << e.n2 << "}";
            else
                os << "{" << e.n2 << "," << e.n1 << "}";
            return os;
        }
};


/// Vector of UEdge
typedef std::vector<UEdge>  UEdgeVec;

/// Vector of DEdge
typedef std::vector<DEdge>  DEdgeVec;

/// Represents an undirected weighted graph, with weights of type T
template<class T> class WeightedGraph : public std::map<UEdge, T> {};

/// Represents an undirected graph
typedef std::set<UEdge>     Graph;


/// Uses Prim's algorithm to construct a minimal spanning tree from the (positively) weighted graph G.
/** Uses implementation in Boost Graph Library.
 */
template<typename T> DEdgeVec MinSpanningTreePrims( const WeightedGraph<T> &G ) {
    DEdgeVec result;
    if( G.size() > 0 ) {
        using namespace boost;
        using namespace std;
        typedef adjacency_list< vecS, vecS, undirectedS, property<vertex_distance_t, int>, property<edge_weight_t, double> > boostGraph;
        typedef pair<size_t, size_t> E;

        set<size_t> nodes;
        vector<E> edges;
        vector<double> weights;
        edges.reserve( G.size() );
        weights.reserve( G.size() );
        for( typename WeightedGraph<T>::const_iterator e = G.begin(); e != G.end(); e++ ) {
            weights.push_back( e->second );
            edges.push_back( E( e->first.n1, e->first.n2 ) );
            nodes.insert( e->first.n1 );
            nodes.insert( e->first.n2 );
        }

        boostGraph g( edges.begin(), edges.end(), weights.begin(), nodes.size() );
        vector< graph_traits< boostGraph >::vertex_descriptor > p( num_vertices(g) );
        prim_minimum_spanning_tree( g, &(p[0]) );

        // Store tree edges in result
        result.reserve( nodes.size() - 1 );
        size_t root = 0;
        for( size_t i = 0; i != p.size(); i++ )
            if( p[i] != i )
                result.push_back( DEdge( p[i], i ) );
            else
                root = i;

        // We have to store the minimum spanning tree in the right
        // order, such that for all (i1, j1), (i2, j2) in result,
        // if j1 == i2 then (i1, j1) comes before (i2, j2) in result.
        // We do this by reordering the contents of result, effectively
        // growing the tree starting at the root. At each step, 
        // result[0..N-1] are the edges already added to the tree,
        // whereas the other elements of result still have to be added.
        // The elements of nodes are the vertices that still have to
        // be added to the tree.

        // Start with the root
        nodes.erase( root );
        size_t N = 0;

        // Iteratively add edges and nodes to the growing tree
        while( N != result.size() ) {
            for( size_t e = N; e != result.size(); e++ ) {
                bool e1_in_tree = !nodes.count( result[e].n1 );
                if( e1_in_tree ) {
                    nodes.erase( result[e].n2 );
                    swap( result[N], result[e] );
                    N++;
                    break;
                }
            }
        }
    }

    return result;
}


/// Use Prim's algorithm to construct a minimal spanning tree from the (positively) weighted graph G.
/** Uses implementation in Boost Graph Library.
 */
template<typename T> DEdgeVec MaxSpanningTreePrims( const WeightedGraph<T> & Graph ) {
    T maxweight = Graph.begin()->second;
    for( typename WeightedGraph<T>::const_iterator it = Graph.begin(); it != Graph.end(); it++ )
        if( it->second > maxweight )
            maxweight = it->second;
    // make a copy of the graph
    WeightedGraph<T> gr( Graph );
    // invoke MinSpanningTreePrims with negative weights
    // (which have to be shifted to satisfy positivity criterion)
    for( typename WeightedGraph<T>::iterator it = gr.begin(); it != gr.end(); it++ )
        it->second = maxweight - it->second;
    return MinSpanningTreePrims( gr );
}


/// Constructs a rooted tree from a tree and a root
DEdgeVec GrowRootedTree( const Graph & T, size_t Root );


/// Constructs a random undirected graph of N nodes, where each node has connectivity d
UEdgeVec RandomDRegularGraph( size_t N, size_t d );


} // end of namespace dai


#endif
