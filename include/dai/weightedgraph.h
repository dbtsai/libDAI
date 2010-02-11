/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/** \file
 *  \brief Defines some utility functions for (weighted) undirected graphs, trees and rooted trees.
 *  \todo Improve general support for graphs and trees.
 */


#ifndef __defined_libdai_weightedgraph_h
#define __defined_libdai_weightedgraph_h


#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <limits>
#include <climits>   // Work-around for bug in boost graph library

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>


namespace dai {


/// Represents a directed edge
class DEdge {
    public:
        /// First node index (source of edge)
        size_t n1;

        /// Second node index (sink of edge)
        size_t n2;

        /// Default constructor
        DEdge() {}

        /// Constructs a directed edge pointing from \a m1 to \a m2
        DEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}

        /// Tests for equality
        bool operator==( const DEdge &x ) const { return ((n1 == x.n1) && (n2 == x.n2)); }

        /// Tests for inequality
        bool operator!=( const DEdge &x ) const { return !(*this == x); }

        /// Smaller-than operator (performs lexicographical comparison)
        bool operator<( const DEdge &x ) const {
            return( (n1 < x.n1) || ((n1 == x.n1) && (n2 < x.n2)) );
        }

        /// Writes a directed edge to an output stream
        friend std::ostream & operator << (std::ostream & os, const DEdge & e) {
            os << "(" << e.n1 << "," << e.n2 << ")";
            return os;
        }
};


/// Represents an undirected edge
class UEdge {
    public:
        /// First node index
        union {
            size_t n1;
            size_t first;   /// alias
        };
        /// Second node index
        union {
            size_t n2;
            size_t second;  /// alias
        };

        /// Default constructor
        UEdge() {}

        /// Constructs an undirected edge between \a m1 and \a m2
        UEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}

        /// Construct from DEdge
        UEdge( const DEdge &e ) : n1(e.n1), n2(e.n2) {}

        /// Tests for inequality (disregarding the ordering of the nodes)
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

        /// Writes an undirected edge to an output stream
        friend std::ostream & operator << (std::ostream & os, const UEdge & e) {
            if( e.n1 < e.n2 )
                os << "{" << e.n1 << "," << e.n2 << "}";
            else
                os << "{" << e.n2 << "," << e.n1 << "}";
            return os;
        }
};


/// Represents an undirected graph, implemented as a std::set of undirected edges
class GraphEL : public std::set<UEdge> {
    public:
        /// Default constructor
        GraphEL() {}

        /// Construct from range of objects that can be cast to DEdge
        template <class InputIterator>
        GraphEL( InputIterator begin, InputIterator end ) {
            insert( begin, end );
        }
};


/// Represents an undirected weighted graph, with weights of type \a T, implemented as a std::map mapping undirected edges to weights
template<class T> class WeightedGraph : public std::map<UEdge, T> {};


/// Represents a rooted tree, implemented as a vector of directed edges
/** By convention, the edges are stored such that they point away from 
 *  the root and such that edges nearer to the root come before edges
 *  farther away from the root.
 */
class RootedTree : public std::vector<DEdge> {
    public:
        /// Default constructor
        RootedTree() {}

        /// Constructs a rooted tree from a tree and a root
        /** \pre T has no cycles and contains node \a Root
         */
        RootedTree( const GraphEL &T, size_t Root );
};


/// Uses Prim's algorithm to construct a minimal spanning tree from the (positively) weighted graph \a G.
/** Uses implementation in Boost Graph Library.
 */
template<typename T> RootedTree MinSpanningTreePrims( const WeightedGraph<T> &G ) {
    RootedTree result;
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
        GraphEL tree;
        size_t root = 0;
        for( size_t i = 0; i != p.size(); i++ )
            if( p[i] != i )
                tree.insert( UEdge( p[i], i ) );
            else
                root = i;
        // Order them to obtain a rooted tree
        result = RootedTree( tree, root );
    }
    return result;
}


/// Use Prim's algorithm to construct a maximal spanning tree from the (positively) weighted graph \a G.
/** Uses implementation in Boost Graph Library.
 */
template<typename T> RootedTree MaxSpanningTreePrims( const WeightedGraph<T> &G ) {
    if( G.size() == 0 )
        return RootedTree();
    else {
        T maxweight = G.begin()->second;
        for( typename WeightedGraph<T>::const_iterator it = G.begin(); it != G.end(); it++ )
            if( it->second > maxweight )
                maxweight = it->second;
        // make a copy of the graph
        WeightedGraph<T> gr( G );
        // invoke MinSpanningTreePrims with negative weights
        // (which have to be shifted to satisfy positivity criterion)
        for( typename WeightedGraph<T>::iterator it = gr.begin(); it != gr.end(); it++ )
            it->second = maxweight - it->second;
        return MinSpanningTreePrims( gr );
    }
}


/// Constructs a random undirected graph of \a N nodes, where each node has connectivity \a d
/** Algorithm 1 in [\ref StW99].
 *  Draws a random graph of size \a N and uniform degree \a d
 *  from an almost uniform probability distribution over these graphs
 *  (which becomes uniform in the limit that \a d is small and \a N goes
 *  to infinity).
 */
GraphEL RandomDRegularGraph( size_t N, size_t d );


} // end of namespace dai


#endif
