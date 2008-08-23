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


#ifndef __WEIGHTEDGRAPH_H__
#define __WEIGHTEDGRAPH_H__


#include <vector>
#include <map>
#include <iostream>
#include <set>


using namespace std;


/// Directed edge
class DEdge {
    public:
        size_t  n1, n2;
    
        DEdge() {}
        DEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}
        bool operator==( const DEdge &x ) const {
            return ((n1 == x.n1) && (n2 == x.n2));
        }
        bool operator!=( const DEdge &x ) const {
            return !(*this == x);
        }
        bool operator<( const DEdge &x ) const {
            return( (n1 < x.n1) || ((n1 == x.n1) && (n2 < x.n2)) );
        }
        friend std::ostream & operator << (std::ostream & os, const DEdge & e) {
            os << "(" << e.n1 << "," << e.n2 << ")";
            return os;
        }
};


/// Undirected edge
class UEdge {
    public:
        size_t  n1, n2;
    
        UEdge() {}
        UEdge( size_t m1, size_t m2 ) : n1(m1), n2(m2) {}
        UEdge( const DEdge & e ) : n1(e.n1), n2(e.n2) {}
        bool operator==( const UEdge &x ) {
            return ((n1 == x.n1) && (n2 == x.n2)) || ((n1 == x.n2) && (n2 == x.n1));
        }
        bool operator<( const UEdge &x ) const {
            size_t s = n1, l = n2;
            if( s > l )
                swap( s, l );
            size_t xs = x.n1, xl = x.n2;
            if( xs > xl )
                swap( xs, xl );
            return( (s < xs) || ((s == xs) && (l < xl)) );
        }
        friend std::ostream & operator << (std::ostream & os, const UEdge & e) {
            if( e.n1 < e.n2 )
                os << "{" << e.n1 << "," << e.n2 << "}";
            else
                os << "{" << e.n2 << "," << e.n1 << "}";
            return os;
        }
};


typedef vector<UEdge>                       UEdgeVec;
typedef vector<DEdge>                       DEdgeVec;
std::ostream & operator << (std::ostream & os, const DEdgeVec & rt);
template<class T> class WeightedGraph : public map<UEdge, T> {};
typedef set<UEdge>                          Graph;


/// Use Prim's algorithm to construct a maximal spanning tree from the weighted graph Graph
template<typename T> DEdgeVec MaxSpanningTreePrim( const WeightedGraph<T> & Graph ) {
    const long verbose = 0;

    DEdgeVec result;
    if( Graph.size() == 0 )
        return result;
    else {
        // Make a copy
        WeightedGraph<T> Gr = Graph;

        // Nodes in the tree
        set<size_t> treeV;

        // Start with one node
        treeV.insert( Gr.begin()->first.n1 );
        
        // Perform Prim's algorithm
        while( Gr.size() ) {
            typename WeightedGraph<T>::iterator largest = Gr.end();
            
            for( typename WeightedGraph<T>::iterator e = Gr.begin(); e != Gr.end(); ) {
                if( verbose >= 1 )
                    cout << "considering edge " << e->first << "...";
                bool e1_in_treeV = treeV.count( e->first.n1 );
                bool e2_in_treeV = treeV.count( e->first.n2 );
                if( e1_in_treeV && e2_in_treeV ) {
                    if( verbose >= 1 )
                        cout << "red";
                    Gr.erase( e++ );    // Nice trick! 
                } else if( e1_in_treeV || e2_in_treeV ) {
                    if( verbose >= 1 )
                        cout << e->second;
                    if( (largest == Gr.end()) || (e->second > largest->second) ) {
                        largest = e;    // largest edge connected to the tree (until now)
                        if( verbose >= 1 )
                            cout << " and largest!";
                    } 
                    e++;
                } else {
                    if( verbose >= 1 )
                        cout << "out of reach";
                    e++;
                }
                if( verbose >= 1 )
                    cout << endl;
            }

            if( largest != Gr.end() ) {
                if( verbose >= 1 )
                    cout << "largest = " << largest->first << endl;
                // Add directed edge, pointing away from the root
                if( treeV.count( largest->first.n1 ) ) {
                    result.push_back( DEdge( largest->first.n1, largest->first.n2 ) );
                    treeV.insert( largest->first.n2 );
                } else {
                    result.push_back( DEdge( largest->first.n2, largest->first.n1 ) );
                    treeV.insert( largest->first.n1 );
                }
                Gr.erase( largest );
            }
        }

        return result;
    }
}


/// Calculate rooted tree from a tree T and a root
DEdgeVec GrowRootedTree( const Graph & T, size_t Root );


UEdgeVec RandomDRegularGraph( size_t N, size_t d );


#endif
