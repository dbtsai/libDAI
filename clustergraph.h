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


#ifndef __CLUSTERGRAPH_H__
#define __CLUSTERGRAPH_H__


#include <set>
#include <vector>
#include "varset.h"


namespace dai {


using namespace std;


/// A ClusterGraph is a hypergraph with VarSets as nodes.
/// It is implemented as a set<VarSet> in which the adjacency
/// relationship is computed realtime. It may be better to
/// implement it as a bipartitegraph, though. The additional
/// functionality compared to a simple set<VarSet> is
/// finding maximal clusters, finding cliques, etc...
class ClusterGraph : public set<VarSet> {
    public:
        /// Default constructor
        ClusterGraph() : set<VarSet>() {}

        /// Construct from vector<VarSet>
        ClusterGraph( const vector<VarSet> & cls ) {
            insert( cls.begin(), cls.end() );
        }
        
        /// Copy constructor
        ClusterGraph(const ClusterGraph &x) : set<VarSet>(x) {}

        /// Assignment operator
        ClusterGraph& operator=( const ClusterGraph &x ) {
            if( this != &x ) {
                set<VarSet>::operator=( x );
            }
            return *this;
        }

        /// Returns true if ns is a maximal member of *this under inclusion (VarSet::operator<<)
        bool isMaximal( const VarSet &ns ) const {
            if( count( ns ) ) {
                // ns is a member
                bool maximal = true;
                for( const_iterator x = begin(); x != end() && maximal; x++ )
                    if( (ns << *x) && (ns != *x) )
                        maximal = false;
                return maximal;
            } else
                return false;
        }

        /// Erase all VarSets that are not maximal
        ClusterGraph& eraseNonMaximal() {
            for( iterator x = begin(); x != end(); )
                if( !isMaximal(*x) )
                    erase(x++);
                else
                    x++;
            return *this;
        }

        /// Return union of all members
        VarSet vars() const {
            VarSet result;
            for( iterator x = begin(); x != end(); x++ )
                result |= *x;
            return result;
        }

        /// Returns true if n1 and n2 are adjacent, i.e.\ by
        /// definition, are both contained in some cluster in *this
        bool adj( const Var& n1, const Var& n2 ) {
            bool result = false;
            for( iterator x = begin(); (x != end()) && (!result); x++ )
                if( (*x && n1) && (*x && n2) )
                    result = true;
            return result;
        }
        
        /// Returns union of clusters that contain n, minus n
        VarSet Delta( const Var& n ) const {
            VarSet result;
            for( iterator x = begin(); x != end(); x++ )
                if( (*x && n) )
                    result |= *x;
            return result;
        }

        /// Returns union of clusters that contain n, minus n
        VarSet delta( const Var& n ) const {
            return Delta( n ) / n;
        }
        
        /// Erases all members that contain n
        ClusterGraph& eraseSubsuming( const Var& n ) {
            for( iterator x = begin(); x != end(); )
                if( (*x && n) )
                    erase(x++);
                else
                    x++;
            return *this;
        }
        
        /// Send to output stream
        friend std::ostream & operator << ( std::ostream & os, const ClusterGraph & cl ) {
            os << "{";
            ClusterGraph::const_iterator x = cl.begin();
            if( x != cl.end() )
                os << *(x++);
            for( ; x != cl.end(); x++ )
                os << ", " << *x;
            os << "}";
            return os;
        }

        /// Convert to vector<VarSet>
        vector<VarSet> toVector() const {
            vector<VarSet> result;
            result.reserve( size() );
            for( const_iterator x = begin(); x != end(); x++ )
                result.push_back( *x );
            return result;
        }

        /// Calculate cost of eliminating variable n
        /// using as a measure "number of added edges in the adjacency graph"
        /// where the adjacency graph has the variables as its nodes and
        /// connects nodes i1 and i2 iff i1 and i2 occur in some common factor I
        size_t eliminationCost( const Var& n ) {
            VarSet d_n = delta( n );
            size_t cost = 0;

            // for each unordered pair {i1,i2} adjacent to n
            for( VarSet::const_iterator i1 = d_n.begin(); i1 != d_n.end(); i1++ ) {
                VarSet d_i1 = delta( *i1 );
                for( VarSet::const_iterator i2 = i1; (++i2) != d_n.end(); ) {
                    // if i1 and i2 are not adjacent, eliminating n would make them adjacent
                    if( !adj(*i1, *i2) )
                        cost++;
                }
            }
            return cost;
        }

        /// Perform Variable Elimination without Probs, i.e. only keeping track of
        /// the interactions that are created along the way.
        /// Input:  a set of outer clusters and an elimination sequence
        /// Output: a set of elimination "cliques"
        ClusterGraph VarElim( const vector<Var> &ElimSeq ) const;

        /// As Taylan does it
        ClusterGraph VarElim_MinFill() const;
};


}


#endif
