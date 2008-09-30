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
/// \brief Defines class ClusterGraph


#ifndef __defined_libdai_clustergraph_h
#define __defined_libdai_clustergraph_h


#include <set>
#include <vector>
#include <dai/varset.h>
#include <dai/bipgraph.h>


namespace dai {


    /// A ClusterGraph is a hypergraph with VarSets as nodes.
    /** It is implemented as bipartite graph with variable (Var) nodes
     *  and cluster (VarSet) nodes. 
     */
    class ClusterGraph {
        public:
            /// Stores the neighborhood structure
            BipartiteGraph       G;

            /// Stores the variables corresponding to the nodes
            std::vector<Var>     vars;

            /// Stores the clusters corresponding to the hyperedges
            std::vector<VarSet>  clusters;

            /// Shorthand for BipartiteGraph::Neighbor
            typedef BipartiteGraph::Neighbor Neighbor;

            /// Shorthand for BipartiteGraph::Edge
            typedef BipartiteGraph::Edge     Edge;

        public:
            /// Default constructor
            ClusterGraph() : G(), vars(), clusters() {}

            /// Construct from vector<VarSet>
            ClusterGraph( const std::vector<VarSet> & cls );
            
            /// Copy constructor
            ClusterGraph( const ClusterGraph &x ) : G(x.G), vars(x.vars), clusters(x.clusters) {}

            /// Assignment operator
            ClusterGraph& operator=( const ClusterGraph &x ) {
                if( this != &x ) {
                    G = x.G;
                    vars = x.vars;
                    clusters = x.clusters;
                }
                return *this;
            }

            /// Returns true if cluster I is not contained in a larger cluster
            bool isMaximal( size_t I ) const {
#ifdef DAI_DEBUG
                assert( I < G.nr2() );
#endif
                const VarSet & clI = clusters[I];
                bool maximal = true;
                // The following may not be optimal, since it may repeatedly test the same cluster *J
                foreach( const Neighbor &i, G.nb2(I) ) {
                    foreach( const Neighbor &J, G.nb1(i) )
                        if( (J != I) && (clI << clusters[J]) ) {
                            maximal = false;
                            break;
                        }
                    if( !maximal )
                        break;
                }
                return maximal;
            }

            /// Erases all VarSets that are not maximal
            ClusterGraph& eraseNonMaximal() {
                for( size_t I = 0; I < G.nr2(); ) {
                    if( !isMaximal(I) ) {
                        clusters.erase( clusters.begin() + I );
                        G.erase2(I);
                    } else
                        I++;
                }
                return *this;
            }

            /// Returns number of clusters
            size_t size() const {
                return G.nr2();
            }

            /// Returns index of variable n
            size_t findVar( const Var &n ) const {
                return find( vars.begin(), vars.end(), n ) - vars.begin();
            }

            /// Returns true if vars with indices i1 and i2 are adjacent, i.e., both contained in the same cluster
            bool adj( size_t i1, size_t i2 ) {
                bool result = false;
                foreach( const Neighbor &I, G.nb1(i1) )
                    if( find( G.nb2(I).begin(), G.nb2(I).end(), i2 ) != G.nb2(I).end() ) {
                        result = true;
                        break;
                    }
                return result;
            }
            
            /// Returns union of clusters that contain the variable with index i
            VarSet Delta( size_t i ) const {
                VarSet result;
                foreach( const Neighbor &I, G.nb1(i) )
                    result |= clusters[I];
                return result;
            }

            /// Inserts a cluster (if it does not already exist)
            void insert( const VarSet &cl ) {
                if( find( clusters.begin(), clusters.end(), cl ) == clusters.end() ) {
                    clusters.push_back( cl );
                    // add variables (if necessary) and calculate neighborhood of new cluster
                    std::vector<size_t> nbs;
                    for( VarSet::const_iterator n = cl.begin(); n != cl.end(); n++ ) {
                        size_t iter = find( vars.begin(), vars.end(), *n ) - vars.begin();
                        nbs.push_back( iter );
                        if( iter == vars.size() ) {
                            G.add1();
                            vars.push_back( *n );
                        }
                    }
                    G.add2( nbs.begin(), nbs.end(), nbs.size() );
                }
            }

            /// Returns union of clusters that contain variable with index i, minus this variable
            VarSet delta( size_t i ) const {
                return Delta( i ) / vars[i];
            }
            
            /// Erases all clusters that contain n where n is the variable with index i
            ClusterGraph& eraseSubsuming( size_t i ) {
                while( G.nb1(i).size() ) {
                    clusters.erase( clusters.begin() + G.nb1(i)[0] );
                    G.erase2( G.nb1(i)[0] );
                }
                return *this;
            }
            
            /// Returns a const reference to the clusters
            const std::vector<VarSet> & toVector() const { return clusters; }

            /// Calculates cost of eliminating the variable with index i.
            /** The cost is measured as "number of added edges in the adjacency graph",
             *  where the adjacency graph has the variables as its nodes and
             *  connects nodes i1 and i2 iff i1 and i2 occur in some common cluster.
             */
            size_t eliminationCost( size_t i ) {
                std::vector<size_t> id_n = G.delta1( i );

                size_t cost = 0;

                // for each unordered pair {i1,i2} adjacent to n
                for( size_t _i1 = 0; _i1 < id_n.size(); _i1++ )
                    for( size_t _i2 = _i1 + 1; _i2 < id_n.size(); _i2++ ) {
                        // if i1 and i2 are not adjacent, eliminating n would make them adjacent
                        if( !adj(id_n[_i1], id_n[_i2]) )
                            cost++;
                    }

                return cost;
            }

            /// Performs Variable Elimination without Probs, i.e. only keeping track of
            /*  the interactions that are created along the way.
             *  \param ElimSeq A set of outer clusters and an elimination sequence
             *  \return A set of elimination "cliques"
             */
            ClusterGraph VarElim( const std::vector<Var> &ElimSeq ) const;

            /// Performs Variable Eliminiation using the MinFill heuristic
            ClusterGraph VarElim_MinFill() const;

            /// Writes a ClusterGraph to an output stream
            friend std::ostream & operator << ( std::ostream & os, const ClusterGraph & cl ) {
                os << cl.toVector();
                return os;
            }
    };


} // end of namespace dai


#endif
