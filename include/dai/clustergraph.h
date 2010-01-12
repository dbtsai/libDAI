/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines class ClusterGraph, which is used by JTree, TreeEP and HAK


#ifndef __defined_libdai_clustergraph_h
#define __defined_libdai_clustergraph_h


#include <set>
#include <vector>
#include <dai/varset.h>
#include <dai/bipgraph.h>


namespace dai {


    /// A ClusterGraph is a hypergraph with variables as nodes, and "clusters" (sets of variables) as hyperedges.
    /** It is implemented as a bipartite graph with variable (Var) nodes and cluster (VarSet) nodes.
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
        /// \name Constructors and destructors
        //@{
            /// Default constructor
            ClusterGraph() : G(), vars(), clusters() {}

            /// Construct from vector of VarSet 's
            ClusterGraph( const std::vector<VarSet> & cls );
        //@}

        /// \name Queries
        //@{
            /// Returns a constant reference to the clusters
            const std::vector<VarSet> & toVector() const { return clusters; }

            /// Returns number of clusters
            size_t size() const {
                return G.nrNodes2();
            }

            /// Returns the index of variable \a n
            size_t findVar( const Var &n ) const {
                return find( vars.begin(), vars.end(), n ) - vars.begin();
            }

            /// Returns union of clusters that contain the \a i 'th variable
            VarSet Delta( size_t i ) const {
                VarSet result;
                foreach( const Neighbor &I, G.nb1(i) )
                    result |= clusters[I];
                return result;
            }

            /// Returns union of clusters that contain the \a i 'th (except this variable itself)
            VarSet delta( size_t i ) const {
                return Delta( i ) / vars[i];
            }

            /// Returns \c true if variables with indices \a i1 and \a i2 are adjacent, i.e., both contained in the same cluster
            bool adj( size_t i1, size_t i2 ) {
                bool result = false;
                foreach( const Neighbor &I, G.nb1(i1) )
                    if( find( G.nb2(I).begin(), G.nb2(I).end(), i2 ) != G.nb2(I).end() ) {
                        result = true;
                        break;
                    }
                return result;
            }

            /// Returns \c true if cluster \a I is not contained in a larger cluster
            bool isMaximal( size_t I ) const {
                DAI_DEBASSERT( I < G.nrNodes2() );
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
        //@}

        /// \name Operations
        //@{
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
                            G.addNode1();
                            vars.push_back( *n );
                        }
                    }
                    G.addNode2( nbs.begin(), nbs.end(), nbs.size() );
                }
            }

            /// Erases all clusters that are not maximal
            ClusterGraph& eraseNonMaximal() {
                for( size_t I = 0; I < G.nrNodes2(); ) {
                    if( !isMaximal(I) ) {
                        clusters.erase( clusters.begin() + I );
                        G.eraseNode2(I);
                    } else
                        I++;
                }
                return *this;
            }

            /// Erases all clusters that contain the \a i 'th variable
            ClusterGraph& eraseSubsuming( size_t i ) {
                while( G.nb1(i).size() ) {
                    clusters.erase( clusters.begin() + G.nb1(i)[0] );
                    G.eraseNode2( G.nb1(i)[0] );
                }
                return *this;
            }
        //@}

        /// \name Input/Ouput
        //@{
            /// Writes a ClusterGraph to an output stream
            friend std::ostream & operator << ( std::ostream & os, const ClusterGraph & cl ) {
                os << cl.toVector();
                return os;
            }
        //@}

        /// \name Variable elimination
        //@{
            /// Calculates cost of eliminating the \a i 'th variable.
            /** The cost is measured as "number of added edges in the adjacency graph",
             *  where the adjacency graph has the variables as its nodes and connects
             *  nodes \a i1 and \a i2 iff \a i1 and \a i2 occur together in some common cluster.
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

            /// Performs Variable Elimination, only keeping track of the interactions that are created along the way.
            /** \param ElimSeq The sequence in which to eliminate the variables
             *  \return A set of elimination "cliques"
             *  \todo Variable elimination should be implemented generically using a function
             *  object that specifies which variable to delete next.
             */
            ClusterGraph VarElim( const std::vector<Var> &ElimSeq ) const;

            /// Performs Variable Eliminiation using the "MinFill" heuristic
            /** The "MinFill" heuristic greedily minimizes the cost of eliminating a variable,
             *  measured with eliminationCost().
             *  \return A set of elimination "cliques"
             */
            ClusterGraph VarElim_MinFill() const;
        //@}
    };


} // end of namespace dai


#endif
