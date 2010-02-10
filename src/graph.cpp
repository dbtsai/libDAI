/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <dai/graph.h>
#include <dai/weightedgraph.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


namespace dai {


using namespace std;


void GraphAL::addEdge( size_t n1, size_t n2, bool check ) {
    DAI_ASSERT( n1 < nrNodes() );
    DAI_ASSERT( n2 < nrNodes() );
    bool exists = false;
    if( check ) {
        // Check whether the edge already exists
        foreach( const Neighbor &n, nb(n1) )
            if( n == n2 ) {
                exists = true;
                break;
            }
    }
    if( !exists ) { // Add edge
        Neighbor nb_1( nb(n1).size(), n2, nb(n2).size() );
        Neighbor nb_2( nb_1.dual, n1, nb_1.iter );
        nb(n1).push_back( nb_1 );
        nb(n2).push_back( nb_2 );
    }
}


void GraphAL::eraseNode( size_t n ) {
    DAI_ASSERT( n < nrNodes() );
    // Erase neighbor entry of node n
    _nb.erase( _nb.begin() + n );
    // Adjust neighbor entries of nodes
    for( size_t n2 = 0; n2 < nrNodes(); n2++ ) {
        for( size_t iter = 0; iter < nb(n2).size(); ) {
            Neighbor &m = nb(n2, iter);
            if( m.node == n ) {
                // delete this entry, because it points to the deleted node
                nb(n2).erase( nb(n2).begin() + iter );
            } else if( m.node > n ) {
                // update this entry and the corresponding dual of the neighboring node
                m.iter = iter;
                m.node--;
                nb( m.node, m.dual ).dual = iter;
                iter++;
            } else {
                // skip
                iter++;
            }
        }
    }
}


void GraphAL::eraseEdge( size_t n1, size_t n2 ) {
    DAI_ASSERT( n1 < nrNodes() );
    DAI_ASSERT( n2 < nrNodes() );
    size_t iter;
    // Search for edge among neighbors of n1
    for( iter = 0; iter < nb(n1).size(); iter++ )
        if( nb(n1, iter).node == n2 ) {
            // Remove it
            nb(n1).erase( nb(n1).begin() + iter );
            break;
        }
    // Change the iter and dual values of the subsequent neighbors
    for( ; iter < nb(n1).size(); iter++ ) {
        Neighbor &m = nb( n1, iter );
        m.iter = iter;
        nb( m.node, m.dual ).dual = iter;
    }
    // Search for edge among neighbors of n2
    for( iter = 0; iter < nb(n2).size(); iter++ )
        if( nb(n2, iter).node == n1 ) {
            // Remove it
            nb(n2).erase( nb(n2).begin() + iter );
            break;
        }
    // Change the iter and node values of the subsequent neighbors
    for( ; iter < nb(n2).size(); iter++ ) {
        Neighbor &m = nb( n2, iter );
        m.iter = iter;
        nb( m.node, m.dual ).dual = iter;
    }
}


bool GraphAL::isConnected() const {
    if( nrNodes() == 0 ) {
        return true;
    } else {
        std::vector<bool> incomponent( nrNodes(), false );

        incomponent[0] = true;
        bool found_new_nodes;
        do {
            found_new_nodes = false;

            // For all nodes, check if they are connected with the (growing) component
            for( size_t n1 = 0; n1 < nrNodes(); n1++ )
                if( !incomponent[n1] ) {
                    foreach( const Neighbor &n2, nb(n1) ) {
                        if( incomponent[n2] ) {
                            found_new_nodes = true;
                            incomponent[n1] = true;
                            break;
                        }
                    }
                }
        } while( found_new_nodes );

        // Check if there are remaining nodes (not in the component)
        bool all_connected = true;
        for( size_t n1 = 0; (n1 < nrNodes()) && all_connected; n1++ )
            if( !incomponent[n1] )
                all_connected = false;

        return all_connected;

        // BGL implementation is slower...
    /*  using namespace boost;
        typedef adjacency_list< vecS, vecS, undirectedS, property<vertex_distance_t, int> > boostGraphAL;
        typedef pair<size_t, size_t> E;

        // Copy graph structure into boostGraphAL object
        vector<E> edges;
        edges.reserve( nrEdges() );
        for( size_t n1 = 0; n1 < nrNodes(); n1++ )
            foreach( const Neighbor &n2, nb(n1) )
                if( n1 < n2 )
                    edges.push_back( E( n1, n2 ) );
        boostGraphAL g( edges.begin(), edges.end(), nrNodes() );

        // Construct connected components using Boost GraphAL Library
        std::vector<int> component( num_vertices( g ) );
        int num_comp = connected_components( g, make_iterator_property_map(component.begin(), get(vertex_index, g)) );

        return (num_comp == 1);
    */
    }
}


bool GraphAL::isTree() const {
    using namespace std;
    vector<levelType> levels;

    bool foundCycle = false;
    size_t Nr = 0;

    if( nrNodes() == 0 )
        return true;
    else {
        levelType newLevel;
        do {
            newLevel.clear();
            if( levels.size() == 0 ) {
                size_t n1 = 0;
                // add n1
                newLevel = vector<size_t>( 1, n1 );
            } else {
                const levelType &prevLevel = levels.back();
                // build newLevel
                foreach( size_t n2, prevLevel ) { // for all n2 in the previous level
                    foreach( const Neighbor &n1, nb(n2) ) { // for all neighbors n1 of n2
                        if( find( prevLevel.begin(), prevLevel.end(), n1 ) == prevLevel.end() ) { // n1 not in previous level
                            if( find( newLevel.begin(), newLevel.end(), n1 ) != newLevel.end() )
                                foundCycle = true; // n1 already in new level: we found a cycle
                            else
                                newLevel.push_back( n1 ); // add n1 to new level
                        }
                        if( foundCycle )
                            break;
                    }
                    if( foundCycle )
                        break;
                }
            }
            levels.push_back( newLevel );
            Nr += newLevel.size();
        } while( (newLevel.size() != 0) && !foundCycle );
        if( Nr == nrNodes() && !foundCycle )
            return true;
        else
            return false;
    }
}


void GraphAL::printDot( std::ostream& os ) const {
    using namespace std;
    os << "graph G {" << endl;
    os << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
    for( size_t n = 0; n < nrNodes(); n++ )
        os << "\tx" << n << ";" << endl;
    for( size_t n1 = 0; n1 < nrNodes(); n1++ )
        foreach( const Neighbor &n2, nb(n1) )
            os << "\tx" << n1 << " -- y" << n2 << ";" << endl;
    os << "}" << endl;
}


void GraphAL::checkConsistency() const {
    size_t N = nrNodes();
    for( size_t n1 = 0; n1 < N; n1++ ) {
        size_t iter = 0;
        foreach( const Neighbor &n2, nb(n1) ) {
            DAI_ASSERT( n2.iter == iter );
            DAI_ASSERT( n2.node < N );
            DAI_ASSERT( n2.dual < nb(n2).size() );
            DAI_ASSERT( nb(n2, n2.dual) == n1 );
            iter++;
        }
    }
    for( size_t n2 = 0; n2 < N; n2++ ) {
        size_t iter = 0;
        foreach( const Neighbor &n1, nb(n2) ) {
            DAI_ASSERT( n1.iter == iter );
            DAI_ASSERT( n1.node < N );
            DAI_ASSERT( n1.dual < nb(n1).size() );
            DAI_ASSERT( nb(n1, n1.dual) == n2 );
            iter++;
        }
    }
}


GraphAL createGraphFull( size_t N ) {
    GraphAL result( N );
    for( size_t i = 0; i < N; i++ )
        for( size_t j = i+1; j < N; j++ )
            result.addEdge( i, j, false );
    return result;
}


GraphAL createGraphGrid( size_t n1, size_t n2, bool periodic ) {
    GraphAL result( n1*n2 );
    for( size_t i1 = 0; i1 < n1; i1++ )
        for( size_t i2 = 0; i2 < n2; i2++ ) {
            if( i1+1 < n1 || periodic )
                result.addEdge( i1*n2 + i2, ((i1+1)%n1)*n2 + i2, n1 <= 2 );
            if( i2+1 < n2 || periodic )
                result.addEdge( i1*n2 + i2, i1*n2 + ((i2+1)%n2), n2 <= 2 );
        }
    return result;
}


GraphAL createGraphGrid3D( size_t n1, size_t n2, size_t n3, bool periodic ) {
    GraphAL result( n1*n2*n3 );
    for( size_t i1 = 0; i1 < n1; i1++ )
        for( size_t i2 = 0; i2 < n2; i2++ )
            for( size_t i3 = 0; i3 < n3; i3++ ) {
                if( i1+1 < n1 || periodic )
                    result.addEdge( i1*n2*n3 + i2*n3 + i3, ((i1+1)%n1)*n2*n3 + i2*n3 + i3, n1 <= 2 );
                if( i2+1 < n2 || periodic )
                    result.addEdge( i1*n2*n3 + i2*n3 + i3, i1*n2*n3 + ((i2+1)%n2)*n3 + i3, n2 <= 2 );
                if( i3+1 < n2 || periodic )
                    result.addEdge( i1*n2*n3 + i2*n3 + i3, i1*n2*n3 + i2*n3 + ((i3+1)%n3), n3 <= 2 );
            }
    return result;
}


GraphAL createGraphLoop( size_t N ) {
    GraphAL result( N );
    for( size_t i = 0; i < N; i++ )
        result.addEdge( i, (i+1)%N, N <= 2 );
    return result;
}


GraphAL createGraphTree( size_t N ) {
    GraphAL result( N );
    for( size_t i = 1; i < N; i++ ) {
        size_t j = rnd_int( 0, i-1 );
        result.addEdge( i, j, false );
    }
    return result;
}


GraphAL createGraphRegular( size_t N, size_t d ) {
    GraphEL g = RandomDRegularGraph( N, d );
    return GraphAL( N, g.begin(), g.end() );
}


} // end of namespace dai
