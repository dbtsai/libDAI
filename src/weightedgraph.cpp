/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <algorithm>
#include <dai/weightedgraph.h>
#include <dai/util.h>
#include <dai/exceptions.h>


namespace dai {


using namespace std;


RootedTree::RootedTree( const GraphEL &T, size_t Root ) {
    if( T.size() != 0 ) {
        // Make a copy
        GraphEL Gr = T;

        // Nodes in the tree
        set<size_t> nodes;

        // Start with the root
        nodes.insert( Root );

        // Keep adding edges until done
        while( !(Gr.empty()) )
            for( GraphEL::iterator e = Gr.begin(); e != Gr.end(); ) {
                bool e1_in_nodes = nodes.count( e->n1 );
                bool e2_in_nodes = nodes.count( e->n2 );
                DAI_ASSERT( !(e1_in_nodes && e2_in_nodes) );
                if( e1_in_nodes ) {
                    // Add directed edge, pointing away from the root
                    push_back( DEdge( e->n1, e->n2 ) );
                    nodes.insert( e->n2 );
                    // Erase the edge
                    Gr.erase( e++ );
                } else if( e2_in_nodes ) {
                    // Add directed edge, pointing away from the root
                    push_back( DEdge( e->n2, e->n1 ) );
                    nodes.insert( e->n1 );
                    // Erase the edge
                    Gr.erase( e++ );
                } else
                    e++;
            }
    }
}


GraphEL RandomDRegularGraph( size_t N, size_t d ) {
    DAI_ASSERT( (N * d) % 2 == 0 );

    bool ready = false;
    std::vector<UEdge> G;

    size_t tries = 0;
    while( !ready ) {
        tries++;

        // Start with N*d points {0,1,...,N*d-1} (N*d even) in N groups.
        // Put U = {0,1,...,N*d-1}. (U denotes the set of unpaired points.)
        vector<size_t> U;
        U.reserve( N * d );
        for( size_t i = 0; i < N * d; i++ )
            U.push_back( i );

        // Repeat the following until no suitable pair can be found: Choose
        // two random points i and j in U, and if they are suitable, pair
        // i with j and delete i and j from U.
        G.clear();
        bool finished = false;
        while( !finished ) {
            random_shuffle( U.begin(), U.end() );
            size_t i1, i2;
            bool suit_pair_found = false;
            for( i1 = 0; i1 < U.size()-1 && !suit_pair_found; i1++ )
                for( i2 = i1+1; i2 < U.size() && !suit_pair_found; i2++ )
                    if( (U[i1] / d) != (U[i2] / d) ) {
                        // they are suitable
                        suit_pair_found = true;
                        G.push_back( UEdge( U[i1] / d, U[i2] / d ) );
                        U.erase( U.begin() + i2 );  // first remove largest
                        U.erase( U.begin() + i1 );  // remove smallest
                    }
            if( !suit_pair_found || U.empty() )
                finished = true;
        }

        if( U.empty() ) {
            // G is a graph with edge from vertex r to vertex s if and only if
            // there is a pair containing points in the r'th and s'th groups.
            // If G is d-regular, output, otherwise return to Step 1.

            vector<size_t> degrees;
            degrees.resize( N, 0 );
            foreach( const UEdge &e, G ) {
                degrees[e.n1]++;
                degrees[e.n2]++;
            }
            ready = true;
            for( size_t n = 0; n < N; n++ )
                if( degrees[n] != d ) {
                    ready = false;
                    break;
                }
        } else
            ready = false;
    }

    return GraphEL( G.begin(), G.end() );
}


} // end of namespace dai
