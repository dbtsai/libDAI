/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/graph.h>
#include <vector>
#include <strstream>


using namespace dai;


#define BOOST_TEST_MODULE GraphALTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    GraphAL G0;
    BOOST_CHECK_EQUAL( G0.nrNodes(), 0 );
    BOOST_CHECK_EQUAL( G0.nrEdges(), 0 );
    BOOST_CHECK_EQUAL( G0.isConnected(), true );
    BOOST_CHECK_EQUAL( G0.isTree(), true );
    G0.checkConsistency();

    GraphAL G2( 2 );
    BOOST_CHECK_EQUAL( G2.nrNodes(), 2 );
    BOOST_CHECK_EQUAL( G2.nrEdges(), 0 );
    BOOST_CHECK_EQUAL( G2.isConnected(), false );
    BOOST_CHECK_EQUAL( G2.isTree(), false );
    G2.checkConsistency();
    
    std::vector<GraphAL::Edge> edges;
    edges.push_back( GraphAL::Edge( 0, 1 ) );
    edges.push_back( GraphAL::Edge( 1, 2 ) );
    edges.push_back( GraphAL::Edge( 2, 1 ) );
    edges.push_back( GraphAL::Edge( 1, 2 ) );
    GraphAL G3( 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G3.nrNodes(), 3 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 2 );
    BOOST_CHECK_EQUAL( G3.isConnected(), true );
    BOOST_CHECK_EQUAL( G3.isTree(), true );
    G3.checkConsistency();
}


BOOST_AUTO_TEST_CASE( NeighborTest ) {
    // check nb() accessor / mutator
    std::vector<GraphAL::Edge> edges;
    edges.push_back( GraphAL::Edge( 0, 1 ) );
    edges.push_back( GraphAL::Edge( 1, 2 ) );
    GraphAL G3( 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G3.nb(0).size(), 1 );
    BOOST_CHECK_EQUAL( G3.nb(1).size(), 2 );
    BOOST_CHECK_EQUAL( G3.nb(2).size(), 1 );
    BOOST_CHECK_EQUAL( G3.nb(0,0).iter, 0 );
    BOOST_CHECK_EQUAL( G3.nb(0,0).node, 1 );
    BOOST_CHECK_EQUAL( G3.nb(0,0).dual, 0 );
    BOOST_CHECK_EQUAL( G3.nb(1,0).iter, 0 );
    BOOST_CHECK_EQUAL( G3.nb(1,0).node, 0 );
    BOOST_CHECK_EQUAL( G3.nb(1,0).dual, 0 );
    BOOST_CHECK_EQUAL( G3.nb(1,1).iter, 1 );
    BOOST_CHECK_EQUAL( G3.nb(1,1).node, 2 );
    BOOST_CHECK_EQUAL( G3.nb(1,1).dual, 0 );
    BOOST_CHECK_EQUAL( G3.nb(2,0).iter, 0 );
    BOOST_CHECK_EQUAL( G3.nb(2,0).node, 1 );
    BOOST_CHECK_EQUAL( G3.nb(2,0).dual, 1 );
}


BOOST_AUTO_TEST_CASE( AddEraseTest ) {
    // check addition and erasure of nodes and edges
    std::vector<GraphAL::Edge> edges;
    edges.push_back( GraphAL::Edge( 0, 1 ) );
    edges.push_back( GraphAL::Edge( 1, 2 ) );
    edges.push_back( GraphAL::Edge( 1, 0 ) );
    GraphAL G3( 2 );
    G3.construct( 3, edges.begin(), edges.end() );
    G3.checkConsistency();
    BOOST_CHECK_EQUAL( G3.nrNodes(), 3 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 2 );
    BOOST_CHECK_EQUAL( G3.addNode(), 3 );
    G3.checkConsistency();
    std::vector<size_t> nbs;
    nbs.push_back( 3 );
    G3.addNode( nbs.begin(), nbs.end() );
    BOOST_CHECK_EQUAL( G3.addNode(), 5 );
    G3.checkConsistency();
    G3.addEdge( 0, 4 );
    G3.checkConsistency();
    G3.addEdge( 0, 5 );
    BOOST_CHECK( G3.isTree() );
    G3.checkConsistency();
    BOOST_CHECK_EQUAL( G3.nrNodes(), 6 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 5 );
    G3.addEdge( 2, 3 );
    BOOST_CHECK( !G3.isTree() );

    G3.addEdge( 5, 3 );
    G3.eraseNode( 0 );
    G3.checkConsistency();
    BOOST_CHECK( G3.isTree() );
    G3.eraseEdge( 0, 1 );
    G3.checkConsistency();
    BOOST_CHECK( !G3.isTree() );
    BOOST_CHECK( !G3.isConnected() );
    G3.eraseNode( 0 );
    G3.checkConsistency();
    BOOST_CHECK( G3.isTree() );
    G3.addEdge( 3, 2 );
    G3.checkConsistency();
    BOOST_CHECK( !G3.isTree() );
    G3.eraseNode( 1 );
    G3.checkConsistency();
    BOOST_CHECK( !G3.isTree() );
    BOOST_CHECK( !G3.isConnected() );
    G3.eraseNode( 2 );
    G3.checkConsistency();
    BOOST_CHECK( !G3.isTree() );
    BOOST_CHECK( !G3.isConnected() );
    G3.addEdge( 1, 0 );
    G3.checkConsistency();
    BOOST_CHECK( G3.isTree() );
    BOOST_CHECK( G3.isConnected() );
    G3.eraseNode( 1 );
    G3.checkConsistency();
    BOOST_CHECK( G3.isTree() );
    BOOST_CHECK( G3.isConnected() );
    G3.eraseNode( 0 );
    BOOST_CHECK( G3.isTree() );
    BOOST_CHECK( G3.isConnected() );
    BOOST_CHECK_EQUAL( G3.nrNodes(), 0 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 0 );
}


BOOST_AUTO_TEST_CASE( QueriesAndCreationTest ) {
    // check queries and createGraph* functions

    // createGraphFull
    for( size_t N = 0; N < 20; N++ ) {
        GraphAL G = createGraphFull( N );
        BOOST_CHECK_EQUAL( G.nrNodes(), N );
        BOOST_CHECK_EQUAL( G.nrEdges(), N > 0 ? N * (N-1) / 2 : 0 );
        BOOST_CHECK( G.isConnected() );
        BOOST_CHECK_EQUAL( G.isTree(), N < 3 );
        for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
            foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                BOOST_CHECK( G.hasEdge( n1, n2 ) );
                BOOST_CHECK( G.hasEdge( n2, n1 ) );
            }
            for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                if( G.hasEdge( n1, n2 ) ) {
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                }
        }
        G.checkConsistency();
    }

    // createGraphGrid
    for( size_t N1 = 0; N1 < 10; N1++ )
        for( size_t N2 = 0; N2 < 10; N2++ ) {
            GraphAL G = createGraphGrid( N1, N2, false );
            BOOST_CHECK_EQUAL( G.nrNodes(), N1 * N2 );
            BOOST_CHECK_EQUAL( G.nrEdges(), (N1 > 0 && N2 > 0) ? 2 * (N1-1) * (N2-1) + (N1-1) + (N2-1) : 0 );
            BOOST_CHECK( G.isConnected() );
            BOOST_CHECK_EQUAL( G.isTree(), (N1 <= 1) || (N2 <= 1) );
            for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
                foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                    BOOST_CHECK( G.hasEdge( n1, n2 ) );
                    BOOST_CHECK( G.hasEdge( n2, n1 ) );
                }
                for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                    if( G.hasEdge( n1, n2 ) ) {
                        BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                        BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                        BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                        BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                    }
            }
            G.checkConsistency();
            
            G = createGraphGrid( N1, N2, true );
            BOOST_CHECK_EQUAL( G.nrNodes(), N1 * N2 );
            if( N1 == 0 || N2 == 0 )
                BOOST_CHECK_EQUAL( G.nrEdges(), 0 );
            else
                BOOST_CHECK_EQUAL( G.nrEdges(), (N1 <= 2 ? (N1-1) : N1) * N2 + N1 * (N2 <= 2 ? (N2-1) : N2) );
            BOOST_CHECK( G.isConnected() );
            BOOST_CHECK_EQUAL( G.isTree(), (G.nrNodes() <= 2) );
            for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
                foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                    BOOST_CHECK( G.hasEdge( n1, n2 ) );
                    BOOST_CHECK( G.hasEdge( n2, n1 ) );
                }
                for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                    if( G.hasEdge( n1, n2 ) ) {
                        BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                        BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                        BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                        BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                    }
            }
            G.checkConsistency();
        }

    // createGraphGrid3D
    for( size_t N1 = 0; N1 < 10; N1++ )
        for( size_t N2 = 0; N2 < 10; N2++ )
            for( size_t N3 = 0; N3 < 10; N3++ ) {
                GraphAL G = createGraphGrid3D( N1, N2, N3, false );
                BOOST_CHECK_EQUAL( G.nrNodes(), N1 * N2 * N3 );
                BOOST_CHECK_EQUAL( G.nrEdges(), (N1 > 0 && N2 > 0 && N3 > 0) ? 3 * (N1-1) * (N2-1) * (N3-1) + 2 * (N1-1) * (N2-1) + 2 * (N1-1) * (N3-1) + 2 *  (N2-1) * (N3-1) + (N1-1) + (N2-1) + (N3-1) : 0 );
                BOOST_CHECK( G.isConnected() );
                BOOST_CHECK_EQUAL( G.isTree(), (G.nrNodes() == 0) || (N1 <= 1 && N2 <= 1) || (N1 <= 1 && N3 <= 1) || (N2 <= 1 && N3 <= 1) );
                for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
                    foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                        BOOST_CHECK( G.hasEdge( n1, n2 ) );
                        BOOST_CHECK( G.hasEdge( n2, n1 ) );
                    }
                    for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                        if( G.hasEdge( n1, n2 ) ) {
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                        }
                }
                G.checkConsistency();
                
                G = createGraphGrid3D( N1, N2, N3, true );
                BOOST_CHECK_EQUAL( G.nrNodes(), N1 * N2 * N3 );
                if( N1 == 0 || N2 == 0 || N3 == 0 )
                    BOOST_CHECK_EQUAL( G.nrEdges(), 0 );
                else
                    BOOST_CHECK_EQUAL( G.nrEdges(), (N1 <= 2 ? (N1-1) : N1) * N2 * N3 + N1 * (N2 <= 2 ? (N2-1) : N2) * N3 + N1 * N2 * (N3 <= 2 ? (N3-1) : N3) );
                BOOST_CHECK( G.isConnected() );
                BOOST_CHECK_EQUAL( G.isTree(), (G.nrNodes() <= 2) );
                for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
                    foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                        BOOST_CHECK( G.hasEdge( n1, n2 ) );
                        BOOST_CHECK( G.hasEdge( n2, n1 ) );
                    }
                    for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                        if( G.hasEdge( n1, n2 ) ) {
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                        }
                }
                G.checkConsistency();
            }

    // createGraphLoop
    for( size_t N = 0; N < 100; N++ ) {
        GraphAL G = createGraphLoop( N );
        BOOST_CHECK_EQUAL( G.nrNodes(), N );
        if( N == 0 )
            BOOST_CHECK_EQUAL( G.nrEdges(), 0 );
        else if( N <= 2 )
            BOOST_CHECK_EQUAL( G.nrEdges(), N-1 );
        else
            BOOST_CHECK_EQUAL( G.nrEdges(), N );
        BOOST_CHECK( G.isConnected() );
        BOOST_CHECK_EQUAL( G.isTree(), N <= 2 );
        for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
            foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                BOOST_CHECK( G.hasEdge( n1, n2 ) );
                BOOST_CHECK( G.hasEdge( n2, n1 ) );
            }
            for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                if( G.hasEdge( n1, n2 ) ) {
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                }
        }
        G.checkConsistency();
    }

    // createGraphTree
    for( size_t N = 0; N < 100; N++ ) {
        GraphAL G = createGraphTree( N );
        BOOST_CHECK_EQUAL( G.nrNodes(), N );
        BOOST_CHECK_EQUAL( G.nrEdges(), N > 0 ? N - 1 : N );
        BOOST_CHECK( G.isConnected() );
        BOOST_CHECK( G.isTree() );
        for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
            foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                BOOST_CHECK( G.hasEdge( n1, n2 ) );
                BOOST_CHECK( G.hasEdge( n2, n1 ) );
            }
            for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                if( G.hasEdge( n1, n2 ) ) {
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                    BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                    BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                }
        }
        G.checkConsistency();
    }

    // createGraphRegular
    for( size_t N = 0; N < 100; N++ ) {
        for( size_t d = 0; d < N && d <= 20; d++ ) {
            if( (N * d) % 2 == 0 ) {
                GraphAL G = createGraphRegular( N, d );
                BOOST_CHECK_EQUAL( G.nrNodes(), N );
                BOOST_CHECK_EQUAL( G.nrEdges(), d * N / 2 );
                for( size_t n1 = 0; n1 < G.nrNodes(); n1++ ) {
                    BOOST_CHECK_EQUAL( G.nb(n1).size(), d );
                    foreach( const GraphAL::Neighbor &n2, G.nb(n1) ) {
                        BOOST_CHECK( G.hasEdge( n1, n2 ) );
                        BOOST_CHECK( G.hasEdge( n2, n1 ) );
                    }
                    for( size_t n2 = 0; n2 < G.nrNodes(); n2++ )
                        if( G.hasEdge( n1, n2 ) ) {
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ), n2 );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ), n1 );
                            BOOST_CHECK_EQUAL( G.nb( n1, G.findNb( n1, n2 ) ).iter, G.findNb( n1, n2 ) );
                            BOOST_CHECK_EQUAL( G.nb( n2, G.findNb( n2, n1 ) ).iter, G.findNb( n2, n1 ) );
                        }
                }
                G.checkConsistency();
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( StreamTest ) {
    // check printDot
    GraphAL G( 4 );
    G.addEdge( 0, 1 );
    G.addEdge( 0, 2 );
    G.addEdge( 1, 3 );
    G.addEdge( 2, 3 );
    G.addEdge( 2, 2 );
    G.addEdge( 3, 2 );

    std::stringstream ss;
    G.printDot( ss );

    std::string s;
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "graph G {" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "node[shape=circle,width=0.4,fixedsize=true];" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx0;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx1;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx2;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx3;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx0 -- x1;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx0 -- x2;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx1 -- x3;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx2 -- x3;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "}" );
}
