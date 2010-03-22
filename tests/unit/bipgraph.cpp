/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/bipgraph.h>
#include <vector>
#include <strstream>


using namespace dai;


#define BOOST_TEST_MODULE BipartiteGraphTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    typedef BipartiteGraph::Edge Edge;

    BipartiteGraph G;
    BOOST_CHECK_EQUAL( G.nrNodes1(), 0 );
    BOOST_CHECK_EQUAL( G.nrNodes2(), 0 );
    BOOST_CHECK_EQUAL( G.nrEdges(), 0 );
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.checkConsistency();

    std::vector<Edge> edges;
    edges.push_back( Edge(0, 0) );
    edges.push_back( Edge(0, 1) );
    edges.push_back( Edge(1, 1) );
    edges.push_back( Edge(1, 2) );
    edges.push_back( Edge(1, 2) );
    BipartiteGraph G2( 2, 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G2.nrNodes1(), 2 );
    BOOST_CHECK_EQUAL( G2.nrNodes2(), 3 );
    BOOST_CHECK_EQUAL( G2.nrEdges(), 4 );
    BOOST_CHECK( G2.isConnected() );
    BOOST_CHECK( G2.isTree() );
    G2.checkConsistency();

    edges.push_back( Edge(1, 0) );
    BipartiteGraph G3( 2, 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G3.nrNodes1(), 2 );
    BOOST_CHECK_EQUAL( G3.nrNodes2(), 3 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 5 );
    BOOST_CHECK( G3.isConnected() );
    BOOST_CHECK( !G3.isTree() );
    G3.checkConsistency();

    BipartiteGraph G4( 3, 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G4.nrNodes1(), 3 );
    BOOST_CHECK_EQUAL( G4.nrNodes2(), 3 );
    BOOST_CHECK_EQUAL( G4.nrEdges(), 5 );
    BOOST_CHECK( !G4.isConnected() );
    BOOST_CHECK( !G4.isTree() );
    G4.checkConsistency();

    G.construct( 3, 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G.nrNodes1(), 3 );
    BOOST_CHECK_EQUAL( G.nrNodes2(), 3 );
    BOOST_CHECK_EQUAL( G.nrEdges(), 5 );
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.checkConsistency();
}


BOOST_AUTO_TEST_CASE( NeighborTest ) {
    // check nb() accessor / mutator
    typedef BipartiteGraph::Edge Edge;
    std::vector<Edge> edges;
    edges.push_back( Edge(0, 0) );
    edges.push_back( Edge(0, 1) );
    edges.push_back( Edge(1, 1) );
    edges.push_back( Edge(1, 2) );
    BipartiteGraph G( 2, 3, edges.begin(), edges.end() );
    BOOST_CHECK_EQUAL( G.nb1(0).size(), 2 );
    BOOST_CHECK_EQUAL( G.nb1(1).size(), 2 );
    BOOST_CHECK_EQUAL( G.nb2(0).size(), 1 );
    BOOST_CHECK_EQUAL( G.nb2(1).size(), 2 );
    BOOST_CHECK_EQUAL( G.nb2(2).size(), 1 );
    BOOST_CHECK_EQUAL( G.nb1(0,0).iter, 0 );
    BOOST_CHECK_EQUAL( G.nb1(0,0).node, 0 );
    BOOST_CHECK_EQUAL( G.nb1(0,0).dual, 0 );
    BOOST_CHECK_EQUAL( G.nb1(0,1).iter, 1 );
    BOOST_CHECK_EQUAL( G.nb1(0,1).node, 1 );
    BOOST_CHECK_EQUAL( G.nb1(0,1).dual, 0 );
    BOOST_CHECK_EQUAL( G.nb1(1,0).iter, 0 );
    BOOST_CHECK_EQUAL( G.nb1(1,0).node, 1 );
    BOOST_CHECK_EQUAL( G.nb1(1,0).dual, 1 );
    BOOST_CHECK_EQUAL( G.nb1(1,1).iter, 1 );
    BOOST_CHECK_EQUAL( G.nb1(1,1).node, 2 );
    BOOST_CHECK_EQUAL( G.nb1(1,1).dual, 0 );
    BOOST_CHECK_EQUAL( G.nb2(0,0).iter, 0 );
    BOOST_CHECK_EQUAL( G.nb2(0,0).node, 0 );
    BOOST_CHECK_EQUAL( G.nb2(0,0).dual, 0 );
    BOOST_CHECK_EQUAL( G.nb2(1,0).iter, 0 );
    BOOST_CHECK_EQUAL( G.nb2(1,0).node, 0 );
    BOOST_CHECK_EQUAL( G.nb2(1,0).dual, 1 );
    BOOST_CHECK_EQUAL( G.nb2(1,1).iter, 1 );
    BOOST_CHECK_EQUAL( G.nb2(1,1).node, 1 );
    BOOST_CHECK_EQUAL( G.nb2(1,1).dual, 0 );
    BOOST_CHECK_EQUAL( G.nb2(2,0).iter, 0 );
    BOOST_CHECK_EQUAL( G.nb2(2,0).node, 1 );
    BOOST_CHECK_EQUAL( G.nb2(2,0).dual, 1 );
}


BOOST_AUTO_TEST_CASE( AddEraseTest ) {
    // check addition and erasure of nodes and edges
    typedef BipartiteGraph::Edge Edge;
    std::vector<Edge> edges;
    edges.push_back( Edge( 0, 0 ) );
    edges.push_back( Edge( 0, 1 ) );
    edges.push_back( Edge( 1, 1 ) );
    BipartiteGraph G( 2, 2, edges.begin(), edges.end() );
    G.checkConsistency();
    BOOST_CHECK_EQUAL( G.nrNodes1(), 2 );
    BOOST_CHECK_EQUAL( G.nrNodes2(), 2 );
    BOOST_CHECK_EQUAL( G.nrEdges(), 3 );
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    BOOST_CHECK_EQUAL( G.addNode1(), 2 );
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    BOOST_CHECK_EQUAL( G.addNode2(), 2 );
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    BOOST_CHECK_EQUAL( G.addNode1(), 3 );
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.checkConsistency();
    std::vector<size_t> nbs;
    nbs.push_back( 2 );
    nbs.push_back( 0 );
    BOOST_CHECK_EQUAL( G.addNode1( nbs.begin(), nbs.end() ), 4 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    BOOST_CHECK_EQUAL( G.addNode2( nbs.begin(), nbs.end() ), 3 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.addEdge( 3, 3 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.addEdge( 1, 3 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    BOOST_CHECK_EQUAL( G.nrNodes1(), 5 );
    BOOST_CHECK_EQUAL( G.nrNodes2(), 4 );
    BOOST_CHECK_EQUAL( G.nrEdges(), 9 );
    G.eraseEdge( 0, 3 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseEdge( 4, 2 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.eraseNode2( 2 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode1( 0 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.addEdge( 1, 0 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode1( 2 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode2( 2 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.addEdge( 1, 1 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode2( 1 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.eraseNode1( 1 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.addEdge( 0, 0 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode2( 0 );
    G.checkConsistency();
    BOOST_CHECK( !G.isConnected() );
    BOOST_CHECK( !G.isTree() );
    G.eraseNode1( 0 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    G.eraseNode1( 0 );
    G.checkConsistency();
    BOOST_CHECK( G.isConnected() );
    BOOST_CHECK( G.isTree() );
    BOOST_CHECK_EQUAL( G.nrNodes1(), 0 );
    BOOST_CHECK_EQUAL( G.nrNodes2(), 0 );
    BOOST_CHECK_EQUAL( G.nrEdges(), 0 );
}


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    // check queries which have not been tested in another test case
    typedef BipartiteGraph::Edge Edge;
    std::vector<Edge> edges;
    edges.push_back( Edge( 0, 0 ) );
    edges.push_back( Edge( 0, 1 ) );
    edges.push_back( Edge( 1, 1 ) );
    BipartiteGraph G( 3, 2, edges.begin(), edges.end() );
    G.checkConsistency();
    std::vector<size_t> v;
    std::vector<size_t> v0; v0.push_back(0);
    std::vector<size_t> v1; v1.push_back(1);
    std::vector<size_t> v01; v01.push_back(0); v01.push_back(1);
    BOOST_CHECK( G.delta1( 0, true ) == v01 );
    BOOST_CHECK( G.delta1( 1, true ) == v01 );
    BOOST_CHECK( G.delta1( 2, true ) == v );
    BOOST_CHECK( G.delta2( 0, true ) == v01 );
    BOOST_CHECK( G.delta2( 1, true ) == v01 );
    BOOST_CHECK( G.delta1( 0, false ) == v1 );
    BOOST_CHECK( G.delta1( 1, false ) == v0 );
    BOOST_CHECK( G.delta1( 2, false ) == v );
    BOOST_CHECK( G.delta2( 0, false ) == v1 );
    BOOST_CHECK( G.delta2( 1, false ) == v0 );
}


BOOST_AUTO_TEST_CASE( StreamTest ) {
    // check printDot
    typedef BipartiteGraph::Edge Edge;
    std::vector<Edge> edges;
    edges.push_back( Edge(0, 0) );
    edges.push_back( Edge(0, 1) );
    edges.push_back( Edge(1, 1) );
    edges.push_back( Edge(1, 2) );
    BipartiteGraph G( 2, 3, edges.begin(), edges.end() );

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
    BOOST_CHECK_EQUAL( s, "node[shape=box,width=0.3,height=0.3,fixedsize=true];" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\ty0;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\ty1;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\ty2;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx0 -- y0;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx0 -- y1;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx1 -- y1;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "\tx1 -- y2;" );
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "}" );
}
