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
    // TODO   
}


BOOST_AUTO_TEST_CASE( AddEraseTest ) {
    // TODO   
}


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    // TODO   
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
