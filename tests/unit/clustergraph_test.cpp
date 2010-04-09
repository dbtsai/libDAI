/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#include <dai/clustergraph.h>
#include <vector>
#include <strstream>


using namespace dai;


const double tol = 1e-8;


#define BOOST_TEST_MODULE ClusterGraphTest


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    ClusterGraph G;
    BOOST_CHECK_EQUAL( G.clusters(), std::vector<VarSet>() );
    BOOST_CHECK( G.bipGraph() == BipartiteGraph() );
    BOOST_CHECK_EQUAL( G.nrVars(), 0 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 0 );
#ifdef DAI_DEBUG
    BOOST_CHECK_THROW( G.var( 0 ), Exception );
    BOOST_CHECK_THROW( G.cluster( 0 ), Exception );
#endif
    BOOST_CHECK_THROW( G.findVar( Var( 0, 2 ) ), Exception );

    Var v0( 0, 2 );
    Var v1( 1, 3 );
    Var v2( 2, 2 );
    Var v3( 3, 4 );
    std::vector<Var> vs;
    vs.push_back( v0 );
    vs.push_back( v1 );
    vs.push_back( v2 );
    vs.push_back( v3 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v03( v0, v3 );
    VarSet v12( v1, v2 );
    VarSet v13( v1, v3 );
    VarSet v23( v2, v3 );
    std::vector<VarSet> cl;
    cl.push_back( v01 );
    cl.push_back( v12 );
    cl.push_back( v23 );
    cl.push_back( v13 );
    ClusterGraph G2( cl );
    BOOST_CHECK_EQUAL( G2.nrVars(), 4 );
    BOOST_CHECK_EQUAL( G2.nrClusters(), 4 );
    BOOST_CHECK_EQUAL( G2.vars(), vs );
    BOOST_CHECK_EQUAL( G2.clusters(), cl );
    BOOST_CHECK_EQUAL( G2.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G2.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G2.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G2.findVar( v3 ), 3 );

    ClusterGraph Gb( G );
    BOOST_CHECK( G.bipGraph() == Gb.bipGraph() );
    BOOST_CHECK( G.vars() == Gb.vars() );
    BOOST_CHECK( G.clusters() == Gb.clusters() );

    ClusterGraph Gc = G;
    BOOST_CHECK( G.bipGraph() == Gc.bipGraph() );
    BOOST_CHECK( G.vars() == Gc.vars() );
    BOOST_CHECK( G.clusters() == Gc.clusters() );

    ClusterGraph G2b( G2 );
    BOOST_CHECK( G2.bipGraph() == G2b.bipGraph() );
    BOOST_CHECK( G2.vars() == G2b.vars() );
    BOOST_CHECK( G2.clusters() == G2b.clusters() );

    ClusterGraph G2c = G2;
    BOOST_CHECK( G2.bipGraph() == G2c.bipGraph() );
    BOOST_CHECK( G2.vars() == G2c.vars() );
    BOOST_CHECK( G2.clusters() == G2c.clusters() );
}


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 3 );
    Var v2( 2, 2 );
    Var v3( 3, 4 );
    Var v4( 4, 2 );
    std::vector<Var> vs;
    vs.push_back( v0 );
    vs.push_back( v1 );
    vs.push_back( v2 );
    vs.push_back( v3 );
    vs.push_back( v4 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v03( v0, v3 );
    VarSet v04( v0, v4 );
    VarSet v12( v1, v2 );
    VarSet v13( v1, v3 );
    VarSet v14( v1, v4 );
    VarSet v23( v2, v3 );
    VarSet v24( v2, v4 );
    VarSet v34( v3, v4 );
    VarSet v123 = v12 | v3;
    std::vector<VarSet> cl;
    cl.push_back( v01 );
    cl.push_back( v12 );
    cl.push_back( v123 );
    cl.push_back( v34 );
    cl.push_back( v04 );
    ClusterGraph G( cl );

    BOOST_CHECK_EQUAL( G.nrVars(), 5 );
    BOOST_CHECK_EQUAL( G.vars(), vs );
    BOOST_CHECK_EQUAL( G.var(0), v0 );
    BOOST_CHECK_EQUAL( G.var(1), v1 );
    BOOST_CHECK_EQUAL( G.var(2), v2 );
    BOOST_CHECK_EQUAL( G.var(3), v3 );
    BOOST_CHECK_EQUAL( G.var(4), v4 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 5 );
    BOOST_CHECK_EQUAL( G.clusters(), cl );
    BOOST_CHECK_EQUAL( G.cluster(0), v01 );
    BOOST_CHECK_EQUAL( G.cluster(1), v12 );
    BOOST_CHECK_EQUAL( G.cluster(2), v123 );
    BOOST_CHECK_EQUAL( G.cluster(3), v34 );
    BOOST_CHECK_EQUAL( G.cluster(4), v04 );
    BOOST_CHECK_EQUAL( G.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G.findVar( v3 ), 3 );
    BOOST_CHECK_EQUAL( G.findVar( v4 ), 4 );
    BipartiteGraph H( 5, 5 );
    H.addEdge( 0, 0 );
    H.addEdge( 1, 0 );
    H.addEdge( 1, 1 );
    H.addEdge( 2, 1 );
    H.addEdge( 1, 2 );
    H.addEdge( 2, 2 );
    H.addEdge( 3, 2 );
    H.addEdge( 3, 3 );
    H.addEdge( 4, 3 );
    H.addEdge( 0, 4 );
    H.addEdge( 4, 4 );
    BOOST_CHECK( G.bipGraph() == H );

    BOOST_CHECK_EQUAL( G.delta( 0 ), v14 );
    BOOST_CHECK_EQUAL( G.delta( 1 ), v02 | v3 );
    BOOST_CHECK_EQUAL( G.delta( 2 ), v13 );
    BOOST_CHECK_EQUAL( G.delta( 3 ), v12 | v4 );
    BOOST_CHECK_EQUAL( G.delta( 4 ), v03 );
    BOOST_CHECK_EQUAL( G.Delta( 0 ), v14 | v0 );
    BOOST_CHECK_EQUAL( G.Delta( 1 ), v01 | v23 );
    BOOST_CHECK_EQUAL( G.Delta( 2 ), v13 | v2 );
    BOOST_CHECK_EQUAL( G.Delta( 3 ), v12 | v34 );
    BOOST_CHECK_EQUAL( G.Delta( 4 ), v03 | v4 );

    BOOST_CHECK( !G.adj( 0, 0 ) );
    BOOST_CHECK(  G.adj( 0, 1 ) );
    BOOST_CHECK( !G.adj( 0, 2 ) );
    BOOST_CHECK( !G.adj( 0, 3 ) );
    BOOST_CHECK(  G.adj( 0, 4 ) );
    BOOST_CHECK(  G.adj( 1, 0 ) );
    BOOST_CHECK( !G.adj( 1, 1 ) );
    BOOST_CHECK(  G.adj( 1, 2 ) );
    BOOST_CHECK(  G.adj( 1, 3 ) );
    BOOST_CHECK( !G.adj( 1, 4 ) );
    BOOST_CHECK( !G.adj( 2, 0 ) );
    BOOST_CHECK(  G.adj( 2, 1 ) );
    BOOST_CHECK( !G.adj( 2, 2 ) );
    BOOST_CHECK(  G.adj( 2, 3 ) );
    BOOST_CHECK( !G.adj( 2, 4 ) );
    BOOST_CHECK( !G.adj( 3, 0 ) );
    BOOST_CHECK(  G.adj( 3, 1 ) );
    BOOST_CHECK(  G.adj( 3, 2 ) );
    BOOST_CHECK( !G.adj( 3, 3 ) );
    BOOST_CHECK(  G.adj( 3, 4 ) );
    BOOST_CHECK(  G.adj( 4, 0 ) );
    BOOST_CHECK( !G.adj( 4, 1 ) );
    BOOST_CHECK( !G.adj( 4, 2 ) );
    BOOST_CHECK(  G.adj( 4, 3 ) );
    BOOST_CHECK( !G.adj( 4, 4 ) );

    BOOST_CHECK(  G.isMaximal( 0 ) );
    BOOST_CHECK( !G.isMaximal( 1 ) );
    BOOST_CHECK(  G.isMaximal( 2 ) );
    BOOST_CHECK(  G.isMaximal( 3 ) );
    BOOST_CHECK(  G.isMaximal( 4 ) );
}


BOOST_AUTO_TEST_CASE( OperationsTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 3 );
    Var v2( 2, 2 );
    Var v3( 3, 4 );
    Var v4( 4, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v03( v0, v3 );
    VarSet v04( v0, v4 );
    VarSet v12( v1, v2 );
    VarSet v13( v1, v3 );
    VarSet v14( v1, v4 );
    VarSet v23( v2, v3 );
    VarSet v24( v2, v4 );
    VarSet v34( v3, v4 );
    VarSet v123 = v12 | v3;
    std::vector<VarSet> cl;
    cl.push_back( v01 );
    cl.push_back( v12 );
    cl.push_back( v123 );
    cl.push_back( v34 );
    cl.push_back( v04 );
    ClusterGraph G( cl );

    BipartiteGraph H( 5, 5 );
    H.addEdge( 0, 0 );
    H.addEdge( 1, 0 );
    H.addEdge( 1, 1 );
    H.addEdge( 2, 1 );
    H.addEdge( 1, 2 );
    H.addEdge( 2, 2 );
    H.addEdge( 3, 2 );
    H.addEdge( 3, 3 );
    H.addEdge( 4, 3 );
    H.addEdge( 0, 4 );
    H.addEdge( 4, 4 );
    BOOST_CHECK( G.bipGraph() == H );

    G.eraseNonMaximal();
    BOOST_CHECK_EQUAL( G.nrClusters(), 4 );
    H.eraseNode2( 1 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseSubsuming( 4 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 2 );
    H.eraseNode2( 2 );
    H.eraseNode2( 2 );
    BOOST_CHECK( G.bipGraph() == H );
    G.insert( v34 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 3 );
    G.insert( v123 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 3 );
    H.addNode2();
    H.addEdge( 3, 2 );
    H.addEdge( 4, 2 );
    BOOST_CHECK( G.bipGraph() == H );
    G.insert( v12 );
    G.insert( v23 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 5 );
    H.addNode2();
    H.addNode2();
    H.addEdge( 1, 3 );
    H.addEdge( 2, 3 );
    H.addEdge( 2, 4 );
    H.addEdge( 3, 4 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseNonMaximal();
    BOOST_CHECK_EQUAL( G.nrClusters(), 3 );
    H.eraseNode2( 3 );
    H.eraseNode2( 3 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseSubsuming( 2 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 2 );
    H.eraseNode2( 1 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseNonMaximal();
    BOOST_CHECK_EQUAL( G.nrClusters(), 2 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseSubsuming( 0 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 1 );
    H.eraseNode2( 0 );
    BOOST_CHECK( G.bipGraph() == H );
    G.eraseSubsuming( 4 );
    BOOST_CHECK_EQUAL( G.nrClusters(), 0 );
    H.eraseNode2( 0 );
    BOOST_CHECK( G.bipGraph() == H );
}


BOOST_AUTO_TEST_CASE( VarElimTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 3 );
    Var v2( 2, 2 );
    Var v3( 3, 4 );
    Var v4( 4, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v03( v0, v3 );
    VarSet v04( v0, v4 );
    VarSet v12( v1, v2 );
    VarSet v13( v1, v3 );
    VarSet v14( v1, v4 );
    VarSet v23( v2, v3 );
    VarSet v24( v2, v4 );
    VarSet v34( v3, v4 );
    VarSet v123 = v12 | v3;
    std::vector<VarSet> cl;
    cl.push_back( v01 );
    cl.push_back( v12 );
    cl.push_back( v123 );
    cl.push_back( v34 );
    cl.push_back( v04 );
    ClusterGraph G( cl );
    ClusterGraph Gorg = G;

    BipartiteGraph H( 5, 5 );
    H.addEdge( 0, 0 );
    H.addEdge( 1, 0 );
    H.addEdge( 1, 1 );
    H.addEdge( 2, 1 );
    H.addEdge( 1, 2 );
    H.addEdge( 2, 2 );
    H.addEdge( 3, 2 );
    H.addEdge( 3, 3 );
    H.addEdge( 4, 3 );
    H.addEdge( 0, 4 );
    H.addEdge( 4, 4 );
    BOOST_CHECK( G.bipGraph() == H );

    BOOST_CHECK_EQUAL( eliminationCost_MinFill( G, 0 ), 1 );
    BOOST_CHECK_EQUAL( eliminationCost_MinFill( G, 1 ), 2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinFill( G, 2 ), 0 );
    BOOST_CHECK_EQUAL( eliminationCost_MinFill( G, 3 ), 2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinFill( G, 4 ), 1 );
    cl.clear();
    cl.push_back( v123 );
    cl.push_back( v01 | v4 );
    cl.push_back( v13 | v4 );
    cl.push_back( v34 );
    cl.push_back( v4 );
    BOOST_CHECK_EQUAL( G.VarElim( greedyVariableElimination( eliminationCost_MinFill ) ).clusters(), cl );

    G = Gorg;
    BOOST_CHECK_EQUAL( eliminationCost_WeightedMinFill( G, 0 ), 2*3 );
    BOOST_CHECK_EQUAL( eliminationCost_WeightedMinFill( G, 1 ), 2*2+2*4 );
    BOOST_CHECK_EQUAL( eliminationCost_WeightedMinFill( G, 2 ), 0 );
    BOOST_CHECK_EQUAL( eliminationCost_WeightedMinFill( G, 3 ), 3*2+2*2 );
    BOOST_CHECK_EQUAL( eliminationCost_WeightedMinFill( G, 4 ), 2*4 );
    cl.clear();
    cl.push_back( v123 );
    cl.push_back( v01 | v4 );
    cl.push_back( v13 | v4 );
    cl.push_back( v34 );
    cl.push_back( v4 );
    BOOST_CHECK_EQUAL( G.VarElim( greedyVariableElimination( eliminationCost_WeightedMinFill ) ).clusters(), cl );

    G = Gorg;
    BOOST_CHECK_EQUAL( eliminationCost_MinNeighbors( G, 0 ), 2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinNeighbors( G, 1 ), 3 );
    BOOST_CHECK_EQUAL( eliminationCost_MinNeighbors( G, 2 ), 2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinNeighbors( G, 3 ), 3 );
    BOOST_CHECK_EQUAL( eliminationCost_MinNeighbors( G, 4 ), 2 );
    cl.clear();
    cl.push_back( v01 | v4 );
    cl.push_back( v123 );
    cl.push_back( v13 | v4 );
    cl.push_back( v34 );
    cl.push_back( v4 );
    BOOST_CHECK_EQUAL( G.VarElim( greedyVariableElimination( eliminationCost_MinNeighbors ) ).clusters(), cl );

    G = Gorg;
    BOOST_CHECK_EQUAL( eliminationCost_MinWeight( G, 0 ), 3*2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinWeight( G, 1 ), 2*2*4 );
    BOOST_CHECK_EQUAL( eliminationCost_MinWeight( G, 2 ), 3*4 );
    BOOST_CHECK_EQUAL( eliminationCost_MinWeight( G, 3 ), 3*2*2 );
    BOOST_CHECK_EQUAL( eliminationCost_MinWeight( G, 4 ), 2*4 );
    cl.clear();
    cl.push_back( v01 | v4 );
    cl.push_back( v123 );
    cl.push_back( v13 | v4 );
    cl.push_back( v14 );
    cl.push_back( v4 );
    BOOST_CHECK_EQUAL( G.VarElim( greedyVariableElimination( eliminationCost_MinWeight ) ).clusters(), cl );

    G = Gorg;
    std::vector<Var> vs;
    vs.push_back( v4 );
    vs.push_back( v3 );
    vs.push_back( v2 );
    vs.push_back( v1 );
    vs.push_back( v0 );
    cl.clear();
    cl.push_back( v03 | v4 );
    cl.push_back( v01 | v23 );
    cl.push_back( v01 | v2 );
    cl.push_back( v01 );
    cl.push_back( v0 );
    BOOST_CHECK_EQUAL( G.VarElim( sequentialVariableElimination( vs ) ).clusters(), cl );
}


BOOST_AUTO_TEST_CASE( IOTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 3 );
    Var v2( 2, 2 );
    Var v3( 3, 4 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v03( v0, v3 );
    VarSet v12( v1, v2 );
    VarSet v13( v1, v3 );
    VarSet v23( v2, v3 );
    std::vector<VarSet> cl;
    cl.push_back( v01 );
    cl.push_back( v12 );
    cl.push_back( v23 );
    cl.push_back( v13 );
    ClusterGraph G( cl );

    std::stringstream ss;
    ss << G;
    std::string s;
    getline( ss, s );
    BOOST_CHECK_EQUAL( s, "({x0, x1}, {x1, x2}, {x2, x3}, {x1, x3})" );
}
