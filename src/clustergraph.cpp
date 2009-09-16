/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <set>
#include <vector>
#include <iostream>
#include <dai/varset.h>
#include <dai/clustergraph.h>


namespace dai {


using namespace std;


ClusterGraph::ClusterGraph( const std::vector<VarSet> & cls ) : G(), vars(), clusters() {
    // construct vars, clusters and edge list
    vector<Edge> edges;
    foreach( const VarSet &cl, cls ) {
        if( find( clusters.begin(), clusters.end(), cl ) == clusters.end() ) {
            // add cluster
            size_t n2 = clusters.size();
            clusters.push_back( cl );
            for( VarSet::const_iterator n = cl.begin(); n != cl.end(); n++ ) {
                size_t n1 = find( vars.begin(), vars.end(), *n ) - vars.begin();
                if( n1 == vars.size() )
                    // add variable
                    vars.push_back( *n );
                edges.push_back( Edge( n1, n2 ) );
            }
        } // disregard duplicate clusters
    }

    // Create bipartite graph
    G.construct( vars.size(), clusters.size(), edges.begin(), edges.end() );
}


ClusterGraph ClusterGraph::VarElim_MinFill() const {
    // Make a copy
    ClusterGraph cl(*this);
    cl.eraseNonMaximal();

    ClusterGraph result;

    // Construct set of variable indices
    set<size_t> varindices;
    for( size_t i = 0; i < vars.size(); ++i )
        varindices.insert( i );

    // Do variable elimination
    while( !varindices.empty() ) {
        set<size_t>::const_iterator lowest = varindices.end();
        size_t lowest_cost = -1UL;
        for( set<size_t>::const_iterator i = varindices.begin(); i != varindices.end(); i++ ) {
            size_t cost = cl.eliminationCost( *i );
            if( lowest == varindices.end() || lowest_cost > cost ) {
                lowest = i;
                lowest_cost = cost;
            }
        }
        size_t i = *lowest;

        result.insert( cl.Delta( i ) );

        cl.insert( cl.delta( i ) );
        cl.eraseSubsuming( i );
        cl.eraseNonMaximal();
        varindices.erase( i );
    }

    return result;
}



ClusterGraph ClusterGraph::VarElim( const std::vector<Var> & ElimSeq ) const {
    // Make a copy
    ClusterGraph cl(*this);
    cl.eraseNonMaximal();

    ClusterGraph result;

    // Do variable elimination
    for( vector<Var>::const_iterator n = ElimSeq.begin(); n != ElimSeq.end(); n++ ) {
        size_t i = cl.findVar( *n );
        assert( i != cl.vars.size() );

        result.insert( cl.Delta(i) );

        cl.insert( cl.delta(i) );
        cl.eraseSubsuming( i );
        cl.eraseNonMaximal();
    }

    return result;
}


} // end of namespace dai
