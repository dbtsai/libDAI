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

        // Check whether the root is in the tree
        bool valid = false;
        for( GraphEL::iterator e = Gr.begin(); e != Gr.end() && !valid; e++ )
            if( e->n1 == Root || e->n2 == Root )
                valid = true;
        if( !valid )
            DAI_THROWE(RUNTIME_ERROR,"Graph does not contain specified root.");

        // Start with the root
        nodes.insert( Root );

        // Keep adding edges until done
        bool done = false;
        while( !done ) {
            bool changed = false;
            for( GraphEL::iterator e = Gr.begin(); e != Gr.end(); ) {
                bool e1_in_nodes = nodes.count( e->n1 );
                bool e2_in_nodes = nodes.count( e->n2 );
                if( e1_in_nodes && e2_in_nodes )
                    DAI_THROWE(RUNTIME_ERROR,"Graph is not acyclic.");
                if( e1_in_nodes ) {
                    // Add directed edge, pointing away from the root
                    push_back( DEdge( e->n1, e->n2 ) );
                    nodes.insert( e->n2 );
                    // Erase the edge
                    Gr.erase( e++ );
                    changed = true;
                } else if( e2_in_nodes ) {
                    // Add directed edge, pointing away from the root
                    push_back( DEdge( e->n2, e->n1 ) );
                    nodes.insert( e->n1 );
                    // Erase the edge
                    Gr.erase( e++ );
                    changed = true;
                } else
                    e++;
            }
            if( Gr.empty() )
                done = true;
            if( !changed && !done )
                DAI_THROWE(RUNTIME_ERROR,"Graph is not connected.");
        }
    }
}


} // end of namespace dai
