/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
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


#include "../factorgraph.h"
#include <iostream>
#include <cstdlib>
#include <string>

using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        cout << "Usage: " << argv[0] << " <in.fg> <out.dot>" << endl << endl;
        cout << "Converts a .fg (FactorGraph) file to a .dot (GraphViz) file for visualization." << endl;
        cout << "The .dot file can be converted to .ps (PostScript by 'neato -T ps out.dot > out.ps'" << endl;
        return 1;
    } else {
        // Read factorgraph
        FactorGraph fg;
        char *infile = argv[1];

        if( fg.ReadFromFile( infile ) ) {
            cerr << "Error reading file " << infile << endl;
            return 2;
        } else {
            if( string( argv[2] ) == "-" ) 
                fg.WriteToDotFile( argv[2] );
            else {
                cout << "graph G {" << endl;
                cout << "graph[size=\"9,9\"];" << endl;
                cout << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
                for( size_t i = 0; i < fg.nrVars(); i++ )
                    cout << "\tx" << fg.var(i).label() << ";" << endl;
                cout << "node[shape=box,style=filled,color=lightgrey,width=0.3,height=0.3,fixedsize=true];" << endl;
                for( size_t I = 0; I < fg.nrFactors(); I++ )
                    cout << "\tp" << I << ";" << endl;
                for( size_t iI = 0; iI < fg.nr_edges(); iI++ )
                    cout << "\tx" << fg.var(fg.edge(iI).first).label() << " -- p" << fg.edge(iI).second << ";" << endl;
                cout << "}" << endl;
            }

            return 0;
        }
    }
}

