/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

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


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <dai/factorgraph.h>


using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        cout << "Usage: " << argv[0] << " <in.fg> <out.dot>" << endl << endl;
        cout << "Converts a .fg (FactorGraph) file to a .dot (GraphViz) file for" << endl;
        cout << "visualization. The .dot file can be converted to .ps (PostScript) by" << endl;
        cout << "'neato -T ps out.dot > out.ps' or by 'dot -T ps out.dot > out.ps'" << endl << endl;
        return 1;
    } else {
        // Read factorgraph
        FactorGraph fg;
        char *infile = argv[1];

        fg.ReadFromFile( infile );

        ostream *os = &cout;
        ofstream outfile;
        if( string( argv[2] ) != "-" ) {
            outfile.open( argv[2] );
            if( !outfile.is_open() ) {
                cerr << "Cannot open " << argv[2] << " for writing" << endl;
                return 1;
            }
            os = &outfile;
        }

        fg.printDot( *os );

        return 0;
    }
}
