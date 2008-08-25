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


#include "factorgraph.h"
#include <iostream>
#include <cstdlib>

using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        cout << "Usage: " << argv[0] << " <in.fg> <out.fg>" << endl << endl;
        cout << "Merges short loops (of length 4 in the factor graph) in <in.fg>" << endl;
        cout << "and writes result to <out.fg>." << endl;
        cout << "Returns 2 in case of error reading <in.fg>, 1 if there are no" << endl;
        cout << "short loops present and 0 if short loops have been merged." << endl;
        return 1;
    } else {
        // Read factorgraph
        FactorGraph fg;
        char *infile = argv[1];

        if( fg.ReadFromFile( infile ) ) {
            cerr << "Error reading file " << infile << endl;
            return 2;
        } else {
            vector<Factor> fg_facs;

            fg_facs.reserve( fg.nrFactors() );
            for( size_t I = 0; I < fg.nrFactors(); I++ )
                fg_facs.push_back( fg.factor(I) );

            if( hasShortLoops( fg_facs) ) {
                RemoveShortLoops( fg_facs );
                FactorGraph newfg( fg_facs );
                newfg.WriteToFile( argv[2] );
                return 0;
            } else {
                fg.WriteToFile( argv[2] );
                return 1;
            }
        }
    }
}
