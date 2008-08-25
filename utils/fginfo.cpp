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

using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if( argc != 2 ) {
        cout << "Usage: " << argv[0] << " <in.fg>" << endl << endl;
        cout << "Reports some characteristics of the .fg network." << endl;
        return 1;
    } else {
        // Read factorgraph
        FactorGraph fg;
        char *infile = argv[1];

        if( fg.ReadFromFile( infile ) ) {
            cerr << "Error reading file " << infile << endl;
            return 2;
        } else {
            cout << "Number of variables: " << fg.nrVars() << endl;
            cout << "Number of factors:   " << fg.nrFactors() << endl;
            cout << "Connected:           " << fg.isConnected() << endl;
//          cout << "Treewidth:           " << endl;

            double cavsum_lcbp = 0.0;
            double cavsum_lcbp2 = 0.0;
            size_t max_Delta_size = 0;
            for( size_t i = 0; i < fg.nrVars(); i++ ) {
                VarSet di = fg.delta(fg.var(i));
                size_t Ds = fg.Delta(fg.var(i)).stateSpace();
                if( Ds > max_Delta_size )
                    max_Delta_size = Ds;
                cavsum_lcbp += di.stateSpace();
                for( VarSet::const_iterator j = di.begin(); j != di.end(); j++ )
                    cavsum_lcbp2 += j->states();
            }
            cout << "Maximum pancake has " << max_Delta_size << " states" << endl;
            cout << "LCBP with full cavities needs " << cavsum_lcbp << " BP runs" << endl;
            cout << "LCBP with only pairinteractions needs " << cavsum_lcbp2 << " BP runs" << endl;

            return 0;
        }
    }
}

