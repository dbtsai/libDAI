/*  Copyright (C) 2009  Charles Vaske  [cvaske at soe dot ucsc dot edu]
    University of California Santa Cruz

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
#include <string>

#include <dai/factorgraph.h>
#include <dai/evidence.h>
#include <dai/alldai.h>


using namespace std;
using namespace dai;


void usage( const string &msg ) {
    cerr << msg << endl;
    cerr << "Usage:" << endl;
    cerr << " testem factorgraph.fg evidence.tab emconfig.em" << endl;
    exit( 1 );
}


int main( int argc, char** argv ) {
    if( argc != 4 )
        usage("Incorrect number of arguments.");
    
    FactorGraph fg;
    fg.ReadFromFile( argv[1] );

    PropertySet infprops;
    infprops.Set( "verbose", (size_t)1 );
    infprops.Set( "updates", string("HUGIN") );
    InfAlg* inf = newInfAlg( "JTREE", fg, infprops );
    inf->init();

    Evidence e;
    ifstream estream( argv[2] );
    e.addEvidenceTabFile( estream, fg );

    cout << "Number of samples: " << e.nrSamples() << endl;
    for( Evidence::iterator ps = e.begin(); ps != e.end(); ps++ )
        cout << "Sample #" << (ps - e.begin()) << " has " << ps->observations().size() << " observations." << endl;

    ifstream emstream( argv[3] );
    EMAlg em(e, *inf, emstream);

    while( !em.hasSatisfiedTermConditions() ) {
        Real l = em.iterate();
        cout << "Iteration " << em.getCurrentIters() << " likelihood: " << l <<endl;
    }

    cout << endl << "Inferred Factor Graph:" << endl << "######################" << endl;
    cout.precision(12);
    cout << inf->fg();

    delete inf;

    return 0;
}
