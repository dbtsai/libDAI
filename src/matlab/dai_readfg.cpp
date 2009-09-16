/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <iostream>
#include "mex.h"
#include <dai/matlab/matlab.h>
#include <dai/factorgraph.h>


using namespace std;
using namespace dai;


/* Input Arguments */

#define FILENAME_IN     prhs[0]
#define NR_IN           1


/* Output Arguments */

#define PSI_OUT         plhs[0]
#define NR_OUT          1


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) {
    char *filename;


    // Check for proper number of arguments
    if ((nrhs != NR_IN) || (nlhs != NR_OUT)) {
        mexErrMsgTxt("Usage: [psi] = dai_readfg(filename);\n\n"
        "\n"
        "INPUT:  filename   = filename of a .fg file\n"
        "\n"
        "OUTPUT: psi        = linear cell array containing the factors\n"
        "                     (psi{i} is a structure with a Member field\n"
        "                     and a P field, like a CPTAB).\n");
    }

    // Get input parameters
    size_t buflen;
    buflen = mxGetN( FILENAME_IN ) + 1;
    filename = (char *)mxCalloc( buflen, sizeof(char) );
    mxGetString( FILENAME_IN, filename, buflen );


    // Read factorgraph
    FactorGraph fg;
    try {
        fg.ReadFromFile( filename );
    } catch( std::exception &e ) {
        mexErrMsgTxt( e.what() );
    }


    // Save factors
    vector<Factor> psi;
    for( size_t I = 0; I < fg.nrFactors(); I++ )
        psi.push_back(fg.factor(I));


    // Hand over results to MATLAB
    PSI_OUT = Factors2mx(psi);


    return;
}
