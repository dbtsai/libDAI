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


/*=================================================================*
 *                                                                 * 
 * This is a MEX-file for MATLAB.                                  *
 *                                                                 * 
 *   dai_writefg(psi, filename);                                   *
 *                                                                 * 
 *=================================================================*/


#include <iostream>
#include "mex.h"
#include <dai/matlab/matlab.h>
#include <dai/factorgraph.h>


using namespace std;
using namespace dai;


/* Input Arguments */

#define PSI_IN          prhs[0]
#define FILENAME_IN     prhs[1]
#define NR_IN           2


/* Output Arguments */

#define NR_OUT          0


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    char *filename;

    
    // Check for proper number of arguments
    if ((nrhs != NR_IN) || (nlhs != NR_OUT)) { 
        mexErrMsgTxt("Usage: dai_writefg(psi,filename);\n\n"
        "\n"
        "INPUT:  psi        = linear cell array containing the factors\n"
        "                     (psi{i} should be a structure with a Member field\n"
        "                     and a P field, like a CPTAB).\n"
        "        filename   = filename of a .fg file\n");
    } 
    
    // Get input parameters
    vector<Factor> factors = mx2Factors(PSI_IN,0);
    
    size_t buflen;
    buflen = mxGetN( FILENAME_IN ) + 1;
    filename = (char *)mxCalloc( buflen, sizeof(char) );
    mxGetString( FILENAME_IN, filename, buflen );
    
    // Construct factorgraph
    FactorGraph fg(factors);
    long nr_v = fg.nrVars();
    long nr_f = fg.nrFactors();

    try {
        fg.WriteToFile( filename );
    } catch( std::exception &e ) {
        mexErrMsgTxt( e.what() );
    }

    return;
}
