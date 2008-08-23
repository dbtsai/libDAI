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


/*=================================================================*
 *                                                                 * 
 * This is a MEX-file for MATLAB.                                  *
 *                                                                 * 
 *   [psi_out] = dai_removeshortloops(psi_in);                     *
 *                                                                 * 
 *=================================================================*/


#include <iostream>
#include "mex.h"
#include "matlab.h"
#include "factorgraph.h"


using namespace std;


/* Input Arguments */

#define PSI_IN          prhs[0]
#define NR_IN           1


/* Output Arguments */

#define PSI_OUT         plhs[0]
#define NR_OUT          1


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) { 
    // Check for proper number of arguments
    if ((nrhs != NR_IN) || (nlhs != NR_OUT)) { 
        mexErrMsgTxt("Usage: [psi_out] = dai_removeshortloops(psi_in);\n\n"
        "\n"
        "INPUT:  psi_in     = linear cell array containing the factors\n"
        "                     (psi{i} is a structure with a Member field\n"
        "                     and a P field, like a CPTAB).\n"
        "\n"
        "OUTPUT: psi_out    = linear cell array containing the factors of psi_in,\n"
        "                     where factors have been merged to prevent short loops\n"
        "                     of length 4 in the factor graph (i.e. loops of type\n"
        "                     var1-factor1-var2-factor2-var1),\n");
    } 
    
    // Get the factors from PSI_IN
    vector<Factor> psi = mx2Factors(PSI_IN, 0);

    // Remove the short loops
    RemoveShortLoops(psi);

    // Hand over results to MATLAB
    PSI_OUT = Factors2mx(psi);

    
    return;
}
