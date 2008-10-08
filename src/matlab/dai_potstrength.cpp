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
 *   N = dai_potstrength(psi,i,j);                                 *
 *                                                                 * 
 *=================================================================*/


#include <iostream>
#include "mex.h"
#include <dai/matlab/matlab.h>
#include <dai/factor.h>


using namespace std;
using namespace dai;


/* Input Arguments */

#define PSI_IN          prhs[0]
#define I_IN            prhs[1]
#define J_IN            prhs[2]
#define NR_IN           3


/* Output Arguments */

#define N_OUT           plhs[0]
#define NR_OUT          1


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    long ilabel, jlabel;

    // Check for proper number of arguments
    if ((nrhs != NR_IN) || (nlhs != NR_OUT)) { 
        mexErrMsgTxt("Usage: N = dai_potstrength(psi,i,j);\n\n"
        "\n"
        "INPUT:  psi        = structure with a Member field and a P field, like a CPTAB.\n"
        "        i          = label of a variable in psi.\n"
        "        j          = label of another variable in psi.\n"
        "\n"
        "OUTPUT: N          = strength of psi in direction i->j.\n");
    } 
    
    // Get input parameters
    Factor psi = mx2Factor(PSI_IN);
    ilabel = (long)*mxGetPr(I_IN);
    jlabel = (long)*mxGetPr(J_IN);

    // Find variable in psi with label ilabel
    Var i;
    for( VarSet::const_iterator n = psi.vars().begin(); n != psi.vars().end(); n++ )
        if( n->label() == ilabel ) {
            i = *n;
            break;
        }
    assert( i.label() == ilabel );

    // Find variable in psi with label jlabel
    Var j;
    for( VarSet::const_iterator n = psi.vars().begin(); n != psi.vars().end(); n++ )
        if( n->label() == jlabel ) {
            j = *n;
            break;
        }
    assert( j.label() == jlabel );

    // Calculate N(psi,i,j);
    double N = psi.strength( i, j );
    
    // Hand over result to MATLAB
    N_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(N_OUT)) = N;
    
    return;
}
