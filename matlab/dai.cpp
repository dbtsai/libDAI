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
 *   [logZ,q,md] = dai(psi,method,opts);                           *
 *                                                                 * 
 *=================================================================*/


#include <iostream>
#include <dai/matlab/matlab.h>
#include "mex.h"
#include <dai/alldai.h>


using namespace std;
using namespace dai;


/* Input Arguments */

#define PSI_IN          prhs[0]
#define METHOD_IN       prhs[1]
#define OPTS_IN         prhs[2]
#define NR_IN           3


/* Output Arguments */

#define LOGZ_OUT        plhs[0]
#define Q_OUT           plhs[1]
#define MD_OUT          plhs[2]
#define NR_OUT          3


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    size_t buflen;
    
    /* Check for proper number of arguments */
    if( (nrhs != NR_IN) || (nlhs != NR_OUT) ) { 
        mexErrMsgTxt("Usage: [logZ,q,md] = dai(psi,method,opts)\n\n"
        "\n"
        "INPUT:  psi        = linear cell array containing the factors \n"
        "                     psi{i} should be a structure with a Member field\n"
        "                     and a P field, like a CPTAB).\n"
        "        method     = name of the method (see README)\n"
        "        opts       = string of options (see README)\n"
        "\n"
        "OUTPUT: logZ       = approximation of the logarithm of the partition sum.\n"
        "        q          = linear cell array containing all final beliefs.\n"
        "        md         = maxdiff (final linf-dist between new and old single node beliefs).\n");
    } 
    
    char *method;
    char *opts;


    // Get psi and construct factorgraph
    vector<Factor> factors = mx2Factors(PSI_IN, 0);
    FactorGraph fg(factors);
    long nr_v = fg.nrVars();

    // Get method
    buflen = mxGetN( METHOD_IN ) + 1;
    method = (char *)mxCalloc( buflen, sizeof(char) );
    mxGetString( METHOD_IN, method, buflen );

    // Get options string
    buflen = mxGetN( OPTS_IN ) + 1;
    opts = (char *)mxCalloc( buflen, sizeof(char) );
    mxGetString( OPTS_IN, opts, buflen );
    // Convert to options object props
    stringstream ss;
    ss << opts;
    Properties props;
    ss >> props;
    
    // Construct InfAlg object, init and run
    InfAlg *obj = newInfAlg( method, fg, props );
    obj->init();
    obj->run();


    // Save logZ
    double logZ = obj->logZ();

    // Save maxdiff
    double maxdiff = obj->MaxDiff();


    // Hand over results to MATLAB
    LOGZ_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(LOGZ_OUT)) = logZ;
    
    Q_OUT = Factors2mx(obj->beliefs());
    
    MD_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(MD_OUT)) = maxdiff;

    
    return;
}
