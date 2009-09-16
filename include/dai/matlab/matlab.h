/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines some utility functions for interfacing with MatLab


#ifndef __defined_libdai_matlab_h
#define __defined_libdai_matlab_h


#include "mex.h"
#include <dai/factor.h>


namespace dai {


#ifdef SMALLMEM
    typedef int mwSize;
    typedef int mwIndex;
#endif


/// Convert vector<Factor> structure to a cell vector of CPTAB-like structs
mxArray *Factors2mx(const std::vector<Factor> &Ps);

/// Convert cell vector of CPTAB-like structs to vector<Factor>
std::vector<Factor> mx2Factors(const mxArray *psi, long verbose);

/// Convert CPTAB-like struct to Factor
Factor mx2Factor(const mxArray *psi);


} // end of namespace dai


#endif
