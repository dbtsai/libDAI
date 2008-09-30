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
