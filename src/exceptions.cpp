/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <dai/exceptions.h>


namespace dai {


    std::string Exception::ErrorStrings[NUM_ERRORS] = {
        "Assertion failed",
        "Feature not implemented",
        "Unknown DAI algorithm",
        "Unknown Property type",
        "Malformed Property",
        "Unknown ENUM value",
        "Cannot read file",
        "Cannot write file",
        "Invalid FactorGraph file",
        "Not all mandatory Properties specified",
        "Multiple undo levels unsupported",
        "FactorGraph is not connected",
        "Impossible typecast",
        "Internal error",
        "Runtime error",
        "Quantity not normalizable",
        "Invalid Evidence file",
        "Invalid Expectation-Maximization file",
        "Unrecognized parameter estimation method",
        "Requested object not found"
    };


}
