/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <dai/exceptions.h>


namespace dai {


    std::string Exception::ErrorStrings[NUM_ERRORS] = {
        "Feature not implemented",
        "Assertion failed",
        "Impossible typecast",
        "Requested object not found",
        "Requested belief not available",
        "Unknown ENUM value",
        "Unknown DAI algorithm",
        "Unrecognized parameter estimation method",
        "Unknown Property type",
        "Unknown Property",
        "Malformed Property",
        "Not all mandatory Properties specified",
        "Invalid alias",
        "Cannot read file",
        "Cannot write file",
        "Invalid FactorGraph file",
        "Invalid Evidence file",
        "Invalid Expectation-Maximization file",
        "Quantity not normalizable",
        "Multiple undo levels unsupported",
        "FactorGraph is not connected",
        "Internal error",
        "Runtime error"
    };


}
