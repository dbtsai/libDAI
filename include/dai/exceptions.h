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


/// \file
/// \brief Defines Exception class and the DAI_THROW macro


#ifndef __defined_libdai_exceptions_h
#define __defined_libdai_exceptions_h


#include <exception>
#include <stdexcept>
#include <string>


/// Used by DAI_THROW
#define DAI_QUOTE(x) #x

/// Used by DAI_THROW
#define DAI_TOSTRING(x) DAI_QUOTE(x)

/// Macro that simplifies throwing an exceptions with a useful error message
/** \param cod Corresponds to one of the enum values in dai::Exception::codes
 *
 * Example:
 *  \code
 *  DAI_THROW(NOT_IMPLEMENTED);
 *  \endcode
 */
#define DAI_THROW(cod) throw dai::Exception(dai::Exception::cod, std::string(__FILE__ ", line " DAI_TOSTRING(__LINE__)))


namespace dai {


/// Represents an exception (based on std::runtime_error)
class Exception : public std::runtime_error {
    public:
        /// Constructor
            Exception(size_t code, const std::string& msg = "") : std::runtime_error(ErrorStrings[code] + " [" +  msg + "]") {}

        /// Enumeration of exceptions used in libDAI
        enum codes {NOT_IMPLEMENTED,
                    UNKNOWN_DAI_ALGORITHM,
                    UNKNOWN_PROPERTY_TYPE,
                    MALFORMED_PROPERTY,
                    UNKNOWN_ENUM_VALUE,
                    CANNOT_READ_FILE,
                    CANNOT_WRITE_FILE,
                    INVALID_FACTORGRAPH_FILE,
                    NOT_ALL_PROPERTIES_SPECIFIED,
                    MULTIPLE_UNDO,
                    FACTORGRAPH_NOT_CONNECTED,
                    IMPOSSIBLE_TYPECAST,
                    INTERNAL_ERROR,
                    NUM_ERRORS};  // NUM_ERRORS should be the last entry

    private:
        /// Error messages corresponding to the exceptions enumerated above
        static std::string ErrorStrings[NUM_ERRORS];
};


}


#endif
