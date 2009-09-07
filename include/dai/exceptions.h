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
/// \todo Improve documentation


#ifndef __defined_libdai_exceptions_h
#define __defined_libdai_exceptions_h


#include <exception>
#include <stdexcept>
#include <string>
#include <iostream>


/// Used by DAI_THROW
#define DAI_QUOTE(x) #x

/// Used by DAI_THROW
#define DAI_TOSTRING(x) DAI_QUOTE(x)

/// Macro that simplifies throwing an exception with a useful error message.
/** \param cod Corresponds to one of the enum values in dai::Exception::codes
 *
 *  Example:
 *  \code
 *  DAI_THROW(NOT_IMPLEMENTED);
 *  \endcode
 */
#define DAI_THROW(cod) throw dai::Exception(dai::Exception::cod, std::string(__FILE__ ", line " DAI_TOSTRING(__LINE__)))

/// Macro that simplifies throwing an exception with a useful error message. It also allows for writing a detailed error message to stderr.
/** \param cod Corresponds to one of the enum values in dai::Exception::codes
 *  \param msg Detailed error message that will be written to std::cerr.
 *
 *  Example:
 *  \code
 *  DAI_THROWE(NOT_IMPLEMENTED,"Detailed error message");
 *  \endcode
 */
#define DAI_THROWE(cod,msg) throw dai::Exception(dai::Exception::cod, std::string(__FILE__ ", line " DAI_TOSTRING(__LINE__)), msg)


namespace dai {


/// Represents an exception (based on std::runtime_error)
class Exception : public std::runtime_error {
    public:
        /// Enumeration of exceptions used in libDAI
        enum Code {NOT_IMPLEMENTED,
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
                   RUNTIME_ERROR,
                   NOT_NORMALIZABLE,
                   INVALID_EVIDENCE_FILE,
                   INVALID_EMALG_FILE,
                   UNKNOWN_PARAMETER_ESTIMATION_METHOD,
                   NUM_ERRORS};  // NUM_ERRORS should be the last entry

        /// Constructor
        Exception( Code _code, const std::string& msg="", const std::string& detailedMsg="" ) : std::runtime_error(ErrorStrings[_code] + " [" +  msg + "]"), errorcode(_code) { 
            if( !detailedMsg.empty() ) 
                std::cerr << "EXCEPTION: " << detailedMsg << std::endl; 
        }
        
        /// Copy constructor
        Exception( const Exception &e ) : std::runtime_error(e), errorcode(e.errorcode) {}

        /// Returns error code of this exception
        Code code() const { return errorcode; }


    private:
        /// Contains the error code of this exception
        Code errorcode;

        /// Error messages corresponding to the exceptions enumerated above
        static std::string ErrorStrings[NUM_ERRORS];
};


}


#endif
