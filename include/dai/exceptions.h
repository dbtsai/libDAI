/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines the Exception class and macros for throwing exceptions and doing assertions


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

/// Macro that simplifies throwing an exception with a useful default error message. 
/** The error message consists of a description of the exception, the source 
 *  code file and line number where the exception has been thrown.
 *  \param cod Corresponds to one of the enum values of dai::Exception::Code
 *
 *  \par Example:
 *  \code
 *  DAI_THROW(NOT_IMPLEMENTED);
 *  \endcode
 */
#define DAI_THROW(cod) throw dai::Exception(dai::Exception::cod, std::string(__FILE__ ", line " DAI_TOSTRING(__LINE__)))

/// Macro that simplifies throwing an exception with a user-defined error message.
/** \param cod Corresponds to one of the enum values of dai::Exception::Code
 *  \param msg Detailed error message that will be written to std::cerr.
 *
 *  \par Example:
 *  \code
 *  DAI_THROWE(NOT_IMPLEMENTED,"Detailed error message");
 *  \endcode
 */
#define DAI_THROWE(cod,msg) throw dai::Exception(dai::Exception::cod, std::string(__FILE__ ", line " DAI_TOSTRING(__LINE__)), msg)

/// Assertion mechanism, similar to the standard assert() macro. It is always active, even if NDEBUG is defined
#define DAI_ASSERT(condition) ((condition) ? ((void)0) : DAI_THROWE(ASSERTION_FAILED, std::string("Assertion \"" #condition "\" failed")))

// Assertion only if DAI_DEBUG is defined
#ifdef DAI_DEBUG
/// Assertion mechanism similar to DAI_ASSERT which is only active if DAI_DEBUG is defined
#define DAI_DEBASSERT(x) do {DAI_ASSERT(x);} while(0)
#else
#define DAI_DEBASSERT(x) do {} while(0)
#endif


namespace dai {


/// Error handling in libDAI is done by throwing an instance of the Exception class.
/** The Exception class inherits from std::runtime_error. It defines several types of exceptions
 *  and corresponding error messages. The recommended way to throw an instance of the Exception
 *  class is by using the #DAI_THROW or #DAI_THROWE macros.
 */
class Exception : public std::runtime_error {
    public:
        /// Enumeration of exceptions used in libDAI
        enum Code {NOT_IMPLEMENTED,
                   ASSERTION_FAILED,
                   IMPOSSIBLE_TYPECAST,
                   OBJECT_NOT_FOUND,
                   BELIEF_NOT_AVAILABLE,
                   UNKNOWN_ENUM_VALUE,
                   UNKNOWN_DAI_ALGORITHM,
                   UNKNOWN_PARAMETER_ESTIMATION_METHOD,
                   UNKNOWN_PROPERTY_TYPE,
                   UNKNOWN_PROPERTY,
                   MALFORMED_PROPERTY,
                   NOT_ALL_PROPERTIES_SPECIFIED,
                   INVALID_ALIAS,
                   CANNOT_READ_FILE,
                   CANNOT_WRITE_FILE,
                   INVALID_FACTORGRAPH_FILE,
                   INVALID_EVIDENCE_FILE,
                   INVALID_EMALG_FILE,
                   NOT_NORMALIZABLE,
                   MULTIPLE_UNDO,
                   FACTORGRAPH_NOT_CONNECTED,
                   INTERNAL_ERROR,
                   RUNTIME_ERROR,
                   NUM_ERRORS};  // NUM_ERRORS should be the last entry

        /// Constructor
        Exception( Code _code, const std::string& msg="", const std::string& detailedMsg="" ) : std::runtime_error(ErrorStrings[_code] + " [" +  msg + "]"), errorcode(_code) {
            if( !detailedMsg.empty() )
                std::cerr << "EXCEPTION: " << detailedMsg << std::endl;
        }

        /// Returns error code of this exception
        Code code() const { return errorcode; }

        /// Returns error message corresponding to an error code
        const std::string &message( const Code c ) const { return ErrorStrings[c]; }


    private:
        /// Contains the error code of this exception
        Code errorcode;

        /// Error messages corresponding to the exceptions enumerated above
        static std::string ErrorStrings[NUM_ERRORS];
};


}


#endif
