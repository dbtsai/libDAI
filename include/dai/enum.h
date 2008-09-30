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
/// \brief Defines the DAI_ENUM macro


#ifndef __defined_libdai_enum_h
#define __defined_libdai_enum_h


#include <cstring>
#include <iostream>
#include <dai/exceptions.h>


/// Extends the C++ enum type by supporting input/output streaming and conversion to and from const char*
/** Example of usage:
 *  \code
 *      DAI_ENUM(colors,RED,GREEN,BLUE)
 *  \endcode
 *  defines a class encapsulating an
 *  \code
 *      enum colors {RED, GREEN, BLUE};
 *  \endcode
 *  It offers additional functionality over the plain "enum" keyword.
 */
#define DAI_ENUM(x,val0,...) class x {\
    public:\
        enum value {val0,__VA_ARGS__};\
\
        x() : v(val0) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
            static char const* labelstring = #val0 "," #__VA_ARGS__ ",";\
            size_t pos_begin = 0;\
            size_t i = 0;\
            for( size_t pos_end = 0; labelstring[pos_end] != '\0'; pos_end++ )\
                if( (labelstring[pos_end] == ',') ) {\
                    if( (strlen( w ) == pos_end - pos_begin) && (strncmp( labelstring + pos_begin, w, pos_end - pos_begin ) == 0) ) {\
                        v = (value)i;\
                        return;\
                    } else {\
                        i++;\
                        pos_begin = pos_end + 1;\
                    }\
                }\
            DAI_THROW(UNKNOWN_ENUM_VALUE);\
        }\
\
        operator value() const { return v; }\
\
        operator size_t() const { return (size_t)v; }\
\
        operator char const*() const {\
            static char labelstring[] = #val0 "," #__VA_ARGS__;\
            size_t pos_begin = 0;\
            size_t i = 0;\
            for( size_t pos_end = 0; ; pos_end++ )\
                if( (labelstring[pos_end] == ',') || (labelstring[pos_end] == '\0') ) {\
                    if( (size_t)v == i ) {\
                        labelstring[pos_end] = '\0';\
                        return labelstring + pos_begin;\
                    } else {\
                        i++;\
                        pos_begin = pos_end + 1;\
                    }\
                }\
        }\
\
        friend std::istream& operator >> (std::istream& is, x& y) {\
            std::string s;\
            is >> s;\
            y = x(s.c_str());\
            return is;\
        }\
\
        friend std::ostream& operator << (std::ostream& os, const x& y) {\
            os << (const char *)y;\
            return os;\
        }\
\
    private:\
        value v;\
};


#endif
