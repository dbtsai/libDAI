/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
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


#ifndef __ENUM_H
#define __ENUM_H


#include <cstring>
#include <iostream>


namespace dai {


// C++ enums are too limited for my purposes. This defines wrapper classes
// that provide much more functionality than a simple enum. The only
// disadvantage is that one wrapper class needs to be written for each
// number of values an enum can take... a better solution is needed.


#define ENUM2(x,a,b) class x {\
    public:\
        enum value {a, b};\
\
        x() : v(a) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
           static char const* labels[] = {#a, #b};\
           size_t i = 0;\
           for( ; i < sizeof(labels) / sizeof(char const *); i++ )\
               if( strcmp( w, labels[i] ) == 0 ) {\
                   v = (value)i;\
                   break;\
               }\
           if( i == sizeof(labels) / sizeof(char const *) )\
               throw "Unknown " #x " value";\
        }\
\
        operator value () const { return v; }\
\
        operator size_t () const { return (size_t)v; }\
\
        operator char const* () const {\
           static char const* labels[] = {#a, #b};\
           return labels[v];\
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


#define ENUM3(x,a,b,c) class x {\
    public:\
        enum value {a, b, c};\
\
        x() : v(a) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
           static char const* labels[] = {#a, #b, #c};\
           size_t i = 0;\
           for( ; i < sizeof(labels) / sizeof(char const *); i++ )\
               if( strcmp( w, labels[i] ) == 0 ) {\
                   v = (value)i;\
                   break;\
               }\
           if( i == sizeof(labels) / sizeof(char const *) )\
               throw "Unknown " #x " value";\
        }\
\
        operator value () const { return v; }\
\
        operator size_t () const { return (size_t)v; }\
\
        operator char const* () const {\
           static char const* labels[] = {#a, #b, #c};\
           return labels[v];\
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


#define ENUM4(x,a,b,c,d) class x {\
    public:\
        enum value {a, b, c, d};\
\
        x() : v(a) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
           static char const* labels[] = {#a, #b, #c, #d};\
           size_t i = 0;\
           for( ; i < sizeof(labels) / sizeof(char const *); i++ )\
               if( strcmp( w, labels[i] ) == 0 ) {\
                   v = (value)i;\
                   break;\
               }\
           if( i == sizeof(labels) / sizeof(char const *) )\
               throw "Unknown " #x " value";\
        }\
\
        operator value () const { return v; }\
\
        operator size_t () const { return (size_t)v; }\
\
        operator char const* () const {\
           static char const* labels[] = {#a, #b, #c, #d};\
           return labels[v];\
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


#define ENUM5(x,a,b,c,d,e) class x {\
    public:\
        enum value {a, b, c, d, e};\
\
        x() : v(a) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
           static char const* labels[] = {#a, #b, #c, #d, #e};\
           size_t i = 0;\
           for( ; i < sizeof(labels) / sizeof(char const *); i++ )\
               if( strcmp( w, labels[i] ) == 0 ) {\
                   v = (value)i;\
                   break;\
               }\
           if( i == sizeof(labels) / sizeof(char const *) )\
               throw "Unknown " #x " value";\
        }\
\
        operator value () const { return v; }\
\
        operator size_t () const { return (size_t)v; }\
\
        operator char const* () const {\
           static char const* labels[] = {#a, #b, #c, #d, #e};\
           return labels[v];\
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


#define ENUM6(x,a,b,c,d,e,f) class x {\
    public:\
        enum value {a, b, c, d, e, f};\
\
        x() : v(a) {}\
\
        x(value w) : v(w) {}\
\
        x(char const *w) {\
           static char const* labels[] = {#a, #b, #c, #d, #e, #f};\
           size_t i = 0;\
           for( ; i < sizeof(labels) / sizeof(char const *); i++ )\
               if( strcmp( w, labels[i] ) == 0 ) {\
                   v = (value)i;\
                   break;\
               }\
           if( i == sizeof(labels) / sizeof(char const *) )\
               throw "Unknown " #x " value";\
        }\
\
        operator value () const { return v; }\
\
        operator size_t () const { return (size_t)v; }\
\
        operator char const* () const {\
           static char const* labels[] = {#a, #b, #c, #d, #e, #f};\
           return labels[v];\
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


}


#endif
