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


#ifndef __defined_libdai_util_h
#define __defined_libdai_util_h


#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cstdio>
#include <boost/foreach.hpp>


#ifdef WINDOWS
    #include <map>
#else
    #include <tr1/unordered_map>
#endif


#define foreach BOOST_FOREACH


#ifdef WINDOWS
    bool isnan( double x );
    double atanh( double x );
    double log1p( double x );
#endif


namespace dai {


#ifdef WINDOWS
    template <typename T, typename U>
        class hash_map : public std::map<T,U> {};
#else
    template <typename T, typename U>
        class hash_map : public std::tr1::unordered_map<T,U> {};
#endif


double toc();
void rnd_seed( size_t seed );
double rnd_uniform();
double rnd_stdnormal();
int rnd_int( int min, int max );


// Output a std::vector
template<class T> 
std::ostream& operator << (std::ostream& os, const std::vector<T> & x) {
    os << "(";
    for( typename std::vector<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << ")";
    return os;
}


// Output a std::set
template<class T> 
std::ostream& operator << (std::ostream& os, const std::set<T> & x) {
    os << "{";
   for( typename std::set<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << "}";
    return os;
}


/// Output a std::map
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::map<T1,T2> & x) {
    os << "{";
    for( typename std::map<T1,T2>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << it->first << "->" << it->second;
    os << "}";
    return os;
}


/// Output a std::pair
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1,T2> & x) {
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}


} // end of namespace dai


#endif
