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
/// \brief Defines general utility functions and adds an abstraction layer for platform-dependent functionality


#ifndef __defined_libdai_util_h
#define __defined_libdai_util_h


#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cstdio>
#include <boost/foreach.hpp>
#include <algorithm>


#ifdef WINDOWS
    #include <map>
#else
    #include <tr1/unordered_map>
#endif


/// An alias to the BOOST_FOREACH macro from the boost::foreach library
#define foreach BOOST_FOREACH


#ifdef WINDOWS
    /// Returns true if argument is NAN (Not A Number)
    bool isnan( double x );

    /// Returns inverse hyperbolic tangent of argument
    double atanh( double x );

    /// Returns log(1+x)
    double log1p( double x );
#endif


namespace dai {


#ifdef WINDOWS
    /// hash_map is an alias for std::map.
    /** Since there is no TR1 unordered_map implementation available yet, we fall back on std::map.
     */
    template <typename T, typename U>
        class hash_map : public std::map<T,U> {};
#else
    /// hash_map is an alias for std::tr1::unordered_map.
    /** We use the (experimental) TR1 unordered_map implementation included in modern GCC distributions.
     */
    template <typename T, typename U>
        class hash_map : public std::tr1::unordered_map<T,U> {};
#endif


/// Returns the time in seconds
double toc();


/// Sets the random seed
void rnd_seed( size_t seed );

/// Returns a real number, distributed uniformly on [0,1)
double rnd_uniform();

/// Returns a real number from a standard-normal distribution
double rnd_stdnormal();

/// Returns a random integer in interval [min, max]
int rnd_int( int min, int max );


/// Writes a std::vector to a std::ostream
template<class T> 
std::ostream& operator << (std::ostream& os, const std::vector<T> & x) {
    os << "(";
    for( typename std::vector<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << ")";
    return os;
}

/// Writes a std::set to a std::ostream
template<class T> 
std::ostream& operator << (std::ostream& os, const std::set<T> & x) {
    os << "{";
   for( typename std::set<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << "}";
    return os;
}

/// Writes a std::map to a std::ostream
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::map<T1,T2> & x) {
    os << "{";
    for( typename std::map<T1,T2>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << it->first << "->" << it->second;
    os << "}";
    return os;
}

/// Writes a std::pair to a std::ostream
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1,T2> & x) {
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}


/// Used to keep track of the progress made by iterative algorithms
class Diffs : public std::vector<double> {
    private:
        size_t _maxsize;
        double _def;
        std::vector<double>::iterator _pos;
        std::vector<double>::iterator _maxpos;
    public:
        /// Constructor
        Diffs(long maxsize, double def) : std::vector<double>(), _maxsize(maxsize), _def(def) { 
            this->reserve(_maxsize); 
            _pos = begin(); 
            _maxpos = begin(); 
        }
        /// Returns maximum difference encountered
        double maxDiff() { 
            if( size() < _maxsize )
                return _def;
            else
                return( *_maxpos ); 
        }
        /// Register new difference x
        void push(double x) {
            if( size() < _maxsize ) {
                push_back(x);
                _pos = end();
                if( size() > 1 ) {
                    if( *_maxpos < back() ) {
                        _maxpos = end();
                        _maxpos--;
                    }
                } else {
                    _maxpos = begin();
                }
            }
            else {
                if( _pos == end() )
                    _pos = begin();
                if( _maxpos == _pos ) {
                    *_pos++ = x; 
                    _maxpos = max_element(begin(),end());
                } else {
                    if( x > *_maxpos )
                        _maxpos = _pos;
                    *_pos++ = x;
                }
            }
        }
        /// Return maximum number of differences stored
        size_t maxSize() { return _maxsize; }
};


} // end of namespace dai


#endif
