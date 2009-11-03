/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines general utility functions and adds an abstraction layer for platform-dependent functionality


#ifndef __defined_libdai_util_h
#define __defined_libdai_util_h


#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>


#if defined(WINDOWS)
    #include <map> // an alternative would be to use boost/tr1/unordered_map.hpp
#elif defined(CYGWIN)
    #include <boost/tr1/unordered_map.hpp> // only present in boost 1.37 and higher
#else
    #include <tr1/unordered_map> // only present in modern GCC distributions
#endif


/// An alias to the BOOST_FOREACH macro from the boost::foreach library
#define foreach BOOST_FOREACH

#ifdef DAI_DEBUG
/// \brief "Print variable". Prints the text of an expression, followed by its value (only if DAI_DEBUG is defined)
/**
 *  Useful debugging macro to see what your code is doing.
 *  Example: \code DAI_PV(3+4) \endcode
 *  Output: \code 3+4= 7 \endcode
 */
#define DAI_PV(x) do {std::cerr << #x "= " << (x) << std::endl;} while(0)
/// "Debugging message": Prints a message (only if DAI_DEBUG is defined)
#define DAI_DMSG(str) do {std::cerr << str << std::endl;} while(0)
#else
#define DAI_PV(x) do {} while(0)
#define DAI_DMSG(str) do {} while(0)
#endif

/// Macro to give error message \a stmt if props.verbose >= \a n
#define DAI_IFVERB(n, stmt) if(props.verbose>=n) { std::cerr << stmt; }


#ifdef WINDOWS
    /// Returns true if argument is NAN (Not A Number)
    bool isnan( double x );

    /// Returns inverse hyperbolic tangent of argument
    double atanh( double x );

    /// Returns log(1+x)
    double log1p( double x );

    /// Define INFINITY
    #define INFINITY (std::numeric_limits<Real>::infinity())
#endif


namespace dai {


/// Real number (alias for \c double, which could be changed to <tt>long double</tt> if necessary)
typedef double Real;

/// Returns logarithm of \a x
inline Real log( Real x ) {
    return std::log(x);
}

/// Returns logarithm of \a x, or 0 if \a x == 0
inline Real log0( Real x ) {
    return x ? std::log(x) : 0;
}

/// Returns exponent of \a x
inline Real exp( Real x ) {
    return std::exp(x);
}


#ifdef WINDOWS
    /// hash_map is an alias for \c std::map.
    /** Since there is no TR1 unordered_map implementation available yet, we fall back on std::map.
     */
    template <typename T, typename U, typename H = boost::hash<T> >
        class hash_map : public std::map<T,U> {};
#else
    /// hash_map is an alias for \c std::tr1::unordered_map.
    /** We use the (experimental) TR1 unordered_map implementation included in modern GCC distributions or in boost versions 1.37 and higher.
     */
    template <typename T, typename U, typename H = boost::hash<T> >
        class hash_map : public std::tr1::unordered_map<T,U,H> {};
#endif


/// Returns wall clock time in seconds
double toc();


/// Returns absolute value of \a t
template<class T>
inline T abs( const T &t ) {
    return (t < 0) ? (-t) : t;
}


/// Sets the random seed
void rnd_seed( size_t seed );

/// Returns a real number, distributed uniformly on [0,1)
Real rnd_uniform();

/// Returns a real number from a standard-normal distribution
Real rnd_stdnormal();

/// Returns a random integer in interval [\a min, \a max]
int rnd_int( int min, int max );

/// Returns a random integer in the half-open interval [0, \a n)
inline int rnd( int n) {
    return rnd_int( 0, n-1 );
}


/// Writes a \c std::vector<> to a \c std::ostream
template<class T>
std::ostream& operator << (std::ostream& os, const std::vector<T> & x) {
    os << "(";
    for( typename std::vector<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << ")";
    return os;
}

/// Writes a \c std::set<> to a \c std::ostream
template<class T>
std::ostream& operator << (std::ostream& os, const std::set<T> & x) {
    os << "{";
    for( typename std::set<T>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << *it;
    os << "}";
    return os;
}

/// Writes a \c std::map<> to a \c std::ostream
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::map<T1,T2> & x) {
    os << "{";
    for( typename std::map<T1,T2>::const_iterator it = x.begin(); it != x.end(); it++ )
        os << (it != x.begin() ? ", " : "") << it->first << "->" << it->second;
    os << "}";
    return os;
}

/// Writes a \c std::pair<> to a \c std::ostream
template<class T1, class T2>
std::ostream& operator << (std::ostream& os, const std::pair<T1,T2> & x) {
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}

/// Concatenates two vectors
template<class T>
std::vector<T> concat( const std::vector<T>& u, const std::vector<T>& v ) {
    std::vector<T> w;
    w.reserve( u.size() + v.size() );
    for( size_t i = 0; i < u.size(); i++ )
        w.push_back( u[i] );
    for( size_t i = 0; i < v.size(); i++ )
        w.push_back( v[i] );
    return w;
}

/// Split a string into tokens delimited by characters in \a delim
void tokenizeString( const std::string& s, std::vector<std::string>& outTokens, const std::string& delim="\t\n" );

/// Used to keep track of the progress made by iterative algorithms.
/** A Diffs object stores an array of fixed size, containing the 
 *  history of certain values (for example, the \f$\ell_\infty\f$ 
 *  differences between all the old beliefs and the new beliefs 
 *  in one pass of belief propagation updates). A new value can be
 *  registered by calling Diffs::push(); the BP algorithm would use
 *  this to register the difference between the old and a newly 
 *  calculated belief. The Diffs object keeps track of the maximum 
 *  value encountered in the last Diffs::maxSize() values registered.
 *  The maximum value can then be queried by Diffs::maxDiff(). The
 *  BP algorithm would use this maximum value to compare it with a
 *  given tolerance, and if the tolerance exceeds the maximum value,
 *  the algorithm has converged.
 */
class Diffs : public std::vector<Real> {
    private:
        size_t _maxsize;
        Real _def;
        std::vector<Real>::iterator _pos;
        std::vector<Real>::iterator _maxpos;
    public:
        /// Constructor
        /** \param maxsize Maximum number of differences to store
         *  \param def Default value
         */
        Diffs(long maxsize, Real def) : std::vector<Real>(), _maxsize(maxsize), _def(def) {
            this->reserve(_maxsize);
            _pos = begin();
            _maxpos = begin();
        }
        /// Returns maximum difference encountered so far
        Real maxDiff() {
            if( size() < _maxsize )
                return _def;
            else
                return( *_maxpos );
        }
        /// Register new difference \a x
        void push(Real x) {
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
            } else {
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
        /// Return maximum number of differences that can be stored
        size_t maxSize() { return _maxsize; }
};


} // end of namespace dai


#endif
