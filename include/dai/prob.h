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
/// \brief Defines TProb<T> and Prob classes


#ifndef __defined_libdai_prob_h
#define __defined_libdai_prob_h


#include <cmath>
#include <vector>
#include <ostream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <functional>
#include <dai/util.h>


namespace dai {


/// Real number (alias for double, could be changed to long double if necessary)
typedef double                  Real;

template<typename T> class      TProb;

/// Represents a probability measure, with entries of type Real.
typedef TProb<Real>             Prob;


/// Represents a probability measure on a finite outcome space (i.e., corresponding to a discrete random variable).
/** It is implemented as a std::vector<T> but adds a convenient interface.
 *  It is not necessarily normalized at all times.
 *  \tparam T Should be castable from and to double.
 */
template <typename T> class TProb {
    private:
        /// The probability measure
        std::vector<T> _p;

    public:
        /// Enumerates different ways of normalizing a probability measure.
        /** 
         *  - NORMPROB means that the sum of all entries should be 1;
         *  - NORMLINF means that the maximum absolute value of all entries should be 1.
         */
        typedef enum { NORMPROB, NORMLINF } NormType;
        /// Enumerates different distance measures between probability measures.
        /** 
         *  - DISTL1 is the L-1 distance (sum of absolute values of pointwise difference);
         *  - DISTLINF is the L-inf distance (maximum absolute value of pointwise difference);
         *  - DISTTV is the Total Variation distance;
         *  - DISTKL is the Kullback-Leibler distance.
         */
        typedef enum { DISTL1, DISTLINF, DISTTV, DISTKL } DistType;
        
        /// Default constructor
        TProb() : _p() {}
        
        /// Construct uniform distribution of given length
        explicit TProb( size_t n ) : _p(std::vector<T>(n, 1.0 / n)) {}
        
        /// Construct from given length and initial value
        TProb( size_t n, Real p ) : _p(n, (T)p) {}
        
        /// Construct from given length and initial array
        TProb( size_t n, const Real* p ) : _p(p, p + n ) {}
        
        /// Returns a const reference to the probability vector
        const std::vector<T> & p() const { return _p; }

        /// Returns a reference to the probability vector
        std::vector<T> & p() { return _p; }
        
        /// Returns a copy of the i'th probability entry
        T operator[]( size_t i ) const { 
#ifdef DAI_DEBUG
            return _p.at(i);
#else
            return _p[i];
#endif
        }
        
        /// Returns a reference to the i'th probability entry
        T& operator[]( size_t i ) { return _p[i]; }

        /// Sets all elements to x
        TProb<T> & fill(T x) { 
            std::fill( _p.begin(), _p.end(), x );
            return *this;
        }

        /// Sets all elements to i.i.d. random numbers from a uniform[0,1) distribution
        TProb<T> & randomize() { 
            std::generate(_p.begin(), _p.end(), rnd_uniform);
            return *this;
        }

        /// Returns number of elements
        size_t size() const {
            return _p.size();
        }

        /// Sets entries that are smaller than epsilon to zero
        TProb<T>& makeZero( Real epsilon ) {
            for( size_t i = 0; i < size(); i++ )
                if( fabs(_p[i]) < epsilon )
                    _p[i] = 0;
            return *this;
        }

        /// Sets entries that are smaller than epsilon to epsilon
        TProb<T>& makePositive (Real epsilon) {
            for( size_t i = 0; i < size(); i++ )
                if( (0 < (Real)_p[i]) && ((Real)_p[i] < epsilon) )
                    _p[i] = epsilon;
            return *this;
        }

        /// Multiplies each entry with x
        TProb<T>& operator*= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::multiplies<T>(), x) );
            return *this;
        }

        /// Returns product of *this with x
        TProb<T> operator* (T x) const {
            TProb<T> prod( *this );
            prod *= x;
            return prod;
        }

        /// Divides each entry by x
        TProb<T>& operator/= (T x) {
#ifdef DAI_DEBUG
            assert( x != 0.0 );
#endif
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::divides<T>(), x ) );
            return *this;
        }

        /// Returns quotient of *this and x
        TProb<T> operator/ (T x) const {
            TProb<T> quot( *this );
            quot /= x;
            return quot;
        }

        /// Adds x to each entry
        TProb<T>& operator+= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::plus<T>(), x ) );
            return *this;
        }

        /// Returns sum of *this and x
        TProb<T> operator+ (T x) const {
            TProb<T> sum( *this );
            sum += x;
            return sum;
        }

        /// Subtracts x from each entry
        TProb<T>& operator-= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::minus<T>(), x ) );
            return *this;
        }

        /// Returns difference of *this and x
        TProb<T> operator- (T x) const {
            TProb<T> diff( *this );
            diff -= x;
            return diff;
        }

        /// Pointwise comparison
        bool operator<= (const TProb<T> & q) const {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            for( size_t i = 0; i < size(); i++ )
                if( !(_p[i] <= q[i]) )
                    return false;
            return true;
        }

        /// Pointwise multiplication with q
        TProb<T>& operator*= (const TProb<T> & q) {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::multiplies<T>() );
            return *this;
        }
        
        /// Return product of *this with q
        TProb<T> operator* (const TProb<T> & q) const {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            TProb<T> prod( *this );
            prod *= q;
            return prod;
        }

        /// Pointwise addition with q
        TProb<T>& operator+= (const TProb<T> & q) {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::plus<T>() );
            return *this;
        }
        
        /// Return sum of *this and q
        TProb<T> operator+ (const TProb<T> & q) const {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            TProb<T> sum( *this );
            sum += q;
            return sum;
        }
        
        /// Pointwise subtraction of q
        TProb<T>& operator-= (const TProb<T> & q) {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::minus<T>() );
            return *this;
        }
        
        /// Return *this minus q
        TProb<T> operator- (const TProb<T> & q) const {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            TProb<T> diff( *this );
            diff -= q;
            return diff;
        }

        /// Pointwise division by q, where division by zero yields zero
        TProb<T>& operator/= (const TProb<T> & q) {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            for( size_t i = 0; i < size(); i++ ) {
                if( q[i] == 0.0 )
                    _p[i] = 0.0;
                else
                    _p[i] /= q[i];
            }
            return *this;
        }
        
        /// Pointwise division by q, where division by zero yields infinity
        TProb<T>& divide (const TProb<T> & q) {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::divides<T>() );
            return *this;
        }
        
        /// Returns quotient of *this with q
        TProb<T> operator/ (const TProb<T> & q) const {
#ifdef DAI_DEBUG
            assert( size() == q.size() );
#endif
            TProb<T> quot( *this );
            quot /= q;
            return quot;
        }

        /// Returns pointwise inverse
        TProb<T> inverse(bool zero = false) const {
            TProb<T> inv;
            inv._p.reserve( size() );
            if( zero )
                for( size_t i = 0; i < size(); i++ )
                    inv._p.push_back( _p[i] == 0.0 ? 0.0 : 1.0 / _p[i] );
            else
                for( size_t i = 0; i < size(); i++ ) {
#ifdef DAI_DEBUG
                    assert( _p[i] != 0.0 );
#endif
                    inv._p.push_back( 1.0 / _p[i] );
                }
            return inv;
        }

        /// Raises elements to the power a
        TProb<T>& operator^= (Real a) {
            if( a != 1.0 )
                std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::ptr_fun<T, Real, T>(std::pow), a) );
            return *this;
        }

        /// Returns *this raised to the power a
        TProb<T> operator^ (Real a) const {
            TProb<T> power(*this);
            power ^= a;
            return power;
        }

        /// Returns pointwise signum
        TProb<T> sgn() const {
            TProb<T> x;
            x._p.reserve( size() );
            for( size_t i = 0; i < size(); i++ ) {
                T s = 0;
                if( _p[i] > 0 )
                    s = 1;
                else if( _p[i] < 0 )
                    s = -1;
                x._p.push_back( s );
            }
            return x;
        }

        /// Returns pointwise absolute value
        TProb<T> abs() const {
            TProb<T> x;
            x._p.reserve( size() );
            for( size_t i = 0; i < size(); i++ )
                x._p.push_back( _p[i] < 0 ? (-p[i]) : p[i] );
            return x;
        }

        /// Applies exp pointwise
        const TProb<T>& takeExp() {
            std::transform( _p.begin(), _p.end(), _p.begin(),  std::ptr_fun<T, T>(std::exp) );
            return *this;
        }

        /// Applies log pointwise
        const TProb<T>& takeLog() {
            std::transform( _p.begin(), _p.end(), _p.begin(),  std::ptr_fun<T, T>(std::log) );
            return *this;
        }

        /// Applies log pointwise (defining log(0)=0)
        const TProb<T>& takeLog0()  {
            for( size_t i = 0; i < size(); i++ )
               _p[i] = ( (_p[i] == 0.0) ? 0.0 : std::log( _p[i] ) );
            return *this;
        }

        /// Returns pointwise exp
        TProb<T> exp() const {
            TProb<T> e(*this);
            e.takeExp();
            return e;
        }

        /// Returns pointwise log
        TProb<T> log() const {
            TProb<T> l(*this);
            l.takeLog();
            return l;
        }

        /// Returns pointwise log (defining log(0)=0)
        TProb<T> log0() const {
            TProb<T> l0(*this);
            l0.takeLog0();
            return l0;
        }

        /// Returns distance of p and q, measured using dt
        friend Real dist( const TProb<T> &p, const TProb<T> &q, DistType dt ) {
#ifdef DAI_DEBUG
            assert( p.size() == q.size() );
#endif
            Real result = 0.0;
            switch( dt ) {
                case DISTL1:
                    for( size_t i = 0; i < p.size(); i++ )
                        result += fabs((Real)p[i] - (Real)q[i]);
                    break;
                    
                case DISTLINF:
                    for( size_t i = 0; i < p.size(); i++ ) {
                        Real z = fabs((Real)p[i] - (Real)q[i]);
                        if( z > result )
                            result = z;
                    }
                    break;

                case DISTTV:
                    for( size_t i = 0; i < p.size(); i++ )
                        result += fabs((Real)p[i] - (Real)q[i]);
                    result *= 0.5;
                    break;

                case DISTKL:
                    for( size_t i = 0; i < p.size(); i++ ) {
                        if( p[i] != 0.0 )
                            result += p[i] * (std::log(p[i]) - std::log(q[i]));
                    }
            }
            return result;
        }

        /// Returns sum of all entries
        T totalSum() const {
            T Z = std::accumulate( _p.begin(),  _p.end(), (T)0 );
            return Z;
        }

        /// Returns maximum absolute value of entries
        T maxAbs() const {
            T Z = 0;
            for( size_t i = 0; i < size(); i++ ) {
                Real mag = fabs( (Real) _p[i] );
                if( mag > Z )
                    Z = mag;
            }
            return Z;
        }

        /// Returns maximum value of entries
        T maxVal() const {
            T Z = *std::max_element( _p.begin(), _p.end() );
            return Z;
        }

        /// Returns minimum value of entries
        T minVal() const {
            T Z = *std::min_element( _p.begin(), _p.end() );
            return Z;
        }

        /// Normalizes using the specified norm
        T normalize( NormType norm = NORMPROB ) {
            T Z = 0.0;
            if( norm == NORMPROB )
                Z = totalSum();
            else if( norm == NORMLINF )
                Z = maxAbs();
#ifdef DAI_DEBUG
            assert( Z != 0.0 );
#endif
            *this /= Z;
            return Z;
        }

        /// Returns normalized copy of *this, using the specified norm
        TProb<T> normalized( NormType norm = NORMPROB ) const {
            TProb<T> result(*this);
            result.normalize( norm );
            return result;
        }
    
        /// Returns true if one or more entries are NaN
        bool hasNaNs() const {
            return (std::find_if( _p.begin(), _p.end(), isnan ) != _p.end());
        }

        /// Returns true if one or more entries are negative
        bool hasNegatives() const {
            return (std::find_if( _p.begin(), _p.end(), std::bind2nd( std::less<Real>(), 0.0 ) ) != _p.end());
        }
        
        /// Returns true if one or more entries are non-positive (causes problems with logscale)
        bool hasNonPositives() const {
            return (std::find_if( _p.begin(), _p.end(), std::bind2nd( std::less_equal<Real>(), 0.0 ) ) != _p.end());
        }

        /// Returns entropy
        Real entropy() const {
            Real S = 0.0;
            for( size_t i = 0; i < size(); i++ )
                S -= xlogx(_p[i]);
            return S;
        }

        /// Writes a TProb<T> to an output stream
        friend std::ostream& operator<< (std::ostream& os, const TProb<T>& P) {
            os << "[";
            std::copy( P._p.begin(), P._p.end(), std::ostream_iterator<T>(os, " ") );
            os << "]";
            return os;
        }

    private:
        /// Returns x*log(x), or 0 if x == 0
        Real xlogx( Real x ) const { return( x == 0.0 ? 0.0 : x * std::log(x)); }
};


/// Returns TProb<T> containing the pointwise minimum of a and b (which should have equal size)
template<typename T> TProb<T> min( const TProb<T> &a, const TProb<T> &b ) {
    assert( a.size() == b.size() );
    TProb<T> result( a.size() );
    for( size_t i = 0; i < a.size(); i++ )
        if( a[i] < b[i] )
            result[i] = a[i];
        else
            result[i] = b[i];
    return result;
}


/// Returns TProb<T> containing the pointwise maximum of a and b (which should have equal size)
template<typename T> TProb<T> max( const TProb<T> &a, const TProb<T> &b ) {
    assert( a.size() == b.size() );
    TProb<T> result( a.size() );
    for( size_t i = 0; i < a.size(); i++ )
        if( a[i] > b[i] )
            result[i] = a[i];
        else
            result[i] = b[i];
    return result;
}


} // end of namespace dai


#endif
