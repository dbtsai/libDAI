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
/// \todo Rename to Vector<T>


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
#include <dai/exceptions.h>


namespace dai {


/// Represents a vector with entries of type \a T.
/** A TProb<T> is a std::vector<T> with an interface designed for dealing with probability mass functions.
 *  It is mainly used for representing measures on a finite outcome space, e.g., the probability 
 *  distribution of a discrete random variable.
 *  \tparam T Should be a scalar that is castable from and to double and should support elementary arithmetic operations.
 */
template <typename T> class TProb {
    private:
        /// The vector
        std::vector<T> _p;

    public:
        /// Iterator over entries
    	typedef typename std::vector<T>::iterator iterator;
        /// Const iterator over entries
    	typedef typename std::vector<T>::const_iterator const_iterator;

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
        
        /// Construct uniform distribution over n outcomes, i.e., a vector of length n with each entry set to 1/n
        explicit TProb( size_t n ) : _p(std::vector<T>(n, 1.0 / n)) {}
        
        /// Construct vector of length n with each entry set to p
        explicit TProb( size_t n, Real p ) : _p(n, (T)p) {}
        
        /// Construct vector from a range
        /** \tparam Iterator Iterates over instances that can be cast to T.
         *  \param begin Points to first instance to be added.
         *  \param end Points just beyond last instance to be added.
         *  \param sizeHint For efficiency, the number of entries can be speficied by sizeHint.
         */
        template <typename Iterator>
        TProb( Iterator begin, Iterator end, size_t sizeHint=0 ) : _p() {
            _p.reserve( sizeHint );
            _p.insert( _p.begin(), begin, end );
        }
        
        /// Returns a const reference to the vector
        const std::vector<T> & p() const { return _p; }

        /// Returns a reference to the vector
        std::vector<T> & p() { return _p; }
        
        /// Returns a copy of the i'th entry
        T operator[]( size_t i ) const { 
#ifdef DAI_DEBUG
            return _p.at(i);
#else
            return _p[i];
#endif
        }
        
        /// Returns reference to the i'th entry
        T& operator[]( size_t i ) { return _p[i]; }
        
        /// Returns iterator pointing to first entry
        iterator begin() { return _p.begin(); }

        /// Returns const iterator pointing to first entry
        const_iterator begin() const { return _p.begin(); }

        /// Returns iterator pointing beyond last entry
        iterator end() { return _p.end(); }

        /// Returns const iterator pointing beyond last entry
        const_iterator end() const { return _p.end(); }

        /// Sets all entries to x
        TProb<T> & fill(T x) { 
            std::fill( _p.begin(), _p.end(), x );
            return *this;
        }

        /// Draws all entries i.i.d. from a uniform distribution on [0,1)
        TProb<T> & randomize() { 
            std::generate(_p.begin(), _p.end(), rnd_uniform);
            return *this;
        }

        /// Returns length of the vector, i.e., the number of entries
        size_t size() const {
            return _p.size();
        }

        /// Sets entries that are smaller than epsilon to 0
        TProb<T>& makeZero( Real epsilon ) {
            for( size_t i = 0; i < size(); i++ )
                if( fabs(_p[i]) < epsilon )
                    _p[i] = 0;
            return *this;
        }
        
        /// Set all entries to 1.0/size()
        TProb<T>& setUniform () {
            fill(1.0/size());
            return *this;
        }

        /// Sets entries that are smaller than epsilon to epsilon
        TProb<T>& makePositive( Real epsilon ) {
            for( size_t i = 0; i < size(); i++ )
                if( (0 < (Real)_p[i]) && ((Real)_p[i] < epsilon) )
                    _p[i] = epsilon;
            return *this;
        }

        /// Multiplies each entry with scalar x
        TProb<T>& operator*= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::multiplies<T>(), x) );
            return *this;
        }

        /// Returns product of *this with scalar x
        TProb<T> operator* (T x) const {
            TProb<T> prod( *this );
            prod *= x;
            return prod;
        }

        /// Divides each entry by scalar x
        TProb<T>& operator/= (T x) {
            DAI_DEBASSERT( x != 0.0 );
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::divides<T>(), x ) );
            return *this;
        }

        /// Returns quotient of *this and scalar x
        TProb<T> operator/ (T x) const {
            TProb<T> quot( *this );
            quot /= x;
            return quot;
        }

        /// Adds scalar x to each entry
        TProb<T>& operator+= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::plus<T>(), x ) );
            return *this;
        }

        /// Returns sum of *this and scalar x
        TProb<T> operator+ (T x) const {
            TProb<T> sum( *this );
            sum += x;
            return sum;
        }

        /// Subtracts scalar x from each entry
        TProb<T>& operator-= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::minus<T>(), x ) );
            return *this;
        }

        /// Returns difference of *this and scalar x
        TProb<T> operator- (T x) const {
            TProb<T> diff( *this );
            diff -= x;
            return diff;
        }

        /// Lexicographical comparison (sizes should be identical)
        bool operator<= (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            for( size_t i = 0; i < size(); i++ )
                if( !(_p[i] <= q[i]) )
                    return false;
            return true;
        }

        /// Pointwise multiplication with q (sizes should be identical)
        TProb<T>& operator*= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::multiplies<T>() );
            return *this;
        }
        
        /// Return product of *this with q (sizes should be identical)
        TProb<T> operator* (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> prod( *this );
            prod *= q;
            return prod;
        }

        /// Pointwise addition with q (sizes should be identical)
        TProb<T>& operator+= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::plus<T>() );
            return *this;
        }
        
        /// Returns sum of *this and q (sizes should be identical)
        TProb<T> operator+ (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> sum( *this );
            sum += q;
            return sum;
        }
        
        /// Pointwise subtraction of q (sizes should be identical)
        TProb<T>& operator-= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::minus<T>() );
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> diff( *this );
            diff -= q;
            return diff;
        }

        /// Pointwise division by q, where division by 0 yields 0 (sizes should be identical)
        TProb<T>& operator/= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            for( size_t i = 0; i < size(); i++ ) {
                if( q[i] == 0.0 )
                    _p[i] = 0.0;
                else
                    _p[i] /= q[i];
            }
            return *this;
        }
        
        /// Pointwise division by q, where division by 0 yields +Inf (sizes should be identical)
        TProb<T>& divide (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::divides<T>() );
            return *this;
        }
        
        /// Returns quotient of *this with q (sizes should be identical)
        TProb<T> operator/ (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> quot( *this );
            quot /= q;
            return quot;
        }

        /// Returns pointwise inverse
        /** If zero==true; uses 1/0==0, otherwise 1/0==Inf.
         */
        TProb<T> inverse(bool zero=true) const {
            TProb<T> inv;
            inv._p.reserve( size() );
            if( zero )
                for( size_t i = 0; i < size(); i++ )
                    inv._p.push_back( _p[i] == 0.0 ? 0.0 : 1.0 / _p[i] );
            else
                for( size_t i = 0; i < size(); i++ )
                    inv._p.push_back( 1.0 / _p[i] );
            return inv;
        }

        /// Raises entries to the power a
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
                x._p.push_back( _p[i] < 0 ? (-_p[i]) : _p[i] );
            return x;
        }

        /// Applies exp pointwise
        const TProb<T>& takeExp() {
            std::transform( _p.begin(), _p.end(), _p.begin(),  std::ptr_fun<T, T>(std::exp) );
            return *this;
        }

        /// Applies log pointwise
        /** If zero==true, uses log(0)==0; otherwise, log(0)==-Inf.
         */
        const TProb<T>& takeLog(bool zero=false) {
            if( zero ) {
                for( size_t i = 0; i < size(); i++ )
                   _p[i] = ( (_p[i] == 0.0) ? 0.0 : std::log( _p[i] ) );
            } else
                std::transform( _p.begin(), _p.end(), _p.begin(),  std::ptr_fun<T, T>(std::log) );
            return *this;
        }

        /// Returns pointwise exp
        TProb<T> exp() const {
            TProb<T> e(*this);
            e.takeExp();
            return e;
        }

        /// Returns pointwise log
        /** If zero==true, uses log(0)==0; otherwise, log(0)==-Inf.
         */
        TProb<T> log(bool zero=false) const {
            TProb<T> l(*this);
            l.takeLog(zero);
            return l;
        }

        /// Returns sum of all entries
        T sum() const {
            T Z = std::accumulate( _p.begin(),  _p.end(), (T)0 );
            return Z;
        }

        /// Return sum of absolute value of all entries
        T sumAbs() const {
            T s = 0;
            for( size_t i = 0; i < size(); i++ )
                s += fabs( (Real) _p[i] );
            return s;
        }

        /// Returns maximum absolute value of all entries
        T maxAbs() const {
            T Z = 0;
            for( size_t i = 0; i < size(); i++ ) {
                Real mag = fabs( (Real) _p[i] );
                if( mag > Z )
                    Z = mag;
            }
            return Z;
        }

        /// Returns maximum value of all entries
        T max() const {
            T Z = *std::max_element( _p.begin(), _p.end() );
            return Z;
        }

        /// Returns minimum value of all entries
        T min() const {
            T Z = *std::min_element( _p.begin(), _p.end() );
            return Z;
        }

        /// Returns {arg,}maximum value
        std::pair<size_t,T> argmax() const {
            T max = _p[0];
            size_t arg = 0;
            for( size_t i = 1; i < size(); i++ ) {
              if( _p[i] > max ) {
                max = _p[i];
                arg = i;
              }
            }
            return std::make_pair(arg,max);
        }

        /// Normalizes vector using the specified norm
        T normalize( NormType norm=NORMPROB ) {
            T Z = 0.0;
            if( norm == NORMPROB )
                Z = sum();
            else if( norm == NORMLINF )
                Z = maxAbs();
            if( Z == 0.0 )
                DAI_THROW(NOT_NORMALIZABLE);
            else
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
            bool foundnan = false;
            for( typename std::vector<T>::const_iterator x = _p.begin(); x != _p.end(); x++ )
                if( isnan( *x ) ) {
                    foundnan = true;
                    break;
                }
            return foundnan;
        }

        /// Returns true if one or more entries are negative
        bool hasNegatives() const {
            return (std::find_if( _p.begin(), _p.end(), std::bind2nd( std::less<Real>(), 0.0 ) ) != _p.end());
        }
        
        /// Returns entropy of *this
        Real entropy() const {
            Real S = 0.0;
            for( size_t i = 0; i < size(); i++ )
                S -= (_p[i] == 0 ? 0 : _p[i] * std::log(_p[i]));
            return S;
        }

        /// Returns a random index, according to the (normalized) distribution described by *this
        size_t draw() {
            double x = rnd_uniform() * sum();
            T s = 0;
            for( size_t i = 0; i < size(); i++ ) {
                s += _p[i];
                if( s > x ) 
                    return i;
            }
            return( size() - 1 );
        }
};


/// Returns distance of p and q (sizes should be identical), measured using distance measure dt
/** \relates TProb
 */
template<typename T> Real dist( const TProb<T> &p, const TProb<T> &q, typename TProb<T>::DistType dt ) {
    DAI_DEBASSERT( p.size() == q.size() );
    Real result = 0.0;
    switch( dt ) {
        case TProb<T>::DISTL1:
            for( size_t i = 0; i < p.size(); i++ )
                result += fabs((Real)p[i] - (Real)q[i]);
            break;
            
        case TProb<T>::DISTLINF:
            for( size_t i = 0; i < p.size(); i++ ) {
                Real z = fabs((Real)p[i] - (Real)q[i]);
                if( z > result )
                    result = z;
            }
            break;

        case TProb<T>::DISTTV:
            for( size_t i = 0; i < p.size(); i++ )
                result += fabs((Real)p[i] - (Real)q[i]);
            result *= 0.5;
            break;

        case TProb<T>::DISTKL:
            for( size_t i = 0; i < p.size(); i++ ) {
                if( p[i] != 0.0 )
                    result += p[i] * (std::log(p[i]) - std::log(q[i]));
            }
    }
    return result;
}


/// Writes a TProb<T> to an output stream
/** \relates TProb
 */
template<typename T> std::ostream& operator<< (std::ostream& os, const TProb<T>& P) {
    os << "[";
    std::copy( P.p().begin(), P.p().end(), std::ostream_iterator<T>(os, " ") );
    os << "]";
    return os;
}


/// Returns the TProb<T> containing the pointwise minimum of a and b (which should have equal size)
/** \relates TProb
 */
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


/// Returns the TProb<T> containing the pointwise maximum of a and b (which should have equal size)
/** \relates TProb
 */
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


/// Represents a vector with entries of type Real.
typedef TProb<Real> Prob;


} // end of namespace dai


#endif
