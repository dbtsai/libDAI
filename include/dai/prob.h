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
/// \brief Defines TProb<> and Prob classes which represent (probability) vectors
/// \todo Rename to Vector<>


#ifndef __defined_libdai_prob_h
#define __defined_libdai_prob_h


#include <cmath>
#include <vector>
#include <ostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <dai/util.h>
#include <dai/exceptions.h>


namespace dai {


/// Represents a vector with entries of type \a T.
/** It is simply a <tt>std::vector</tt><<em>T</em>> with an interface designed for dealing with probability mass functions.
 *
 *  It is mainly used for representing measures on a finite outcome space, for example, the probability
 *  distribution of a discrete random variable. However, entries are not necessarily non-negative; it is also used to
 *  represent logarithms of probability mass functions.
 *
 *  \tparam T Should be a scalar that is castable from and to dai::Real and should support elementary arithmetic operations.
 */
template <typename T> class TProb {
    private:
        /// The vector
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
         *  - DISTL1 is the \f$\ell_1\f$ distance (sum of absolute values of pointwise difference);
         *  - DISTLINF is the \f$\ell_\infty\f$ distance (maximum absolute value of pointwise difference);
         *  - DISTTV is the total variation distance (half of the \f$\ell_1\f$ distance);
         *  - DISTKL is the Kullback-Leibler distance (\f$\sum_i p_i (\log p_i - \log q_i)\f$).
         */
        typedef enum { DISTL1, DISTLINF, DISTTV, DISTKL } DistType;

    /// \name Constructors and destructors
    //@{
        /// Default constructor (constructs empty vector)
        TProb() : _p() {}

        /// Construct uniform probability distribution over \a n outcomes (i.e., a vector of length \a n with each entry set to \f$1/n\f$)
        explicit TProb( size_t n ) : _p(std::vector<T>(n, (T)1 / n)) {}

        /// Construct vector of length \a n with each entry set to \a p
        explicit TProb( size_t n, T p ) : _p(n, p) {}

        /// Construct vector from a range
        /** \tparam TIterator Iterates over instances that can be cast to \a T
         *  \param begin Points to first instance to be added.
         *  \param end Points just beyond last instance to be added.
         *  \param sizeHint For efficiency, the number of entries can be speficied by \a sizeHint.
         */
        template <typename TIterator>
        TProb( TIterator begin, TIterator end, size_t sizeHint=0 ) : _p() {
            _p.reserve( sizeHint );
            _p.insert( _p.begin(), begin, end );
        }

        /// Construct vector from another vector
        /** \tparam S type of elements in \a v (should be castable to type \a T)
         *  \param v vector used for initialization
         */
        template <typename S>
        TProb( const std::vector<S> &v ) : _p() {
            _p.reserve( v.size() );
            _p.insert( _p.begin(), v.begin(), v.end() );
        }
    //@}

        /// Constant iterator over the elements
        typedef typename std::vector<T>::const_iterator const_iterator;
        /// Iterator over the elements
        typedef typename std::vector<T>::iterator iterator;
        /// Constant reverse iterator over the elements
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
        /// Reverse iterator over the elements
        typedef typename std::vector<T>::reverse_iterator reverse_iterator;

    /// @name Iterator interface
    //@{
        /// Returns iterator that points to the first element
        iterator begin() { return _p.begin(); }
        /// Returns constant iterator that points to the first element
        const_iterator begin() const { return _p.begin(); }

        /// Returns iterator that points beyond the last element
        iterator end() { return _p.end(); }
        /// Returns constant iterator that points beyond the last element
        const_iterator end() const { return _p.end(); }

        /// Returns reverse iterator that points to the last element
        reverse_iterator rbegin() { return _p.rbegin(); }
        /// Returns constant reverse iterator that points to the last element
        const_reverse_iterator rbegin() const { return _p.rbegin(); }

        /// Returns reverse iterator that points beyond the first element
        reverse_iterator rend() { return _p.rend(); }
        /// Returns constant reverse iterator that points beyond the first element
        const_reverse_iterator rend() const { return _p.rend(); }
    //@}

    /// \name Queries
    //@{
        /// Returns a const reference to the wrapped vector
        const std::vector<T> & p() const { return _p; }

        /// Returns a reference to the wrapped vector
        std::vector<T> & p() { return _p; }

        /// Returns a copy of the \a i 'th entry
        T operator[]( size_t i ) const {
#ifdef DAI_DEBUG
            return _p.at(i);
#else
            return _p[i];
#endif
        }

        /// Returns reference to the \a i 'th entry
        T& operator[]( size_t i ) { return _p[i]; }

        /// Returns length of the vector (i.e., the number of entries)
        size_t size() const { return _p.size(); }

        /// Returns the Shannon entropy of \c *this, \f$-\sum_i p_i \log p_i\f$
        T entropy() const {
            T S = 0;
            for( size_t i = 0; i < size(); i++ )
                S -= (_p[i] == 0 ? 0 : _p[i] * dai::log(_p[i]));
            return S;
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

        /// Returns a pair consisting of the index of the maximum value and the maximum value itself
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

        /// Returns sum of all entries
        T sum() const {
            T Z = std::accumulate( _p.begin(),  _p.end(), (T)0 );
            return Z;
        }

        /// Return sum of absolute value of all entries
        T sumAbs() const {
            T s = 0;
            for( size_t i = 0; i < size(); i++ )
                s += dai::abs(_p[i]);
            return s;
        }

        /// Returns maximum absolute value of all entries
        T maxAbs() const {
            T Z = 0;
            for( size_t i = 0; i < size(); i++ ) {
                T mag = dai::abs(_p[i]);
                if( mag > Z )
                    Z = mag;
            }
            return Z;
        }

        /// Returns \c true if one or more entries are NaN
        bool hasNaNs() const {
            bool foundnan = false;
            for( typename std::vector<T>::const_iterator x = _p.begin(); x != _p.end(); x++ )
                if( isnan( *x ) ) {
                    foundnan = true;
                    break;
                }
            return foundnan;
        }

        /// Returns \c true if one or more entries are negative
        bool hasNegatives() const {
            return (std::find_if( _p.begin(), _p.end(), std::bind2nd( std::less<T>(), (T)0 ) ) != _p.end());
        }

        /// Returns a random index, according to the (normalized) distribution described by *this
        size_t draw() {
            Real x = rnd_uniform() * sum();
            T s = 0;
            for( size_t i = 0; i < size(); i++ ) {
                s += _p[i];
                if( s > x )
                    return i;
            }
            return( size() - 1 );
        }

        /// Lexicographical comparison
        /** \pre <tt>this->size() == q.size()</tt>
         */
        bool operator<= (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            for( size_t i = 0; i < size(); i++ )
                if( !(_p[i] <= q[i]) )
                    return false;
            return true;
        }
    //@}

    /// \name Unary transformations
    //@{
        // OBSOLETE
        /// Returns pointwise signum
        /** \note Obsolete, to be removed soon
         */
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
                x._p.push_back( abs(_p[i]) );
            return x;
        }

        /// Returns pointwise exponent
        TProb<T> exp() const {
            TProb<T> e(*this);
            e.takeExp();
            return e;
        }

        /// Returns pointwise logarithm
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        TProb<T> log(bool zero=false) const {
            TProb<T> l(*this);
            l.takeLog(zero);
            return l;
        }

        /// Returns pointwise inverse
        /** If \a zero == \c true, uses <tt>1/0==0</tt>; otherwise, <tt>1/0==Inf</tt>.
         */
        TProb<T> inverse(bool zero=true) const {
            TProb<T> inv;
            inv._p.reserve( size() );
            if( zero )
                for( size_t i = 0; i < size(); i++ )
                    inv._p.push_back( _p[i] == (T)0 ? (T)0 : (T)1 / _p[i] );
            else
                for( size_t i = 0; i < size(); i++ )
                    inv._p.push_back( (T)1 / _p[i] );
            return inv;
        }

        /// Returns normalized copy of \c *this, using the specified norm
        TProb<T> normalized( NormType norm = NORMPROB ) const {
            TProb<T> result(*this);
            result.normalize( norm );
            return result;
        }
    //@}

    /// \name Unary operations
    //@{
        /// Draws all entries i.i.d. from a uniform distribution on [0,1)
        TProb<T>& randomize() {
            std::generate( _p.begin(), _p.end(), rnd_uniform );
            return *this;
        }

        /// Sets all entries to \f$1/n\f$ where \a n is the length of the vector
        TProb<T>& setUniform () {
            fill( (T)1 / size() );
            return *this;
        }

        /// Applies exponent pointwise
        const TProb<T>& takeExp() {
            for( size_t i = 0; i < size(); i++ )
               _p[i] = dai::exp(_p[i]);
            return *this;
        }

        /// Applies logarithm pointwise
        /** If \a zero == \c true, uses <tt>log(0)==0</tt>; otherwise, <tt>log(0)==-Inf</tt>.
         */
        const TProb<T>& takeLog(bool zero=false) {
            if( zero ) {
                for( size_t i = 0; i < size(); i++ )
                   _p[i] = ( (_p[i] == 0.0) ? 0.0 : dai::log(_p[i]) );
            } else
                for( size_t i = 0; i < size(); i++ )
                   _p[i] = dai::log(_p[i]);
            return *this;
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
    //@}

    /// \name Operations with scalars
    //@{
        /// Sets all entries to \a x
        TProb<T> & fill(T x) {
            std::fill( _p.begin(), _p.end(), x );
            return *this;
        }

        // OBSOLETE
        /// Sets entries that are smaller (in absolute value) than \a epsilon to 0
        /** \note Obsolete, to be removed soon
         */
        TProb<T>& makeZero( T epsilon ) {
            for( size_t i = 0; i < size(); i++ )
                if( (_p[i] < epsilon) && (_p[i] > -epsilon) )
                    _p[i] = 0;
            return *this;
        }
        
        // OBSOLETE
        /// Sets entries that are smaller than \a epsilon to \a epsilon
        /** \note Obsolete, to be removed soon
         */
        TProb<T>& makePositive( T epsilon ) {
            for( size_t i = 0; i < size(); i++ )
                if( (0 < _p[i]) && (_p[i] < epsilon) )
                    _p[i] = epsilon;
            return *this;
        }

        /// Adds scalar \a x to each entry
        TProb<T>& operator+= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::plus<T>(), x ) );
            return *this;
        }

        /// Subtracts scalar \a x from each entry
        TProb<T>& operator-= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::minus<T>(), x ) );
            return *this;
        }

        /// Multiplies each entry with scalar \a x
        TProb<T>& operator*= (T x) {
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::multiplies<T>(), x) );
            return *this;
        }

        /// Divides each entry by scalar \a x
        TProb<T>& operator/= (T x) {
            DAI_DEBASSERT( x != 0 );
            std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::divides<T>(), x ) );
            return *this;
        }

        /// Raises entries to the power \a x
        TProb<T>& operator^= (T x) {
            if( x != (T)1 )
                std::transform( _p.begin(), _p.end(), _p.begin(), std::bind2nd( std::ptr_fun<T, T, T>(std::pow), x) );
            return *this;
        }
    //@}

    /// \name Transformations with scalars
    //@{
        /// Returns sum of \c *this and scalar \a x
        TProb<T> operator+ (T x) const {
            TProb<T> sum( *this );
            sum += x;
            return sum;
        }

        /// Returns difference of \c *this and scalar \a x
        TProb<T> operator- (T x) const {
            TProb<T> diff( *this );
            diff -= x;
            return diff;
        }

        /// Returns product of \c *this with scalar \a x
        TProb<T> operator* (T x) const {
            TProb<T> prod( *this );
            prod *= x;
            return prod;
        }

        /// Returns quotient of \c *this and scalar \a x
        TProb<T> operator/ (T x) const {
            TProb<T> quot( *this );
            quot /= x;
            return quot;
        }

        /// Returns \c *this raised to the power \a x
        TProb<T> operator^ (T x) const {
            TProb<T> power(*this);
            power ^= x;
            return power;
        }
    //@}

    /// \name Operations with other equally-sized vectors
    //@{
        /// Pointwise addition with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T>& operator+= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::plus<T>() );
            return *this;
        }

        /// Pointwise subtraction of \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T>& operator-= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::minus<T>() );
            return *this;
        }

        /// Pointwise multiplication with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T>& operator*= (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::multiplies<T>() );
            return *this;
        }

        /// Pointwise division by \a q, where division by 0 yields 0
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see divide(const TProb<T> &)
         */
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

        /// Pointwise division by \a q, where division by 0 yields +Inf
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see operator/=(const TProb<T> &)
         */
        TProb<T>& divide (const TProb<T> & q) {
            DAI_DEBASSERT( size() == q.size() );
            std::transform( _p.begin(), _p.end(), q._p.begin(), _p.begin(), std::divides<T>() );
            return *this;
        }
    //@}

    /// \name Transformations with other equally-sized vectors
    //@{
        /// Returns sum of \c *this and \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T> operator+ (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> sum( *this );
            sum += q;
            return sum;
        }

        /// Return \c *this minus \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T> operator- (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> diff( *this );
            diff -= q;
            return diff;
        }

        /// Return product of \c *this with \a q
        /** \pre <tt>this->size() == q.size()</tt>
         */
        TProb<T> operator* (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> prod( *this );
            prod *= q;
            return prod;
        }

        /// Returns quotient of \c *this with \a q, where division by 0 yields +Inf
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see divided_by(const TProb<T> &)
         */
        TProb<T> operator/ (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> quot( *this );
            quot /= q;
            return quot;
        }

        /// Pointwise division by \a q, where division by 0 yields 0
        /** \pre <tt>this->size() == q.size()</tt>
         *  \see operator/(const TProb<T> &)
         */
        TProb<T> divided_by (const TProb<T> & q) const {
            DAI_DEBASSERT( size() == q.size() );
            TProb<T> quot( *this );
            quot.divide(q);
            return quot;
        }
    //@}
};


/// Returns distance between \a p and \a q, measured using distance measure \a dt
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> T dist( const TProb<T> &p, const TProb<T> &q, typename TProb<T>::DistType dt ) {
    DAI_DEBASSERT( p.size() == q.size() );
    T result = 0;
    switch( dt ) {
        case TProb<T>::DISTL1:
            for( size_t i = 0; i < p.size(); i++ )
                result += abs(p[i] - q[i]);
            break;

        case TProb<T>::DISTLINF:
            for( size_t i = 0; i < p.size(); i++ ) {
                T z = abs(p[i] - q[i]);
                if( z > result )
                    result = z;
            }
            break;

        case TProb<T>::DISTTV:
            for( size_t i = 0; i < p.size(); i++ )
                result += abs(p[i] - q[i]);
            result /= 2;
            break;

        case TProb<T>::DISTKL:
            for( size_t i = 0; i < p.size(); i++ ) {
                if( p[i] != 0.0 )
                    result += p[i] * (dai::log(p[i]) - dai::log(q[i]));
            }
    }
    return result;
}


/// Writes a TProb<T> to an output stream
/** \relates TProb
 */
template<typename T> std::ostream& operator<< (std::ostream& os, const TProb<T>& p) {
    os << "[";
    std::copy( p.p().begin(), p.p().end(), std::ostream_iterator<T>(os, " ") );
    os << "]";
    return os;
}


/// Returns the pointwise minimum of \a a and \a b
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> TProb<T> min( const TProb<T> &a, const TProb<T> &b ) {
    DAI_ASSERT( a.size() == b.size() );
    TProb<T> result( a.size() );
    for( size_t i = 0; i < a.size(); i++ )
        if( a[i] < b[i] )
            result[i] = a[i];
        else
            result[i] = b[i];
    return result;
}


/// Returns the pointwise maximum of \a a and \a b
/** \relates TProb
 *  \pre <tt>this->size() == q.size()</tt>
 */
template<typename T> TProb<T> max( const TProb<T> &a, const TProb<T> &b ) {
    DAI_ASSERT( a.size() == b.size() );
    TProb<T> result( a.size() );
    for( size_t i = 0; i < a.size(); i++ )
        if( a[i] > b[i] )
            result[i] = a[i];
        else
            result[i] = b[i];
    return result;
}


/// Represents a vector with entries of type dai::Real.
typedef TProb<Real> Prob;


} // end of namespace dai


#endif
