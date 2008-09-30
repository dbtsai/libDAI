/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

    Copyright (C) 2002  Martijn Leisink  [martijn@mbfys.kun.nl]
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


/// \file
/// \brief Defines TFactor<T> and Factor classes


#ifndef __defined_libdai_factor_h
#define __defined_libdai_factor_h


#include <iostream>
#include <cmath>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>


namespace dai {


// predefine TFactor<T> class
template<typename T> class      TFactor;


/// Represents a factor with probability entries represented as Real
typedef TFactor<Real>           Factor;


/// Represents a probability factor.
/** A \e factor is a function of the Cartesian product of the state
 *  spaces of some set of variables to the nonnegative real numbers.
 *  More formally, if \f$x_i \in X_i\f$ for all \f$i\f$, then a factor
 *  depending on the variables \f$\{x_i\}\f$ is a function defined
 *  on \f$\prod_i X_i\f$ with values in \f$[0,\infty)\f$.
 *  
 *  A Factor has two components: a VarSet, defining the set of variables
 *  that the factor depends on, and a TProb<T>, containing the values of
 *  the factor for all possible joint states of the variables.
 *
 *  \tparam T Should be castable from and to double.
 */
template <typename T> class TFactor {
    private:
        VarSet      _vs;
        TProb<T>    _p;

    public:
        /// Construct Factor with empty VarSet
        TFactor ( Real p = 1.0 ) : _vs(), _p(1,p) {}

        /// Construct Factor from VarSet
        TFactor( const VarSet& ns ) : _vs(ns), _p(_vs.nrStates()) {}
        
        /// Construct Factor from VarSet and initial value
        TFactor( const VarSet& ns, Real p ) : _vs(ns), _p(_vs.nrStates(),p) {}
        
        /// Construct Factor from VarSet and initial array
        TFactor( const VarSet& ns, const Real *p ) : _vs(ns), _p(_vs.nrStates(),p) {}

        /// Construct Factor from VarSet and TProb<T>
        TFactor( const VarSet& ns, const TProb<T>& p ) : _vs(ns), _p(p) {
#ifdef DAI_DEBUG
            assert( _vs.nrStates() == _p.size() );
#endif
        }
        
        /// Construct Factor from Var
        TFactor( const Var& n ) : _vs(n), _p(n.states()) {}

        /// Copy constructor
        TFactor( const TFactor<T> &x ) : _vs(x._vs), _p(x._p) {}
        
        /// Assignment operator
        TFactor<T> & operator= (const TFactor<T> &x) {
            if( this != &x ) {
                _vs = x._vs;
                _p  = x._p;
            }
            return *this;
        }

        /// Returns const reference to probability entries
        const TProb<T> & p() const { return _p; }
        /// Returns reference to probability entries
        TProb<T> & p() { return _p; }

        /// Returns const reference to variables
        const VarSet & vars() const { return _vs; }

        /// Returns the number of possible joint states of the variables
        size_t states() const { return _p.size(); }

        /// Returns a copy of the i'th probability value
        T operator[] (size_t i) const { return _p[i]; }

        /// Returns a reference to the i'th probability value
        T& operator[] (size_t i) { return _p[i]; }

        /// Sets all probability entries to p
        TFactor<T> & fill (T p) { _p.fill( p ); return(*this); }

        /// Fills all probability entries with random values
        TFactor<T> & randomize () { _p.randomize(); return(*this); }

        /// Returns product of *this with x
        TFactor<T> operator* (T x) const {
            Factor result = *this;
            result.p() *= x;
            return result;
        }

        /// Multiplies each probability entry with x
        TFactor<T>& operator*= (T x) {
            _p *= x;
            return *this;
        }

        /// Returns quotient of *this with x
        TFactor<T> operator/ (T x) const {
            Factor result = *this;
            result.p() /= x;
            return result;
        }

        /// Divides each probability entry by x
        TFactor<T>& operator/= (T x) {
            _p /= x;
            return *this;
        }

        /// Returns product of *this with another Factor
        TFactor<T> operator* (const TFactor<T>& Q) const;

        /// Returns quotient of *this with another Factor
        TFactor<T> operator/ (const TFactor<T>& Q) const;

        /// Multiplies *this with another Factor
        TFactor<T>& operator*= (const TFactor<T>& Q) { return( *this = (*this * Q) ); }

        /// Divides *this by another Factor
        TFactor<T>& operator/= (const TFactor<T>& Q) { return( *this = (*this / Q) ); }

        /// Returns sum of *this and another Factor (their vars() should be identical)
        TFactor<T> operator+ (const TFactor<T>& Q) const {
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            TFactor<T> sum(*this); 
            sum._p += Q._p; 
            return sum; 
        }

        /// Returns difference of *this and another Factor (their vars() should be identical)
        TFactor<T> operator- (const TFactor<T>& Q) const {
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            TFactor<T> sum(*this); 
            sum._p -= Q._p; 
            return sum; 
        }

        /// Adds another Factor to *this (their vars() should be identical)
        TFactor<T>& operator+= (const TFactor<T>& Q) { 
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            _p += Q._p;
            return *this;
        }

        /// Subtracts another Factor from *this (their vars() should be identical)
        TFactor<T>& operator-= (const TFactor<T>& Q) { 
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            _p -= Q._p;
            return *this;
        }

        /// Adds scalar to *this
        TFactor<T>& operator+= (T q) { 
            _p += q;
            return *this;
        }

        /// Subtracts scalar from *this
        TFactor<T>& operator-= (T q) { 
            _p -= q;
            return *this;
        }

        /// Returns sum of *this and a scalar
        TFactor<T> operator+ (T q) const {
            TFactor<T> result(*this); 
            result._p += q; 
            return result; 
        }

        /// Returns difference of *this with a scalar
        TFactor<T> operator- (T q) const {
            TFactor<T> result(*this); 
            result._p -= q; 
            return result; 
        }

        /// Returns *this raised to some power
        TFactor<T> operator^ (Real a) const { TFactor<T> x; x._vs = _vs; x._p = _p^a; return x; }

        /// Raises *this to some power
        TFactor<T>& operator^= (Real a) { _p ^= a; return *this; }

        /// Sets all entries that are smaller than epsilon to zero
        TFactor<T>& makeZero( Real epsilon ) {
            _p.makeZero( epsilon );
            return *this;
        }

        /// Sets all entries that are smaller than epsilon to epsilon
        TFactor<T>& makePositive( Real epsilon ) {
            _p.makePositive( epsilon );
            return *this;
        }
            
        /// Returns inverse of *this
        TFactor<T> inverse() const { 
            TFactor<T> inv; 
            inv._vs = _vs; 
            inv._p = _p.inverse(true);  // FIXME
            return inv; 
        }

        /// Returns *this divided by another Factor
        TFactor<T> divided_by( const TFactor<T>& denom ) const { 
#ifdef DAI_DEBUG
            assert( denom._vs == _vs );
#endif
            TFactor<T> quot(*this); 
            quot._p /= denom._p; 
            return quot; 
        }

        /// Divides *this by another Factor
        TFactor<T>& divide( const TFactor<T>& denom ) {
#ifdef DAI_DEBUG
            assert( denom._vs == _vs );
#endif
            _p /= denom._p;
            return *this;
        }

        /// Returns exp of *this
        TFactor<T> exp() const { 
            TFactor<T> e; 
            e._vs = _vs; 
            e._p = _p.exp(); 
            return e; 
        }

        /// Returns absolute value of *this
        TFactor<T> abs() const { 
            TFactor<T> e; 
            e._vs = _vs; 
            e._p = _p.abs(); 
            return e; 
        }

        /// Returns logarithm of *this
        TFactor<T> log() const {
            TFactor<T> l; 
            l._vs = _vs; 
            l._p = _p.log(); 
            return l; 
        }

        /// Returns logarithm of *this (defining log(0)=0)
        TFactor<T> log0() const {
            TFactor<T> l0; 
            l0._vs = _vs; 
            l0._p = _p.log0(); 
            return l0; 
        }

        /// Normalizes *this Factor
        T normalize( typename Prob::NormType norm = Prob::NORMPROB ) { return _p.normalize( norm ); }

        /// Returns a normalized copy of *this
        TFactor<T> normalized( typename Prob::NormType norm = Prob::NORMPROB ) const { 
            TFactor<T> result;
            result._vs = _vs;
            result._p = _p.normalized( norm );
            return result;
        }

        /// Returns a slice of this factor, where the subset ns is in state ns_state
        Factor slice( const VarSet & ns, size_t ns_state ) const {
            assert( ns << _vs );
            VarSet nsrem = _vs / ns;
            Factor result( nsrem, 0.0 );
            
            // OPTIMIZE ME
            IndexFor i_ns (ns, _vs);
            IndexFor i_nsrem (nsrem, _vs);
            for( size_t i = 0; i < states(); i++, ++i_ns, ++i_nsrem )
                if( (size_t)i_ns == ns_state )
                    result._p[i_nsrem] = _p[i];

            return result;
        }

        /// Returns unnormalized marginal; ns should be a subset of vars()
        TFactor<T> partSum(const VarSet & ns) const;

        /// Returns (normalized by default) marginal; ns should be a subset of vars()
        TFactor<T> marginal(const VarSet & ns, bool normed = true) const { if(normed) return partSum(ns).normalized(); else return partSum(ns); }

        /// Sums out all variables except those in ns
        TFactor<T> notSum(const VarSet & ns) const { return partSum(vars() ^ ns); }

        /// Embeds this factor in a larger VarSet
        TFactor<T> embed(const VarSet & ns) const { 
            VarSet vs = vars();
            assert( ns >> vs );
            if( vs == ns )
                return *this;
            else
                return (*this) * Factor(ns / vs, 1.0);
        }

        /// Returns true if *this has NANs
        bool hasNaNs() const { return _p.hasNaNs(); }

        /// Returns true if *this has negative entries
        bool hasNegatives() const { return _p.hasNegatives(); }

        /// Returns total sum of probability entries
        T totalSum() const { return _p.totalSum(); }

        /// Returns maximum absolute value of probability entries
        T maxAbs() const { return _p.maxAbs(); }

        /// Returns maximum value of probability entries
        T maxVal() const { return _p.maxVal(); }

        /// Returns minimum value of probability entries
        T minVal() const { return _p.minVal(); }

        /// Returns entropy of *this
        Real entropy() const { return _p.entropy(); }

        /// Returns strength of *this, between variables i and j, using (52) of [\ref MoK07b]
        T strength( const Var &i, const Var &j ) const;
};


template<typename T> TFactor<T> TFactor<T>::partSum(const VarSet & ns) const {
#ifdef DAI_DEBUG
    assert( ns << _vs );
#endif

    TFactor<T> res( ns, 0.0 );

    IndexFor i_res( ns, _vs );
    for( size_t i = 0; i < _p.size(); i++, ++i_res )
        res._p[i_res] += _p[i];

    return res;
}


template<typename T> TFactor<T> TFactor<T>::operator* (const TFactor<T>& Q) const {
    TFactor<T> prod( _vs | Q._vs, 0.0 );

    IndexFor i1(_vs, prod._vs);
    IndexFor i2(Q._vs, prod._vs);

    for( size_t i = 0; i < prod._p.size(); i++, ++i1, ++i2 )
        prod._p[i] += _p[i1] * Q._p[i2];

    return prod;
}


template<typename T> TFactor<T> TFactor<T>::operator/ (const TFactor<T>& Q) const {
    TFactor<T> quot( _vs + Q._vs, 0.0 );

    IndexFor i1(_vs, quot._vs);
    IndexFor i2(Q._vs, quot._vs);

    for( size_t i = 0; i < quot._p.size(); i++, ++i1, ++i2 )
        quot._p[i] += _p[i1] / Q._p[i2];

    return quot;
}


template<typename T> T TFactor<T>::strength( const Var &i, const Var &j ) const {
#ifdef DAI_DEBUG
    assert( _vs.contains( i ) );
    assert( _vs.contains( j ) );
    assert( i != j );
#endif
    VarSet ij(i, j);

    T max = 0.0;
    for( size_t alpha1 = 0; alpha1 < i.states(); alpha1++ )
        for( size_t alpha2 = 0; alpha2 < i.states(); alpha2++ )
            if( alpha2 != alpha1 )
                for( size_t beta1 = 0; beta1 < j.states(); beta1++ ) 
                    for( size_t beta2 = 0; beta2 < j.states(); beta2++ )
                        if( beta2 != beta1 ) {
                            size_t as = 1, bs = 1;
                            if( i < j )
                                bs = i.states();
                            else
                                as = j.states();
                            T f1 = slice( ij, alpha1 * as + beta1 * bs ).p().divide( slice( ij, alpha2 * as + beta1 * bs ).p() ).maxVal();
                            T f2 = slice( ij, alpha2 * as + beta2 * bs ).p().divide( slice( ij, alpha1 * as + beta2 * bs ).p() ).maxVal();
                            T f = f1 * f2;
                            if( f > max )
                                max = f;
                        }
    
    return std::tanh( 0.25 * std::log( max ) );
}


/// Writes a Factor to an output stream
template<typename T> std::ostream& operator<< (std::ostream& os, const TFactor<T>& P) {
    os << "(" << P.vars() << " <";
    for( size_t i = 0; i < P.states(); i++ )
        os << P[i] << " ";
    os << ">)";
    return os;
}


/// Returns distance between two Factors (with identical vars())
template<typename T> Real dist( const TFactor<T> & x, const TFactor<T> & y, Prob::DistType dt ) {
    if( x.vars().empty() || y.vars().empty() )
        return -1;
    else {
#ifdef DAI_DEBUG
        assert( x.vars() == y.vars() );
#endif
        return dist( x.p(), y.p(), dt );
    }
}


/// Returns the pointwise maximum of two Factors
template<typename T> TFactor<T> max( const TFactor<T> & P, const TFactor<T> & Q ) {
    assert( P._vs == Q._vs );
    return TFactor<T>( P._vs, min( P.p(), Q.p() ) );
}


/// Returns the pointwise minimum of two Factors
template<typename T> TFactor<T> min( const TFactor<T> & P, const TFactor<T> & Q ) {
    assert( P._vs == Q._vs );
    return TFactor<T>( P._vs, max( P.p(), Q.p() ) );
}


/// Calculates the mutual information between the two variables in P
template<typename T> Real MutualInfo(const TFactor<T> & P) {
    assert( P.vars().size() == 2 );
    VarSet::const_iterator it = P.vars().begin();
    Var i = *it; it++; Var j = *it;
    TFactor<T> projection = P.marginal(i) * P.marginal(j);
    return real( dist( P.normalized(), projection, Prob::DISTKL ) );
}


} // end of namespace dai


#endif
