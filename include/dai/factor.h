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
#include <functional>
#include <cmath>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>
#include <dai/util.h>


namespace dai {


// Function object similar to std::divides(), but different in that dividing by zero results in zero
template<typename T> struct divides0 : public std::binary_function<T, T, T> {
    // Returns (j == 0 ? 0 : (i/j))
    T operator()(const T& i, const T& j) const {
        if( j == (T)0 )
            return (T)0;
        else
            return i / j;
    }
};


/// Represents a (probability) factor.
/** Mathematically, a \e factor is a function mapping joint states of some
 *  variables to the nonnegative real numbers.
 *  More formally, denoting a discrete variable with label \f$l\f$ by
 *  \f$x_l\f$ and its state space by \f$X_l = \{0,1,\dots,S_l-1\}\f$,
 *  then a factor depending on the variables \f$\{x_l\}_{l\in L}\f$ is 
 *  a function \f$f_L : \prod_{l\in L} X_l \to [0,\infty)\f$.
 *  
 *  In libDAI, a factor is represented by a TFactor<\a T> object, which has two
 *  components: 
 *  \arg a VarSet, corresponding with the set of variables \f$\{x_l\}_{l\in L}\f$ 
 *  that the factor depends on;
 *  \arg a TProb<\a T>, a vector containing the value of the factor for each possible 
 *  joint state of the variables.
 *
 *  The factor values are stored in the entries of the TProb<\a T> in a particular
 *  ordering, which is defined by the one-to-one correspondence of a joint state 
 *  in \f$\prod_{l\in L} X_l\f$ with a linear index in 
 *  \f$\{0,1,\dots,\prod_{l\in L} S_l-1\}\f$ according to the mapping \f$\sigma\f$
 *  induced by VarSet::calcState(const std::map<Var,size_t> &).
 *
 *  \tparam T Should be a scalar that is castable from and to double and should support elementary arithmetic operations.
 *  \todo Define a better fileformat for .fg files (maybe using XML)?
 *  \todo Add support for sparse factors.
 */
template <typename T> class TFactor {
    private:
        VarSet      _vs;
        TProb<T>    _p;

    public:
        /// Iterator over factor entries
		typedef typename TProb<T>::iterator iterator;

        /// Const iterator over factor entries
		typedef typename TProb<T>::const_iterator const_iterator;

        /// Constructs TFactor depending on no variables, with value p
        TFactor ( Real p = 1.0 ) : _vs(), _p(1,p) {}

        /// Constructs TFactor depending on variables in ns, with uniform distribution
        TFactor( const VarSet& ns ) : _vs(ns), _p(_vs.nrStates()) {}
        
        /// Constructs TFactor depending on variables in ns, with all values set to p
        TFactor( const VarSet& ns, Real p ) : _vs(ns), _p(_vs.nrStates(),p) {}
        
        /// Constructs TFactor depending on variables in ns, copying the values from the range starting at begin
        /** \param ns contains the variables that the new TFactor should depend on.
         *  \tparam Iterator Iterates over instances of type T; should support addition of size_t.
         *  \param begin Points to first element to be added.
         */
        template<typename TIterator>
        TFactor( const VarSet& ns, TIterator begin ) : _vs(ns), _p(begin, begin + _vs.nrStates(), _vs.nrStates()) {}

        /// Constructs TFactor depending on variables in ns, with values set to the TProb p
        TFactor( const VarSet& ns, const TProb<T>& p ) : _vs(ns), _p(p) {
            DAI_DEBASSERT( _vs.nrStates() == _p.size() );
        }
        
        /// Constructs TFactor depending on the variable n, with uniform distribution
        TFactor( const Var& n ) : _vs(n), _p(n.states()) {}

        /// Returns const reference to value vector
        const TProb<T> & p() const { return _p; }
        /// Returns reference to value vector
        TProb<T> & p() { return _p; }

        /// Returns const reference to variable set
        const VarSet & vars() const { return _vs; }

        /// Returns the number of possible joint states of the variables
        /** \note This is equal to the length of the value vector.
         */
        size_t states() const { return _p.size(); }

        /// Returns a copy of the i'th entry of the value vector
        T operator[] (size_t i) const { return _p[i]; }

        /// Returns a reference to the i'th entry of the value vector
        T& operator[] (size_t i) { return _p[i]; }
        
        /// Returns iterator pointing to first entry
        iterator begin() { return _p.begin(); }
        /// Returns const iterator pointing to first entry
		const_iterator begin() const { return _p.begin(); }
		/// Returns iterator pointing beyond last entry
		iterator end() { return _p.end(); }
		/// Returns const iterator pointing beyond last entry
		const_iterator end() const { return _p.end(); }

        /// Sets all values to p
        TFactor<T> & fill (T p) { _p.fill( p ); return(*this); }

        /// Draws all values i.i.d. from a uniform distribution on [0,1)
        TFactor<T> & randomize () { _p.randomize(); return(*this); }


        /// Multiplies *this with scalar t
        TFactor<T>& operator*= (T t) {
            _p *= t;
            return *this;
        }

        /// Divides *this by scalar t
        TFactor<T>& operator/= (T t) {
            _p /= t;
            return *this;
        }

        /// Adds scalar t to *this
        TFactor<T>& operator+= (T t) { 
            _p += t;
            return *this;
        }

        /// Subtracts scalar t from *this
        TFactor<T>& operator-= (T t) { 
            _p -= t;
            return *this;
        }

        /// Raises *this to the power a
        TFactor<T>& operator^= (Real a) { _p ^= a; return *this; }


        /// Returns product of *this with scalar t
        TFactor<T> operator* (T t) const {
            TFactor<T> result = *this;
            result.p() *= t;
            return result;
        }

        /// Returns quotient of *this with scalar t
        TFactor<T> operator/ (T t) const {
            TFactor<T> result = *this;
            result.p() /= t;
            return result;
        }

        /// Returns sum of *this and scalar t
        TFactor<T> operator+ (T t) const {
            TFactor<T> result(*this); 
            result._p += t; 
            return result; 
        }

        /// Returns *this minus scalar t
        TFactor<T> operator- (T t) const {
            TFactor<T> result(*this); 
            result._p -= t; 
            return result; 
        }

        /// Returns *this raised to the power a
        TFactor<T> operator^ (Real a) const { 
            TFactor<T> x; 
            x._vs = _vs; 
            x._p = _p^a; 
            return x; 
        }

        /// Multiplies *this with the TFactor f
        TFactor<T>& operator*= (const TFactor<T>& f) { 
            if( f._vs == _vs ) // optimize special case
                _p *= f._p;
            else
                *this = (*this * f); 
            return *this;
        }

        /// Divides *this by the TFactor f
        TFactor<T>& operator/= (const TFactor<T>& f) { 
            if( f._vs == _vs ) // optimize special case
                _p /= f._p;
            else
                *this = (*this / f); 
            return *this;
        }

        /// Adds the TFactor f to *this
        TFactor<T>& operator+= (const TFactor<T>& f) { 
            if( f._vs == _vs ) // optimize special case
                _p += f._p;
            else
                *this = (*this + f); 
            return *this;
        }

        /// Subtracts the TFactor f from *this
        TFactor<T>& operator-= (const TFactor<T>& f) { 
            if( f._vs == _vs ) // optimize special case
                _p -= f._p;
            else
                *this = (*this - f); 
            return *this;
        }

        /// Returns product of *this with the TFactor f
        /** The product of two factors is defined as follows: if 
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[fg : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) g(x_M).\f]
         */
        TFactor<T> operator* (const TFactor<T>& f) const {
            return pointwiseOp(*this,f,std::multiplies<T>());
        }

        /// Returns quotient of *this by the TFactor f
        /** The quotient of two factors is defined as follows: if 
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[\frac{f}{g} : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto \frac{f(x_L)}{g(x_M)}.\f]
         */
        TFactor<T> operator/ (const TFactor<T>& f) const {
            return pointwiseOp(*this,f,divides0<T>());
        }

        /// Returns sum of *this and the TFactor f
        /** The sum of two factors is defined as follows: if 
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f+g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) + g(x_M).\f]
         */
        TFactor<T> operator+ (const TFactor<T>& f) const {
            return pointwiseOp(*this,f,std::plus<T>());
        }

        /// Returns *this minus the TFactor f
        /** The difference of two factors is defined as follows: if 
         *  \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$g : \prod_{m\in M} X_m \to [0,\infty)\f$, then
         *  \f[f-g : \prod_{l\in L\cup M} X_l \to [0,\infty) : x \mapsto f(x_L) - g(x_M).\f]
         */
        TFactor<T> operator- (const TFactor<T>& f) const {
            return pointwiseOp(*this,f,std::minus<T>());
        }


        /// Sets all values that are smaller than epsilon to 0
        TFactor<T>& makeZero( T epsilon ) {
            _p.makeZero( epsilon );
            return *this;
        }

        /// Sets all values that are smaller than epsilon to epsilon
        TFactor<T>& makePositive( T epsilon ) {
            _p.makePositive( epsilon );
            return *this;
        }
            
        /// Returns pointwise inverse of *this.
        /** If zero == true, uses 1 / 0 == 0; otherwise 1 / 0 == Inf.
         */
        TFactor<T> inverse(bool zero=true) const { 
            TFactor<T> inv; 
            inv._vs = _vs; 
            inv._p = _p.inverse(zero);
            return inv; 
        }

        /// Returns pointwise exp of *this
        TFactor<T> exp() const { 
            TFactor<T> e; 
            e._vs = _vs; 
            e._p = _p.exp(); 
            return e; 
        }

        /// Returns pointwise logarithm of *this
        /** If zero==true, uses log(0)==0; otherwise, log(0)=-Inf.
         */
        TFactor<T> log(bool zero=false) const {
            TFactor<T> l; 
            l._vs = _vs; 
            l._p = _p.log(zero); 
            return l; 
        }

        /// Returns pointwise absolute value of *this
        TFactor<T> abs() const { 
            TFactor<T> e; 
            e._vs = _vs; 
            e._p = _p.abs(); 
            return e; 
        }

        /// Normalizes *this TFactor according to the specified norm
        T normalize( typename Prob::NormType norm=Prob::NORMPROB ) { return _p.normalize( norm ); }

        /// Returns a normalized copy of *this, according to the specified norm
        TFactor<T> normalized( typename Prob::NormType norm=Prob::NORMPROB ) const { 
            TFactor<T> result;
            result._vs = _vs;
            result._p = _p.normalized( norm );
            return result;
        }

        /// Returns a slice of this TFactor, where the subset ns is in state nsState
        /** \pre \a ns sould be a subset of vars()
         *  \pre \a nsState < ns.states()
         *
         *  The result is a TFactor that depends on the variables in this->vars() except those in \a ns,
         *  obtained by setting the variables in \a ns to the joint state specified by the linear index 
         *  \a nsState. Formally, if *this corresponds with the factor \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$,
         *  \f$M \subset L\f$ corresponds with \a ns and \a nsState corresponds with a mapping \f$s\f$ that
         *  maps a variable \f$x_m\f$ with \f$m\in M\f$ to its state \f$s(x_m) \in X_m\f$, then the slice 
         *  returned corresponds with the factor \f$g : \prod_{l \in L \setminus M} X_l \to [0,\infty)\f$
         *  defined by \f$g(\{x_l\}_{l\in L \setminus M}) = f(\{x_l\}_{l\in L \setminus M}, \{s(x_m)\}_{m\in M})\f$.
         */
        TFactor<T> slice( const VarSet& ns, size_t nsState ) const {
            assert( ns << _vs );
            VarSet nsrem = _vs / ns;
            TFactor<T> result( nsrem, T(0) );
            
            // OPTIMIZE ME
            IndexFor i_ns (ns, _vs);
            IndexFor i_nsrem (nsrem, _vs);
            for( size_t i = 0; i < states(); i++, ++i_ns, ++i_nsrem )
                if( (size_t)i_ns == nsState )
                    result._p[i_nsrem] = _p[i];

            return result;
        }

        /// Returns marginal on ns, obtained by summing out all variables except those in ns, and normalizing the result if normed==true
        TFactor<T> marginal(const VarSet & ns, bool normed=true) const;

        /// Returns max-marginal on ns, obtained by maximizing all variables except those in ns, and normalizing the result if normed==true
        TFactor<T> maxMarginal(const VarSet & ns, bool normed=true) const;

        /// Embeds this factor in a larger VarSet
        /** \pre vars() should be a subset of ns
         *
         *  If *this corresponds with \f$f : \prod_{l\in L} X_l \to [0,\infty)\f$ and \f$L \subset M\f$, then
         *  the embedded factor corresponds with \f$g : \prod_{m\in M} X_m \to [0,\infty) : x \mapsto f(x_L)\f$.
         */
        TFactor<T> embed(const VarSet & ns) const { 
            assert( ns >> _vs );
            if( _vs == ns )
                return *this;
            else
                return (*this) * TFactor<T>(ns / _vs, (T)1);
        }

        /// Returns true if *this has NaN values
        bool hasNaNs() const { return _p.hasNaNs(); }

        /// Returns true if *this has negative values
        bool hasNegatives() const { return _p.hasNegatives(); }

        /// Returns total sum of values
        T sum() const { return _p.sum(); }

        /// Returns maximum absolute value
        T maxAbs() const { return _p.maxAbs(); }

        /// Returns maximum value
        T max() const { return _p.max(); }

        /// Returns minimum value
        T min() const { return _p.min(); }

        /// Returns entropy of *this TFactor
        Real entropy() const { return _p.entropy(); }

        /// Returns strength of *this TFactor (between variables i and j), as defined in eq. (52) of [\ref MoK07b]
        T strength( const Var &i, const Var &j ) const;
};


template<typename T> TFactor<T> TFactor<T>::marginal(const VarSet & ns, bool normed) const {
    VarSet res_ns = ns & _vs;

    TFactor<T> res( res_ns, 0.0 );

    IndexFor i_res( res_ns, _vs );
    for( size_t i = 0; i < _p.size(); i++, ++i_res )
        res._p[i_res] += _p[i];

    if( normed )
        res.normalize( Prob::NORMPROB );

    return res;
}


template<typename T> TFactor<T> TFactor<T>::maxMarginal(const VarSet & ns, bool normed) const {
    VarSet res_ns = ns & _vs;

    TFactor<T> res( res_ns, 0.0 );

    IndexFor i_res( res_ns, _vs );
    for( size_t i = 0; i < _p.size(); i++, ++i_res )
        if( _p[i] > res._p[i_res] )
            res._p[i_res] = _p[i];

    if( normed )
        res.normalize( Prob::NORMPROB );

    return res;
}


/// Apply binary operator pointwise on two factors
template<typename T, typename binaryOp> TFactor<T> pointwiseOp( const TFactor<T> &f, const TFactor<T> &g, binaryOp op ) {
    if( f.vars() == g.vars() ) { // optimizate special case
        TFactor<T> result(f); 
        for( size_t i = 0; i < result.states(); i++ )
            result[i] = op( result[i], g[i] );
        return result; 
    } else {
        TFactor<T> result( f.vars() | g.vars(), 0.0 );

        IndexFor i1(f.vars(), result.vars());
        IndexFor i2(g.vars(), result.vars());

        for( size_t i = 0; i < result.states(); i++, ++i1, ++i2 )
            result[i] = op( f[i1], g[i2] );

        return result;
    }
}


template<typename T> T TFactor<T>::strength( const Var &i, const Var &j ) const {
    DAI_DEBASSERT( _vs.contains( i ) );
    DAI_DEBASSERT( _vs.contains( j ) );
    DAI_DEBASSERT( i != j );
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
                            T f1 = slice( ij, alpha1 * as + beta1 * bs ).p().divide( slice( ij, alpha2 * as + beta1 * bs ).p() ).max();
                            T f2 = slice( ij, alpha2 * as + beta2 * bs ).p().divide( slice( ij, alpha1 * as + beta2 * bs ).p() ).max();
                            T f = f1 * f2;
                            if( f > max )
                                max = f;
                        }
    
    return std::tanh( 0.25 * std::log( max ) );
}


/// Writes a TFactor to an output stream
/** \relates TFactor
 */
template<typename T> std::ostream& operator<< (std::ostream& os, const TFactor<T>& P) {
    os << "(" << P.vars() << ", (";
    for( size_t i = 0; i < P.states(); i++ )
        os << (i == 0 ? "" : ", ") << P[i];
    os << "))";
    return os;
}


/// Returns distance between two TFactors f and g, according to the distance measure dt
/** \relates TFactor
 *  \pre f.vars() == g.vars()
 */
template<typename T> Real dist( const TFactor<T> &f, const TFactor<T> &g, Prob::DistType dt ) {
    if( f.vars().empty() || g.vars().empty() )
        return -1;
    else {
        DAI_DEBASSERT( f.vars() == g.vars() );
        return dist( f.p(), g.p(), dt );
    }
}


/// Returns the pointwise maximum of two TFactors
/** \relates TFactor
 *  \pre f.vars() == g.vars()
 */
template<typename T> TFactor<T> max( const TFactor<T> &f, const TFactor<T> &g ) {
    assert( f._vs == g._vs );
    return TFactor<T>( f._vs, min( f.p(), g.p() ) );
}


/// Returns the pointwise minimum of two TFactors
/** \relates TFactor
 *  \pre f.vars() == g.vars()
 */
template<typename T> TFactor<T> min( const TFactor<T> &f, const TFactor<T> &g ) {
    assert( f._vs == g._vs );
    return TFactor<T>( f._vs, max( f.p(), g.p() ) );
}


/// Calculates the mutual information between the two variables that f depends on, under the distribution given by f
/** \relates TFactor
 *  \pre f.vars().size() == 2
 */ 
template<typename T> Real MutualInfo(const TFactor<T> &f) {
    assert( f.vars().size() == 2 );
    VarSet::const_iterator it = f.vars().begin();
    Var i = *it; it++; Var j = *it;
    TFactor<T> projection = f.marginal(i) * f.marginal(j);
    return dist( f.normalized(), projection, Prob::DISTKL );
}


/// Represents a factor with values of type Real.
typedef TFactor<Real> Factor;


} // end of namespace dai


#endif
