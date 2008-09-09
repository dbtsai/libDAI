/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
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


#ifndef __defined_libdai_factor_h
#define __defined_libdai_factor_h


#include <iostream>
#include <cmath>
#include <dai/prob.h>
#include <dai/varset.h>
#include <dai/index.h>


namespace dai {


template<typename T> class      TFactor;
typedef TFactor<Real>           Factor;
typedef TFactor<Complex>        CFactor;


// predefine friends
template<typename T> Real           dist( const TFactor<T> & x, const TFactor<T> & y, Prob::DistType dt );
template<typename T> Complex        KL_dist( const TFactor<T> & p, const TFactor<T> & q );
template<typename T> std::ostream&  operator<< (std::ostream& os, const TFactor<T>& P);

        
// T should be castable from and to double and to complex
template <typename T> class TFactor {
    protected:
        VarSet      _vs;
        TProb<T>    _p;

    public:
        // Default constructor
        TFactor () : _vs(), _p(1,1.0) {}
        
        // Construct Factor from VarSet
        TFactor( const VarSet& ns ) : _vs(ns), _p(_vs.states()) {}
        
        // Construct Factor from VarSet and initial value
        TFactor( const VarSet& ns, Real p ) : _vs(ns), _p(_vs.states(),p) {}
        
        // Construct Factor from VarSet and initial array
        TFactor( const VarSet& ns, const Real* p ) : _vs(ns), _p(_vs.states(),p) {}

        // Construct Factor from VarSet and TProb<T>
        TFactor( const VarSet& ns, const TProb<T> p ) : _vs(ns), _p(p) {
#ifdef DAI_DEBUG
            assert( _vs.states() == _p.size() );
#endif
        }
        
        // Construct Factor from Var
        TFactor( const Var& n ) : _vs(n), _p(n.states()) {}

        // Copy constructor
        TFactor( const TFactor<T> &x ) : _vs(x._vs), _p(x._p) {}
        
        // Assignment operator
        TFactor<T> & operator= (const TFactor<T> &x) {
            if( this != &x ) {
                _vs = x._vs;
                _p  = x._p;
            }
            return *this;
        }

        const TProb<T> & p() const { return _p; }
        TProb<T> & p() { return _p; }
        const VarSet & vars() const { return _vs; }
        size_t states() const { 
#ifdef DAI_DEBUG
            assert( _vs.states() == _p.size() );
#endif
            return _p.size();
        }

        T operator[] (size_t i) const { return _p[i]; }
        T& operator[] (size_t i) { return _p[i]; }
        TFactor<T> & fill (T p)
            { _p.fill( p ); return(*this); }
        TFactor<T> & randomize ()
            { _p.randomize(); return(*this); }
        TFactor<T> operator* (T x) const {
            Factor result = *this;
            result.p() *= x;
            return result;
        }
        TFactor<T>& operator*= (T x) {
            _p *= x;
            return *this;
        }
        TFactor<T> operator/ (T x) const {
            Factor result = *this;
            result.p() /= x;
            return result;
        }
        TFactor<T>& operator/= (T x) {
            _p /= x;
            return *this;
        }
        TFactor<T> operator* (const TFactor<T>& Q) const;
        TFactor<T>& operator*= (const TFactor<T>& Q) { return( *this = (*this * Q) ); }
        TFactor<T> operator+ (const TFactor<T>& Q) const {
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            TFactor<T> sum(*this); 
            sum._p += Q._p; 
            return sum; 
        }
        TFactor<T> operator- (const TFactor<T>& Q) const {
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            TFactor<T> sum(*this); 
            sum._p -= Q._p; 
            return sum; 
        }
        TFactor<T>& operator+= (const TFactor<T>& Q) { 
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            _p += Q._p;
            return *this;
        }
        TFactor<T>& operator-= (const TFactor<T>& Q) { 
#ifdef DAI_DEBUG
            assert( Q._vs == _vs );
#endif
            _p -= Q._p;
            return *this;
        }

        TFactor<T> operator^ (Real a) const { TFactor<T> x; x._vs = _vs; x._p = _p^a; return x; }
        TFactor<T>& operator^= (Real a) { _p ^= a; return *this; }

        TFactor<T>& makeZero( Real epsilon ) {
            _p.makeZero( epsilon );
            return *this;
        }
            
        TFactor<T> inverse() const { 
            TFactor<T> inv; 
            inv._vs = _vs; 
            inv._p = _p.inverse(true);  // FIXME
            return inv; 
        }

        TFactor<T> divided_by( const TFactor<T>& denom ) const { 
#ifdef DAI_DEBUG
            assert( denom._vs == _vs );
#endif
            TFactor<T> quot(*this); 
            quot._p /= denom._p; 
            return quot; 
        }

        TFactor<T>& divide( const TFactor<T>& denom ) {
#ifdef DAI_DEBUG
            assert( denom._vs == _vs );
#endif
            _p /= denom._p;
            return *this;
        }

        TFactor<T> exp() const { 
            TFactor<T> e; 
            e._vs = _vs; 
            e._p = _p.exp(); 
            return e; 
        }

        TFactor<T> log() const {
            TFactor<T> l; 
            l._vs = _vs; 
            l._p = _p.log(); 
            return l; 
        }

        TFactor<T> log0() const {
            TFactor<T> l0; 
            l0._vs = _vs; 
            l0._p = _p.log0(); 
            return l0; 
        }

        CFactor clog0() const {
            CFactor l0; 
            l0._vs = _vs; 
            l0._p = _p.clog0(); 
            return l0; 
        }

        T normalize( typename Prob::NormType norm ) { return _p.normalize( norm ); }
        TFactor<T> normalized( typename Prob::NormType norm ) const { 
            TFactor<T> result;
            result._vs = _vs;
            result._p = _p.normalized( norm );
            return result;
        }

        // returns slice of this factor where the subset ns is in state ns_state
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

        // returns unnormalized marginal
        TFactor<T> part_sum(const VarSet & ns) const;
        // returns normalized marginal
        TFactor<T> marginal(const VarSet & ns) const { return part_sum(ns).normalized( Prob::NORMPROB ); }

        bool hasNaNs() const { return _p.hasNaNs(); }
        bool hasNegatives() const { return _p.hasNegatives(); }
        T totalSum() const { return _p.totalSum(); }
        T maxAbs() const { return _p.maxAbs(); }
        T maxVal() const { return _p.maxVal(); }
        Complex entropy() const { return _p.entropy(); }
        T strength( const Var &i, const Var &j ) const;

        friend Real dist( const TFactor<T> & x, const TFactor<T> & y, Prob::DistType dt ) {
            if( x._vs.empty() || y._vs.empty() )
                return -1;
            else {
#ifdef DAI_DEBUG
                assert( x._vs == y._vs );
#endif
                return dist( x._p, y._p, dt );
            }
        }
        friend Complex KL_dist <> (const TFactor<T> & p, const TFactor<T> & q);
        template<class U> friend std::ostream& operator<< (std::ostream& os, const TFactor<U>& P);
};


template<typename T> TFactor<T> TFactor<T>::part_sum(const VarSet & ns) const {
#ifdef DAI_DEBUG
    assert( ns << _vs );
#endif

    TFactor<T> res( ns, 0.0 );

    IndexFor i_res( ns, _vs );
    for( size_t i = 0; i < _p.size(); i++, ++i_res )
        res._p[i_res] += _p[i];

    return res;
}


template<typename T> std::ostream& operator<< (std::ostream& os, const TFactor<T>& P) {
    os << "(" << P.vars() << " <";
    for( size_t i = 0; i < P._p.size(); i++ )
        os << P._p[i] << " ";
    os << ">)";
    return os;
}


template<typename T> TFactor<T> TFactor<T>::operator* (const TFactor<T>& Q) const {
    TFactor<T> prod( _vs | Q._vs, 0.0 );

    IndexFor i1(_vs, prod._vs);
    IndexFor i2(Q._vs, prod._vs);

    for( size_t i = 0; i < prod._p.size(); i++, ++i1, ++i2 )
        prod._p[i] += _p[i1] * Q._p[i2];

    return prod;
}


template<typename T> Complex KL_dist(const TFactor<T> & P, const TFactor<T> & Q) {
    if( P._vs.empty() || Q._vs.empty() )
        return -1;
    else {
#ifdef DAI_DEBUG
        assert( P._vs == Q._vs );
#endif
        return KL_dist( P._p, Q._p );
    }
}


// calculate N(psi, i, j)
template<typename T> T TFactor<T>::strength( const Var &i, const Var &j ) const {
#ifdef DAI_DEBUG
    assert( _vs && i );
    assert( _vs && j );
    assert( i != j );
#endif
    VarSet ij = i | j;

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


template<typename T> TFactor<T> RemoveFirstOrderInteractions( const TFactor<T> & psi ) {
    TFactor<T> result = psi;

    VarSet vars = psi.vars();
    for( size_t iter = 0; iter < 100; iter++ ) {
        for( VarSet::const_iterator n = vars.begin(); n != vars.end(); n++ )
            result = result * result.part_sum(*n).inverse();
        result.normalize( Prob::NORMPROB );
    }

    return result;
}


} // end of namespace dai


#endif
