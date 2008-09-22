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


#ifndef __defined_libdai_varset_h
#define __defined_libdai_varset_h


#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <dai/var.h>
#include <dai/util.h>


namespace dai {


/// A VarSet represents a set of variables.
/**
 *  It is implemented as an ordered std::vector<Var> for efficiency reasons
 *  (indeed, it was found that a std::set<Var> usually has more overhead). 
 *  In addition, it provides an interface for common set-theoretic operations.
 */
class VarSet {
    private:
        /// The variables in this set
        std::vector<Var> _vars;

        /// Product of number of states of all contained variables
        size_t _states;

    public:
        /// Default constructor
        VarSet() : _vars(), _states(1) {};

        /// Construct a VarSet from one variable
        VarSet( const Var &n ) : _vars(), _states( n.states() ) { 
            _vars.push_back( n );
        }

        /// Construct a VarSet from two variables
        VarSet( const Var &n1, const Var &n2 ) { 
            if( n1 < n2 ) {
                _vars.push_back( n1 );
                _vars.push_back( n2 );
            } else if( n1 > n2 ) {
                _vars.push_back( n2 );
                _vars.push_back( n1 );
            } else
                _vars.push_back( n1 );
            calcStates();
        }

        /// Construct from a range of iterators
        /** The value_type of the VarIterator should be Var.
         *  For efficiency, the number of variables can be
         *  speficied by sizeHint.
         */
        template <typename VarIterator>
        VarSet( VarIterator begin, VarIterator end, size_t sizeHint=0 ) {
            _vars.reserve( sizeHint );
            _vars.insert( _vars.begin(), begin, end );
            std::sort( _vars.begin(), _vars.end() );
            std::vector<Var>::iterator new_end = std::unique( _vars.begin(), _vars.end() );
            _vars.erase( new_end, _vars.end() );
            calcStates();
        }

        /// Copy constructor
        VarSet( const VarSet &x ) : _vars( x._vars ), _states( x._states ) {}

        /// Assignment operator
        VarSet & operator=( const VarSet &x ) {
            if( this != &x ) {
                _vars = x._vars;
                _states = x._states;
            }
            return *this;
        }
        

        /// Returns the product of the number of states of each variable in this set
        size_t states() const { 
            return _states; 
        }
        

        /// Setminus operator (result contains all variables in *this, except those in ns)
        VarSet operator/ ( const VarSet& ns ) const {
            VarSet res;
            std::set_difference( _vars.begin(), _vars.end(), ns._vars.begin(), ns._vars.end(), inserter( res._vars, res._vars.begin() ) );
            res.calcStates();
            return res;
        }

        /// Set-union operator (result contains all variables in *this, plus those in ns)
        VarSet operator| ( const VarSet& ns ) const {
            VarSet res;
            std::set_union( _vars.begin(), _vars.end(), ns._vars.begin(), ns._vars.end(), inserter( res._vars, res._vars.begin() ) );
            res.calcStates();
            return res;
        }

        /// Set-intersection operator (result contains all variables in *this that are also contained in ns)
        VarSet operator& ( const VarSet& ns ) const {
            VarSet res;
            std::set_intersection( _vars.begin(), _vars.end(), ns._vars.begin(), ns._vars.end(), inserter( res._vars, res._vars.begin() ) );
            res.calcStates();
            return res;
        }
        
        /// Erases from *this all variables in ns
        VarSet& operator/= ( const VarSet& ns ) {
            return (*this = (*this / ns));
        }

        /// Erase one variable
        VarSet& operator/= ( const Var& n ) { 
            std::vector<Var>::iterator pos = lower_bound( _vars.begin(), _vars.end(), n );
            if( pos != _vars.end() )
                if( *pos == n ) { // found variable, delete it
                    _vars.erase( pos ); 
                    _states /= n.states();
                }
            return *this; 
        }

        /// Adds to *this all variables in ns
        VarSet& operator|= ( const VarSet& ns ) {
            return( *this = (*this | ns) );
        }

        /// Add one variable
        VarSet& operator|= ( const Var& n ) {
            std::vector<Var>::iterator pos = lower_bound( _vars.begin(), _vars.end(), n );
            if( pos == _vars.end() || *pos != n ) { // insert it
                _vars.insert( pos, n );
                _states *= n.states();
            }
            return *this;
        }


        /// Erases from *this all variables not in ns
        VarSet& operator&= ( const VarSet& ns ) { 
            return (*this = (*this & ns));
        }


        /// Returns true if *this is a subset of ns
        bool operator<< ( const VarSet& ns ) const { 
            return std::includes( ns._vars.begin(), ns._vars.end(), _vars.begin(), _vars.end() );
        }

        /// Returns true if ns is a subset of *this
        bool operator>> ( const VarSet& ns ) const { 
            return std::includes( _vars.begin(), _vars.end(), ns._vars.begin(), ns._vars.end() );
        }

        /// Returns true if *this and ns contain common variables
        bool intersects( const VarSet& ns ) const { 
            return( (*this & ns).size() > 0 ); 
        }

        /// Returns true if *this contains the variable n
        bool contains( const Var& n ) const { 
            return std::binary_search( _vars.begin(), _vars.end(), n );
        }

        /// Sends a VarSet to an output stream
        friend std::ostream& operator<< (std::ostream & os, const VarSet& ns) {
            foreach( const Var &n, ns._vars )
                os << n;
            return( os );
        }

        /// Constant iterator over Vars
        typedef std::vector<Var>::const_iterator const_iterator;
        /// Iterator over Vars
        typedef std::vector<Var>::iterator iterator;
        /// Constant reverse iterator over Vars
        typedef std::vector<Var>::const_reverse_iterator const_reverse_iterator;
        /// Reverse iterator over Vars
        typedef std::vector<Var>::reverse_iterator reverse_iterator;
        
        /// Returns iterator that points to the first variable
        iterator begin() { return _vars.begin(); }
        /// Returns constant iterator that points to the first variable
        const_iterator begin() const { return _vars.begin(); }

        /// Returns iterator that points beyond the last variable
        iterator end() { return _vars.end(); }
        /// Returns constant iterator that points beyond the last variable
        const_iterator end() const { return _vars.end(); }

        /// Returns reverse iterator that points to the last variable
        reverse_iterator rbegin() { return _vars.rbegin(); }
        /// Returns constant reverse iterator that points to the last variable
        const_reverse_iterator rbegin() const { return _vars.rbegin(); }

        /// Returns reverse iterator that points beyond the first variable
        reverse_iterator rend() { return _vars.rend(); }
        /// Returns constant reverse iterator that points beyond the first variable
        const_reverse_iterator rend() const { return _vars.rend(); }


        /// Returns number of variables
        std::vector<Var>::size_type size() const { return _vars.size(); }


        /// Returns whether the VarSet is empty
        bool empty() const { return _vars.size() == 0; }


        /// Test for equality of variable labels
        friend bool operator==( const VarSet &a, const VarSet &b ) {
            return (a._vars == b._vars);
        }

        /// Test for inequality of variable labels
        friend bool operator!=( const VarSet &a, const VarSet &b ) {
            return !(a._vars == b._vars);
        }

        /// Lexicographical comparison of variable labels
        friend bool operator<( const VarSet &a, const VarSet &b ) {
            return a._vars < b._vars;
        }

        /// calcState calculates the linear index of this VarSet that corresponds
        /// to the states of the variables given in states, implicitly assuming
        /// states[m] = 0 for all m in this VarSet which are not in states.
        size_t calcState( const std::map<Var, size_t> &states ) const {
            size_t prod = 1;
            size_t state = 0;
            foreach( const Var &n, *this ) {
                std::map<Var, size_t>::const_iterator m = states.find( n );
                if( m != states.end() )
                    state += prod * m->second;
                prod *= n.states();
            }
            return state;
        }

    private:
        /// Calculates the number of states
        size_t calcStates() {
            _states = 1;
            foreach( Var &i, _vars )
                _states *= i.states();
            return _states;
        }
};


/// For two Vars n1 and n2, the expression n1 | n2 gives the Varset containing n1 and n2
inline VarSet operator| (const Var& n1, const Var& n2) {
    return( VarSet(n1, n2) );
}


} // end of namespace dai


#endif
