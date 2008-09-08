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


#ifndef __defined_libdai_varset_h
#define __defined_libdai_varset_h


#include <set>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <dai/var.h>


namespace dai {


/// VarSet represents a set of variables and is a descendant of set<Var>. 
/// In addition, it provides an easy interface for set-theoretic operations
/// by operator overloading.
class VarSet : private std::set<Var> {
    protected:
        /// Product of number of states of all contained variables
        size_t _statespace;

        /// Check whether ns is a subset
        bool includes( const VarSet& ns ) const {
            return std::includes( begin(), end(), ns.begin(), ns.end() );
        }

        /// Calculate statespace
        size_t calcStateSpace() {
            _statespace = 1;
            for( const_iterator i = begin(); i != end(); ++i )
                _statespace *= i->states();
            return _statespace;
        }


    public:
        /// Default constructor
        VarSet() : _statespace(0) {};

        /// Construct a VarSet with one variable
        VarSet( const Var &n ) : _statespace( n.states() ) { 
            insert( n ); 
        }

        /// Construct a VarSet with two variables
        VarSet( const Var &n1, const Var &n2 ) { 
            insert( n1 ); 
            insert( n2 ); 
            calcStateSpace();
        }

        /// Construct from a set<Var>
        VarSet( const std::set<Var> &ns ) {
            std::set<Var>::operator=( ns );
            calcStateSpace();
        }

        /// Copy constructor
        VarSet( const VarSet &x ) : std::set<Var>( x ), _statespace( x._statespace ) {}

        /// Assignment operator
        VarSet & operator=( const VarSet &x ) {
            if( this != &x ) {
                std::set<Var>::operator=( x );
                _statespace = x._statespace;
            }
            return *this;
        }
        

        /// Return statespace, i.e. the product of the number of states of each variable
        size_t states() const { 
#ifdef DAI_DEBUG
            size_t x = 1;
            for( const_iterator i = begin(); i != end(); ++i )
                x *= i->states();
            assert( x == _statespace );
#endif
            return _statespace; 
        }
        

        /// Erase one variable
        VarSet& operator/= (const Var& n) { 
            erase( n ); 
            calcStateSpace();
            return *this; 
        }

        /// Add one variable
        VarSet& operator|= (const Var& n) {
            insert( n ); 
            calcStateSpace();
            return *this;
        }

        /// Setminus operator (result contains all variables except those in ns)
        VarSet operator/ (const VarSet& ns) const {
            VarSet res;
            std::set_difference( begin(), end(), ns.begin(), ns.end(), inserter( res, res.begin() ) );
            res.calcStateSpace();
            return res;
        }

        /// Set-union operator (result contains all variables plus those in ns)
        VarSet operator| (const VarSet& ns) const {
            VarSet res;
            std::set_union( begin(), end(), ns.begin(), ns.end(), inserter( res, res.begin() ) );
            res.calcStateSpace();
            return res;
        }

        /// Set-intersection operator (result contains all variables that are also contained in ns)
        VarSet operator& (const VarSet& ns) const {
            VarSet res;
            std::set_intersection( begin(), end(), ns.begin(), ns.end(), inserter( res, res.begin() ) );
            res.calcStateSpace();
            return res;
        }
        
        /// Erases from *this all variables in ns
        VarSet& operator/= (const VarSet& ns) {
            return (*this = (*this / ns));
        }

        /// Adds to *this all variables in ns
        VarSet& operator|= (const VarSet& ns) {
            return (*this = (*this | ns));
        }

        /// Erases from *this all variables not in ns
        VarSet& operator&= (const VarSet& ns) { 
            return (*this = (*this & ns)); 
        }
        

        /// Returns false if both *this and ns are empty
        bool operator|| (const VarSet& ns) const { 
            return !( this->empty() && ns.empty() );
        }

        /// Returns true if *this and ns contain common variables
        bool operator&& (const VarSet& ns) const { 
            return !( (*this & ns).empty() ); 
        }

        /// Returns true if *this is a subset of ns
        bool operator<< (const VarSet& ns) const { 
            return ns.includes( *this ); 
        }

        /// Returns true if ns is a subset of *this
        bool operator>> (const VarSet& ns) const { 
            return includes( ns ); 
        }

        /// Returns true if *this contains the variable n
        bool operator&& (const Var& n) const { 
            return( find( n ) == end() ? false : true ); 
        }

        
        /// Sends a VarSet to an output stream
        friend std::ostream& operator<< (std::ostream & os, const VarSet& ns) {
            for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++)
                os << *n;
            return( os );
        }

        
/*      The following makes part of the public interface of set<Var> available.
 *      It is important to note that insert functions have to be overloaded,
 *      because they have to recalculate the statespace. A different approach
 *      would be to publicly inherit from set<Var> and only overload the insert
 *      methods.
 */
        
        // Additional interface from set<Var> that has to be provided
        using std::set<Var>::const_iterator;
        using std::set<Var>::iterator;
        using std::set<Var>::const_reference;
        using std::set<Var>::begin;
        using std::set<Var>::end;
        using std::set<Var>::size;
        using std::set<Var>::empty;

        /// Copy of set<Var>::insert which additionally calculates the new statespace
        std::pair<iterator, bool> insert( const Var& x ) {
            std::pair<iterator, bool> result = std::set<Var>::insert( x );
            calcStateSpace();
            return result;
        }

        /// Copy of set<Var>::insert which additionally calculates the new statespace
        iterator insert( iterator pos, const value_type& x ) {
            iterator result = std::set<Var>::insert( pos, x );
            calcStateSpace();
            return result;
        }

        /// Test for equality (ignore _statespace member)
        friend bool operator==( const VarSet &a, const VarSet &b ) {
            return operator==( (std::set<Var>)a, (std::set<Var>)b );
        }

        /// Test for inequality (ignore _statespace member)
        friend bool operator!=( const VarSet &a, const VarSet &b ) {
            return operator!=( (std::set<Var>)a, (std::set<Var>)b );
        }

        friend bool operator<( const VarSet &a, const VarSet &b ) {
            return operator<( (std::set<Var>)a, (std::set<Var>)b );
        }
};


/// For two Vars n1 and n2, the expression n1 | n2 gives the Varset containing n1 and n2
inline VarSet operator| (const Var& n1, const Var& n2) {
    return( VarSet(n1, n2) );
}


} // end of namespace dai


#endif
