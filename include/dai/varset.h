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
/// \brief Defines VarSet class


#ifndef __defined_libdai_varset_h
#define __defined_libdai_varset_h


#include <vector>
#include <map>
#include <ostream>
#include <dai/var.h>
#include <dai/util.h>
#include <dai/smallset.h>


namespace dai {


/// Represents a set of variables.
/** \note A VarSet is implemented using a std::vector<Var> instead
 *  of the more natural std::set<Var> because of efficiency reasons.
 */
class VarSet : public smallSet<Var> {
    public:
        /// Default constructor
        VarSet() : smallSet<Var>() {}

        /// Copy constructor
        VarSet( const VarSet &x ) : smallSet<Var>(x) {}

        /// Assignment operator
        VarSet& operator=( const VarSet &x ) {
            if( this != &x ) {
                smallSet<Var>::operator=( x );
            }
            return *this;
        }
        
        /// Construct from smallSet<Var>
        VarSet( const smallSet<Var> &x ) : smallSet<Var>(x) {}

        /// Calculates the product of the number of states of all variables in this VarSet.
        size_t nrStates() {
            size_t states = 1;
            for( VarSet::const_iterator n = begin(); n != end(); n++ )
                states *= n->states();
            return states;
        }

        /// Construct a VarSet with one element
        VarSet( const Var &n ) : smallSet<Var>(n) {}

        /// Construct a VarSet with two elements
        VarSet( const Var &n1, const Var &n2 ) : smallSet<Var>(n1,n2) {} 

        /// Construct a VarSet from a range of iterators.
        /** \tparam VarIterator Iterator with value_type Var.
         *  \param begin Points to first Var to be added.
         *  \param end Points just beyond last Var to be added.
         *  \param sizeHint For efficiency, the number of elements can be speficied by sizeHint.
         */
        template <typename VarIterator>
        VarSet( VarIterator begin, VarIterator end, size_t sizeHint=0 ) : smallSet<Var>(begin,end,sizeHint) {}

        /// Calculates the linear index in the cartesian product of the variables in *this, which corresponds to a particular joint assignment of the variables.
        /** \param states Specifies the states of some variables.
         *  \return The linear index in the cartesian product of the variables in *this
         *  corresponding with the joint assignment specified by \c states (where it is
         *  assumed that states[m] == 0 for all m in vars which are not in states).
         */
        size_t calcState( const std::map<Var, size_t> &states ) {
            size_t prod = 1;
            size_t state = 0;
            for( VarSet::const_iterator n = begin(); n != end(); n++ ) {
                std::map<Var, size_t>::const_iterator m = states.find( *n );
                if( m != states.end() )
                    state += prod * m->second;
                prod *= n->states();
            }
            return state;
        }

        /// Writes a VarSet to an output stream
        friend std::ostream& operator<< (std::ostream &os, const VarSet& ns)  {
            for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
                os << *n;
            return( os );
        }
};


} // end of namespace dai


#endif
