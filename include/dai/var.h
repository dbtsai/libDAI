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
/// \brief Defines class Var


#ifndef __defined_libdai_var_h
#define __defined_libdai_var_h


#include <iostream>


namespace dai {


/// Represents a discrete random variable.
/** It contains the \a label of the variable (an integer-valued unique ID) 
 *  and the number of possible values (\a states) of the variable.
 */
class Var {
    private:
        /// Label of the variable (its unique ID)
        long    _label;

        /// Number of possible values
        size_t  _states;
        
    public:
        /// Default constructor
        Var() : _label(-1), _states(0) {}
        /// Constructor
        Var( long label, size_t states ) : _label(label), _states(states) {}

        /// Returns the label
        long label() const { return _label; }
        /// Returns reference to label
        long & label() { return _label; }

        /// Returns the number of states
        size_t states () const { return _states; }
        /// Returns reference to number of states
        size_t& states () { return _states; }

        /// Smaller-than operator (only compares labels)
        bool operator < ( const Var& n ) const { return( _label <  n._label ); }
        /// Larger-than operator (only compares labels)
        bool operator > ( const Var& n ) const { return( _label >  n._label ); }
        /// Smaller-than-or-equal-to operator (only compares labels)
        bool operator <= ( const Var& n ) const { return( _label <= n._label ); }
        /// Larger-than-or-equal-to operator (only compares labels)
        bool operator >= ( const Var& n ) const { return( _label >= n._label ); }
        /// Not-equal-to operator (only compares labels)
        bool operator != ( const Var& n ) const { return( _label != n._label ); }
        /// Equal-to operator (only compares labels)
        bool operator == ( const Var& n ) const { return( _label == n._label ); }

        /// Writes a Var to an output stream
        friend std::ostream& operator << ( std::ostream& os, const Var& n ) {
            return( os << "[" << n.label() << "]" );
        }
};


} // end of namespace dai


#endif
