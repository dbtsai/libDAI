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
/// \brief Defines SmallSet<T> class


#ifndef __defined_libdai_smallset_h
#define __defined_libdai_smallset_h


#include <vector>
#include <algorithm>


namespace dai {


/// Represents a set; the implementation is optimized for a small number of elements.
/** SmallSet uses an ordered std::vector<T> to represent a set; this is faster than
 *  using a std::set<T> if the number of elements is small.
 *  \tparam T Should be less-than-comparable.
 */
template <typename T>
class SmallSet {
    private:
        /// The elements in this set
        std::vector<T> _elements;

    public:
        /// Default constructor (construct an empty set)
        SmallSet() : _elements() {}

        /// Construct a set with one element
        SmallSet( const T &t ) : _elements() {
            _elements.push_back( t );
        }

        /// Construct a set with two elements
        SmallSet( const T &t1, const T &t2 ) {
            if( t1 < t2 ) {
                _elements.push_back( t1 );
                _elements.push_back( t2 );
            } else if( t2 < t1 ) {
                _elements.push_back( t2 );
                _elements.push_back( t1 );
            } else
                _elements.push_back( t1 );
        }

        /// Construct a SmallSet from a range of elements.
        /** \tparam TIterator Iterates over instances of type T.
         *  \param begin Points to first element to be added.
         *  \param end Points just beyond last element to be added.
         *  \param sizeHint For efficiency, the number of elements can be speficied by sizeHint.
         */
        template <typename TIterator>
        SmallSet( TIterator begin, TIterator end, size_t sizeHint=0 ) {
            _elements.reserve( sizeHint );
            _elements.insert( _elements.begin(), begin, end );
            std::sort( _elements.begin(), _elements.end() );
            typename std::vector<T>::iterator new_end = std::unique( _elements.begin(), _elements.end() );
            _elements.erase( new_end, _elements.end() );
        }

        /// Set-minus operator: returns all elements in *this, except those in x
        SmallSet operator/ ( const SmallSet& x ) const {
            SmallSet res;
            std::set_difference( _elements.begin(), _elements.end(), x._elements.begin(), x._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }

        /// Set-union operator: returns all elements in *this, plus those in x
        SmallSet operator| ( const SmallSet& x ) const {
            SmallSet res;
            std::set_union( _elements.begin(), _elements.end(), x._elements.begin(), x._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }

        /// Set-intersection operator: returns all elements in *this that are also contained in x
        SmallSet operator& ( const SmallSet& x ) const {
            SmallSet res;
            std::set_intersection( _elements.begin(), _elements.end(), x._elements.begin(), x._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }

        /// Erases from *this all elements in x
        SmallSet& operator/= ( const SmallSet& x ) {
            return (*this = (*this / x));
        }

        /// Erases one element
        SmallSet& operator/= ( const T &t ) {
            typename std::vector<T>::iterator pos = lower_bound( _elements.begin(), _elements.end(), t );
            if( pos != _elements.end() )
                if( *pos == t ) // found element, delete it
                    _elements.erase( pos );
            return *this;
        }

        /// Adds to *this all elements in x
        SmallSet& operator|= ( const SmallSet& x ) {
            return( *this = (*this | x) );
        }

        /// Adds one element
        SmallSet& operator|= ( const T& t ) {
            typename std::vector<T>::iterator pos = lower_bound( _elements.begin(), _elements.end(), t );
            if( pos == _elements.end() || *pos != t ) // insert it
                _elements.insert( pos, t );
            return *this;
        }

        /// Erases from *this all elements not in x
        SmallSet& operator&= ( const SmallSet& x ) {
            return (*this = (*this & x));
        }

        /// Returns true if *this is a subset of x
        bool operator<< ( const SmallSet& x ) const {
            return std::includes( x._elements.begin(), x._elements.end(), _elements.begin(), _elements.end() );
        }

        /// Returns true if x is a subset of *this
        bool operator>> ( const SmallSet& x ) const {
            return std::includes( _elements.begin(), _elements.end(), x._elements.begin(), x._elements.end() );
        }

        /// Returns true if *this and x have elements in common
        bool intersects( const SmallSet& x ) const {
            return( (*this & x).size() > 0 );
        }

        /// Returns true if *this contains the element t
        bool contains( const T &t ) const {
            return std::binary_search( _elements.begin(), _elements.end(), t );
        }

        /// Constant iterator over the elements
        typedef typename std::vector<T>::const_iterator const_iterator;
        /// Iterator over the elements
        typedef typename std::vector<T>::iterator iterator;
        /// Constant reverse iterator over the elements
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
        /// Reverse iterator over the elements
        typedef typename std::vector<T>::reverse_iterator reverse_iterator;

        /// Returns iterator that points to the first element
        iterator begin() { return _elements.begin(); }
        /// Returns constant iterator that points to the first element
        const_iterator begin() const { return _elements.begin(); }

        /// Returns iterator that points beyond the last element
        iterator end() { return _elements.end(); }
        /// Returns constant iterator that points beyond the last element
        const_iterator end() const { return _elements.end(); }

        /// Returns reverse iterator that points to the last element
        reverse_iterator rbegin() { return _elements.rbegin(); }
        /// Returns constant reverse iterator that points to the last element
        const_reverse_iterator rbegin() const { return _elements.rbegin(); }

        /// Returns reverse iterator that points beyond the first element
        reverse_iterator rend() { return _elements.rend(); }
        /// Returns constant reverse iterator that points beyond the first element
        const_reverse_iterator rend() const { return _elements.rend(); }

        /// Returns number of elements
        typename std::vector<T>::size_type size() const { return _elements.size(); }

        /// Returns whether the SmallSet is empty
        bool empty() const { return _elements.size() == 0; }

        /// Returns true if a and b are identical
        friend bool operator==( const SmallSet &a, const SmallSet &b ) {
            return (a._elements == b._elements);
        }

        /// Returns true if a and b are not identical
        friend bool operator!=( const SmallSet &a, const SmallSet &b ) {
            return !(a._elements == b._elements);
        }

        /// Lexicographical comparison of elements
        friend bool operator<( const SmallSet &a, const SmallSet &b ) {
            return a._elements < b._elements;
        }
};


} // end of namespace dai


#endif
