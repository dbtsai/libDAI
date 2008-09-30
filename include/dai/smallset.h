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
/// \brief Defines smallSet<T> class


#ifndef __defined_libdai_smallset_h
#define __defined_libdai_smallset_h


#include <vector>
#include <algorithm>


namespace dai {


/// Represents a set (optimized for a small number of elements).
/** For sets consisting of a small number of elements, an implementation using
 *  an ordered std::vector<T> is faster than an implementation using std::set<T>.
 *  The elements should be less-than-comparable.
 */
template <typename T>
class smallSet {
    private:
        /// The elements in this set
        std::vector<T> _elements;

    public:
        /// Default constructor
        smallSet() : _elements() {}

        /// Construct a smallSet with one element
        smallSet( const T &n ) : _elements() { 
            _elements.push_back( n );
        }

        /// Construct a smallSet with two elements
        smallSet( const T &n1, const T &n2 ) { 
            if( n1 < n2 ) {
                _elements.push_back( n1 );
                _elements.push_back( n2 );
            } else if( n2 < n1 ) {
                _elements.push_back( n2 );
                _elements.push_back( n1 );
            } else
                _elements.push_back( n1 );
        }

        /// Construct a smallSet from a range of iterators.
        /** \tparam Iterator Iterator with value_type T.
         *  \param begin Points to first element to be added.
         *  \param end Points just beyond last element to be added.
         *  \param sizeHint For efficiency, the number of elements can be speficied by sizeHint.
         */
        template <typename Iterator>
        smallSet( Iterator begin, Iterator end, size_t sizeHint=0 ) {
            _elements.reserve( sizeHint );
            _elements.insert( _elements.begin(), begin, end );
            std::sort( _elements.begin(), _elements.end() );
            typename std::vector<T>::iterator new_end = std::unique( _elements.begin(), _elements.end() );
            _elements.erase( new_end, _elements.end() );
        }

        /// Copy constructor
        smallSet( const smallSet &x ) : _elements( x._elements ) {}

        /// Assignment operator
        smallSet & operator=( const smallSet &x ) {
            if( this != &x ) {
                _elements = x._elements;
            }
            return *this;
        }
        
        /// Setminus operator: returns all elements in *this, except those in ns
        smallSet operator/ ( const smallSet& ns ) const {
            smallSet res;
            std::set_difference( _elements.begin(), _elements.end(), ns._elements.begin(), ns._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }

        /// Set-union operator: returns all elements in *this, plus those in ns
        smallSet operator| ( const smallSet& ns ) const {
            smallSet res;
            std::set_union( _elements.begin(), _elements.end(), ns._elements.begin(), ns._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }

        /// Set-intersection operator: returns all elements in *this that are also contained in ns
        smallSet operator& ( const smallSet& ns ) const {
            smallSet res;
            std::set_intersection( _elements.begin(), _elements.end(), ns._elements.begin(), ns._elements.end(), inserter( res._elements, res._elements.begin() ) );
            return res;
        }
        
        /// Erases from *this all elements in ns
        smallSet& operator/= ( const smallSet& ns ) {
            return (*this = (*this / ns));
        }

        /// Erases one element
        smallSet& operator/= ( const T& n ) { 
            typename std::vector<T>::iterator pos = lower_bound( _elements.begin(), _elements.end(), n );
            if( pos != _elements.end() )
                if( *pos == n ) // found element, delete it
                    _elements.erase( pos ); 
            return *this; 
        }

        /// Adds to *this all elements in ns
        smallSet& operator|= ( const smallSet& ns ) {
            return( *this = (*this | ns) );
        }

        /// Adds one element
        smallSet& operator|= ( const T& n ) {
            typename std::vector<T>::iterator pos = lower_bound( _elements.begin(), _elements.end(), n );
            if( pos == _elements.end() || *pos != n ) // insert it
                _elements.insert( pos, n );
            return *this;
        }

        /// Erases from *this all elements not in ns
        smallSet& operator&= ( const smallSet& ns ) { 
            return (*this = (*this & ns));
        }

        /// Returns true if *this is a subset of ns
        bool operator<< ( const smallSet& ns ) const { 
            return std::includes( ns._elements.begin(), ns._elements.end(), _elements.begin(), _elements.end() );
        }

        /// Returns true if ns is a subset of *this
        bool operator>> ( const smallSet& ns ) const { 
            return std::includes( _elements.begin(), _elements.end(), ns._elements.begin(), ns._elements.end() );
        }

        /// Returns true if *this and ns contain common elements
        bool intersects( const smallSet& ns ) const { 
            return( (*this & ns).size() > 0 ); 
        }

        /// Returns true if *this contains the element n
        bool contains( const T& n ) const { 
            return std::binary_search( _elements.begin(), _elements.end(), n );
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

        /// Returns whether the smallSet is empty
        bool empty() const { return _elements.size() == 0; }

        /// Returns true if the two sets are identical
        friend bool operator==( const smallSet &a, const smallSet &b ) {
            return (a._elements == b._elements);
        }

        /// Returns true if the two sets are not identical
        friend bool operator!=( const smallSet &a, const smallSet &b ) {
            return !(a._elements == b._elements);
        }

        /// Lexicographical comparison of elements
        friend bool operator<( const smallSet &a, const smallSet &b ) {
            return a._elements < b._elements;
        }
};


} // end of namespace dai


#endif
