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
/// \brief Defines the IndexFor, MultiFor, Permute and State classes


#ifndef __defined_libdai_index_h
#define __defined_libdai_index_h


#include <vector>
#include <algorithm>
#include <map>
#include <cassert>
#include <dai/varset.h>


namespace dai {


/// Tool for looping over the states of several variables.
/** The class IndexFor is an important tool for indexing Factor entries.
 *  Its usage can best be explained by an example.
 *  Assume indexVars, forVars are both VarSets.
 *  Then the following code:
 *  \code
 *      IndexFor i( indexVars, forVars );
 *      for( ; i >= 0; ++i ) {
 *          // use long(i)
 *      }
 *  \endcode
 *  loops over all joint states of the variables in forVars,
 *  and (long)i is equal to the linear index of the corresponding
 *  state of indexVars, where the variables in indexVars that are
 *  not in forVars assume their zero'th value.
 */
class IndexFor {
    private:
        /// The current linear index corresponding to the state of indexVars
        long                _index;

        /// For each variable in forVars, the amount of change in _index
        std::vector<long>   _sum;

        /// For each variable in forVars, the current state
        std::vector<size_t> _count;
        
        /// For each variable in forVars, its number of possible values
        std::vector<size_t> _dims;

    public:
        /// Default constructor
        IndexFor() { 
            _index = -1; 
        }

        /// Constructor
        IndexFor( const VarSet& indexVars, const VarSet& forVars ) : _count( forVars.size(), 0 ) {
            long sum = 1;

            _dims.reserve( forVars.size() );
            _sum.reserve( forVars.size() );

            VarSet::const_iterator j = forVars.begin();
            for( VarSet::const_iterator i = indexVars.begin(); i != indexVars.end(); ++i ) {
                for( ; j != forVars.end() && *j <= *i; ++j ) {
                    _dims.push_back( j->states() );
                    _sum.push_back( (*i == *j) ? sum : 0 );
                }
                sum *= i->states();
            }
            for( ; j != forVars.end(); ++j ) {
                _dims.push_back( j->states() );
                _sum.push_back( 0 );
            }
            _index = 0;
        }

        /// Copy constructor
        IndexFor( const IndexFor & ind ) : _index(ind._index), _sum(ind._sum), _count(ind._count), _dims(ind._dims) {}

        /// Assignment operator
        IndexFor& operator=( const IndexFor &ind ) {
            if( this != &ind ) {
                _index = ind._index;
                _sum = ind._sum;
                _count = ind._count;
                _dims = ind._dims;
            }
            return *this;
        }

        /// Sets the index back to zero
        IndexFor& clear() {
            fill( _count.begin(), _count.end(), 0 );
            _index = 0;
            return( *this );
        }

        /// Conversion to long
        operator long () const { 
            return( _index ); 
        }

        /// Pre-increment operator
        IndexFor& operator++ () {
            if( _index >= 0 ) {
                size_t i = 0;

                while( i < _count.size() ) {
                    _index += _sum[i];
                    if( ++_count[i] < _dims[i] )
                        break;
                    _index -= _sum[i] * _dims[i];
                    _count[i] = 0;
                    i++;
                }

                if( i == _count.size() ) 
                    _index = -1;
            }
            return( *this );
        }
};


/// MultiFor makes it easy to perform a dynamic number of nested for loops.
/** An example of the usage is as follows:
 *  \code
 *  std::vector<size_t> dims;
 *  dims.push_back( 3 );
 *  dims.push_back( 4 );
 *  dims.push_back( 5 );
 *  for( MultiFor s(dims); s.valid(); ++s )
 *      cout << "linear index: " << (size_t)s << " corresponds to indices " << s[0] << ", " << s[1] << ", " << s[2] << endl;
 *  \endcode
 *  which would be equivalent to:
 *  \code
 *  size_t s = 0;
 *  for( size_t s0 = 0; s0 < 3; s0++ )
 *      for( size_t s1 = 0; s1 < 4; s1++ )
 *          for( size_t s2 = 0; s2 < 5; s++, s2++ )
 *              cout << "linear index: " << (size_t)s << " corresponds to indices " << s0 << ", " << s1 << ", " << s2 << endl;
 *  \endcode
 */
class MultiFor {
    private:
        std::vector<size_t>  _dims;
        std::vector<size_t>  _states;
        long                 _state;

    public:
        /// Default constructor
        MultiFor() : _dims(), _states(), _state(0) {}

        /// Initialize from vector of index dimensions
        MultiFor( const std::vector<size_t> &d ) : _dims(d), _states(d.size(),0), _state(0) {}

        /// Copy constructor
        MultiFor( const MultiFor &x ) : _dims(x._dims), _states(x._states), _state(x._state) {}

        /// Assignment operator
        MultiFor& operator=( const MultiFor & x ) {
            if( this != &x ) {
                _dims   = x._dims;
                _states = x._states;
                _state  = x._state;
            }
            return *this;
        }

        /// Return linear state
        operator size_t() const { 
            assert( valid() );
            return( _state );
        }

        /// Return k'th index
        size_t operator[]( size_t k ) const {
            assert( valid() );
            assert( k < _states.size() );
            return _states[k];
        }

        /// Prefix increment operator
        MultiFor & operator++() {
            if( valid() ) {
                _state++;
                size_t i;
                for( i = 0; i != _states.size(); i++ ) {
                    if( ++(_states[i]) < _dims[i] )
                        break;
                    _states[i] = 0;
                }
                if( i == _states.size() )
                    _state = -1;
            }
            return *this;
        }

        /// Postfix increment operator
        void operator++( int ) {
            operator++();
        }

        /// Returns true if the current state is valid
        bool valid() const {
            return( _state >= 0 );
        }
};


/// Tool for calculating permutations of multiple indices.
class Permute {
    private:
        std::vector<size_t>  _dims;
        std::vector<size_t>  _sigma;

    public:
        /// Default constructor
        Permute() : _dims(), _sigma() {}

        /// Initialize from vector of index dimensions and permutation sigma
        Permute( const std::vector<size_t> &d, const std::vector<size_t> &sigma ) : _dims(d), _sigma(sigma) {
            assert( _dims.size() == _sigma.size() );
        }

        /// Copy constructor
        Permute( const Permute &x ) : _dims(x._dims), _sigma(x._sigma) {}

        /// Assignment operator
        Permute& operator=( const Permute &x ) {
            if( this != &x ) {
                _dims  = x._dims;
                _sigma = x._sigma;
            }
            return *this;
        }

        /// Converts the linear index li to a vector index
        /// corresponding with the dimensions in _dims,
        /// permutes it according to sigma, 
        /// and converts it back to a linear index
        /// according to the permuted dimensions.
        size_t convert_linear_index( size_t li ) {
            size_t N = _dims.size();

            // calculate vector index corresponding to linear index
            std::vector<size_t> vi;
            vi.reserve( N );
            size_t prod = 1;
            for( size_t k = 0; k < N; k++ ) {
                vi.push_back( li % _dims[k] );
                li /= _dims[k];
                prod *= _dims[k];
            }

            // convert permuted vector index to corresponding linear index
            prod = 1;
            size_t sigma_li = 0;
            for( size_t k = 0; k < N; k++ ) {
                sigma_li += vi[_sigma[k]] * prod;
                prod *= _dims[_sigma[k]];
            }

            return sigma_li;
        }
};


/// Contains the joint state of variables within a VarSet and useful things to do with this information.
/** This is very similar to a MultiFor, but tailored for Vars and Varsets.
 */
class State {
    private:
        typedef std::map<Var, size_t> states_type;

        long                          state;
        states_type                   states;
        
    public:
        /// Default constructor
        State() : state(0), states() {}

        /// Initialize from VarSet
        State( const VarSet &vs ) : state(0) {
            for( VarSet::const_iterator v = vs.begin(); v != vs.end(); v++ )
                states[*v] = 0;
        }

        /// Copy constructor
        State( const State & x ) : state(x.state), states(x.states) {}

        /// Assignment operator
        State& operator=( const State &x ) {
            if( this != &x ) {
                state  = x.state;
                states = x.states;
            }
            return *this;
        }

        /// Return linear state
        operator size_t() const { 
            assert( valid() );
            return( state );
        }

        /// Return state of variable n,
        /// or zero if n is not in this State
        size_t operator() ( const Var &n ) const {
            assert( valid() );
            states_type::const_iterator entry = states.find( n );
            if( entry == states.end() )
                return 0;
            else
                return entry->second;
        }

        /// Return linear state of variables in varset,
        /// setting them to zero if they are not in this State
        size_t operator() ( const VarSet &vs ) const {
            assert( valid() );
            size_t vs_state = 0;
            size_t prod = 1;
            for( VarSet::const_iterator v = vs.begin(); v != vs.end(); v++ ) {
                states_type::const_iterator entry = states.find( *v );
                if( entry != states.end() )
                    vs_state += entry->second * prod; 
                prod *= v->states();
            }
            return vs_state;
        }

        /// Postfix increment operator
        void operator++( int ) {
            if( valid() ) {
                state++;
                states_type::iterator entry = states.begin();
                while( entry != states.end() ) {
                    if( ++(entry->second) < entry->first.states() )
                        break;
                    entry->second = 0;
                    entry++;
                }
                if( entry == states.end() )
                    state = -1;
            }
        }

        /// Returns true if the current state is valid
        bool valid() const {
            return( state >= 0 );
        }
};


} // end of namespace dai


#endif
