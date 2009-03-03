/*  Copyright (C) 2008  Frederik Eaton  [frederik at ofb dot net],
                        Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]

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
/// \brief Defines class Gibbs
/// \todo Improve documentation


#ifndef __defined_libdai_gibbs_h
#define __defined_libdai_gibbs_h


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


/// Approximate inference algorithm "Gibbs sampling"
class Gibbs : public DAIAlgFG {
    private:
        typedef std::vector<size_t> _count_t;
        typedef std::vector<size_t> _state_t;

        size_t _sample_count;
        std::vector<_count_t> _var_counts;
        std::vector<_count_t> _factor_counts;
        _state_t _state;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Number of iterations
            size_t iters;

            /// Verbosity
            size_t verbose;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        Gibbs() : DAIAlgFG(), _sample_count(0), _var_counts(), _factor_counts(), _state() {}

        /// Construct from FactorGraph fg and PropertySet opts
        Gibbs( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), _sample_count(0), _var_counts(), _factor_counts(), _state() {
            setProperties( opts );
            construct();
        }


        /// @name General InfAlg interface
        //@{
        virtual Gibbs* clone() const { return new Gibbs(*this); }
        virtual Gibbs* create() const { return new Gibbs(); }
        virtual std::string identify() const { return std::string(Name) + printProperties(); }
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual double run();
        virtual double maxDiff() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual size_t Iterations() const { return props.iters; }
        //@}


        /// @name Additional interface specific for Gibbs
        //@{
        Factor beliefV( size_t i ) const;
        Factor beliefF( size_t I ) const;
        void randomizeState();
        //@}

        /// Return reference to current state vector
        std::vector<size_t>& state() { return _state; }

        /// Return const reference to current state vector
        const std::vector<size_t>& state() const { return _state; }

    private:
        void updateCounts();
        Prob getVarDist( size_t i );
        void resampleVar( size_t i );
        size_t getFactorEntry( size_t I );
        size_t getFactorEntryDiff( size_t I, size_t i );

        void construct();
        /// Set Props according to the PropertySet opts, where the values can be stored as std::strings or as the type of the corresponding Props member
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
