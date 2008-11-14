/*  Copyright (C) 2008  Frederik Eaton [frederik at ofb dot net]

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


#ifndef __defined_libdai_gibbs_h
#define __defined_libdai_gibbs_h


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


class Gibbs : public DAIAlgFG {
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

    protected:
        typedef std::vector<size_t> _count_t;
        size_t _sample_count;
        std::vector<_count_t> _var_counts;
        std::vector<_count_t> _factor_counts;

        typedef std::vector<size_t> _state_t;
        void update_counts(_state_t &st);
        void randomize_state(_state_t &st);
        Prob get_var_dist(_state_t &st, size_t i);
        void resample_var(_state_t &st, size_t i);
        size_t get_factor_entry(const _state_t &st, int factor);

    public:
        // default constructor
        Gibbs() : DAIAlgFG() {}
        // copy constructor
        Gibbs(const Gibbs & x) : DAIAlgFG(x), _sample_count(x._sample_count), _var_counts(x._var_counts), _factor_counts(x._factor_counts) {}
        // construct Gibbs object from FactorGraph
        Gibbs( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg) {
            setProperties( opts );
            construct();
        }
        // assignment operator
        Gibbs & operator=(const Gibbs & x) {
            if(this!=&x) {
                DAIAlgFG::operator=(x);
                _sample_count = x._sample_count;
                _var_counts = x._var_counts;
                _factor_counts = x._factor_counts;
            }
            return *this;
        }
        
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

        Factor beliefV( size_t i ) const;
        Factor beliefF( size_t I ) const;

        void construct();
        /// Set Props according to the PropertySet opts, where the values can be stored as std::strings or as the type of the corresponding Props member
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
