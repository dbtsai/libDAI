/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2008  Frederik Eaton  [frederik at ofb dot net]
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


    /// \name General InfAlg interface
    //@{
        virtual Gibbs* clone() const { return new Gibbs(*this); }
        virtual std::string identify() const { return std::string(Name) + printProperties(); }
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual Factor beliefF( size_t I ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual Real run();
        virtual Real maxDiff() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual size_t Iterations() const { return props.iters; }
    //@}


    /// \name Additional interface specific for Gibbs
    //@{
        void randomizeState();
        /// Return reference to current state vector
        std::vector<size_t>& state() { return _state; }

        /// Return const reference to current state vector
        const std::vector<size_t>& state() const { return _state; }
    //@}

    /// \name Managing parameters (which are stored in Gibbs::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding Gibbs::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        void updateCounts();
        Prob getVarDist( size_t i );
        void resampleVar( size_t i );
        size_t getFactorEntry( size_t I );
        size_t getFactorEntryDiff( size_t I, size_t i );

        void construct();
};


} // end of namespace dai


#endif
