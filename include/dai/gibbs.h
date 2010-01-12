/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2008  Frederik Eaton  [frederik at ofb dot net]
 */


/// \file
/// \brief Defines class Gibbs, which implements Gibbs sampling


#ifndef __defined_libdai_gibbs_h
#define __defined_libdai_gibbs_h


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


/// Approximate inference algorithm "Gibbs sampling"
/** \author Frederik Eaton
 */
class Gibbs : public DAIAlgFG {
    private:
        /// Type used to store the counts of various states
        typedef std::vector<size_t> _count_t;
        /// Type used to store the joint state of all variables
        typedef std::vector<size_t> _state_t;
        /// Number of samples counted so far (excluding burn-in)
        size_t _sample_count;
        /// State counts for each variable
        std::vector<_count_t> _var_counts;
        /// State counts for each factor
        std::vector<_count_t> _factor_counts;
        /// Current joint state of all variables
        _state_t _state;

    public:
        /// Parameters for Gibbs
        struct Properties {
            /// Total number of iterations
            size_t iters;

            /// Number of "burn-in" iterations
            size_t burnin;

            /// Verbosity (amount of output sent to stderr)
            size_t verbose;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        Gibbs() : DAIAlgFG(), _sample_count(0), _var_counts(), _factor_counts(), _state() {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param opts Parameters @see Properties
         */
        Gibbs( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), _sample_count(0), _var_counts(), _factor_counts(), _state() {
            setProperties( opts );
            construct();
        }


    /// \name General InfAlg interface
    //@{
        virtual Gibbs* clone() const { return new Gibbs(*this); }
        virtual std::string identify() const { return std::string(Name) + printProperties(); }
        virtual Factor belief( const Var &v ) const { return beliefV( findVar( v ) ); }
        virtual Factor belief( const VarSet &vs ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual Factor beliefF( size_t I ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual Real run();
        virtual Real maxDiff() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual size_t Iterations() const { return props.iters; }
        virtual void setProperties( const PropertySet &opts );
        virtual PropertySet getProperties() const;
        virtual std::string printProperties() const;
    //@}


    /// \name Additional interface specific for Gibbs
    //@{
        /// Draw the current joint state of all variables from a uniform random distribution
        void randomizeState();
        /// Return reference to current state of all variables
        std::vector<size_t>& state() { return _state; }
        /// Return constant reference to current state of all variables
        const std::vector<size_t>& state() const { return _state; }
    //@}

    private:
        /// Helper function for constructors
        void construct();
        /// Updates all counts (_sample_count, _var_counts, _factor_counts) based on current state
        void updateCounts();
        /// Calculate conditional distribution of variable \a i, given the current state
        Prob getVarDist( size_t i );
        /// Draw state of variable \a i randomly from its conditional distribution and update the current state
        void resampleVar( size_t i );
        /// Calculates linear index into factor \a I corresponding to the current state
        size_t getFactorEntry( size_t I );
        /// Calculates the differences between linear indices into factor \a I corresponding with a state change of variable \a i
        size_t getFactorEntryDiff( size_t I, size_t i );
};


/// Runs Gibbs sampling for \a iters iterations (of which \a burnin for burn-in) on FactorGraph \a fg, and returns the resulting state
/** \relates Gibbs
 */
std::vector<size_t> getGibbsState( const FactorGraph &fg, size_t iters );


} // end of namespace dai


/** \example example_sprinkler_gibbs.cpp
 *  This example shows how to use the Gibbs class.
 */


#endif
