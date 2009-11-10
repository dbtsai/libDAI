/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2009  Frederik Eaton [frederik at ofb dot net]
 */


/// \file
/// \brief Defines class CBP [\ref EaG09]
/// \author Frederik Eaton
/// \todo Improve documentation


#ifndef __defined_libdai_cbp_h
#define __defined_libdai_cbp_h


#include <fstream>
#include <boost/shared_ptr.hpp>

#include <dai/daialg.h>
#include <dai/bbp.h>


namespace dai {


/// Find a variable to clamp using BBP (goes with maximum adjoint)
/// \see BBP
std::pair<size_t, size_t> bbpFindClampVar( const InfAlg &in_bp, bool clampingVar, const PropertySet &bbp_props, const BBPCostFunction &cfn, Real *maxVarOut );


/// Class for CBP (Clamped Belief Propagation)
/** This algorithm uses configurable heuristics to choose a variable
 *  x_i and a state x_i*. Inference is done with x_i "clamped" to x_i*
 *  (i.e., conditional on x_i == x_i*), and also with the negation of this
 *  condition. Clamping is done recursively up to a fixed number of
 *  levels (other stopping criteria are also implemented, see
 *  \a recursion property). The resulting approximate marginals are
 *  combined using logZ estimates.
 *
 *  \author Frederik Eaton
 */
class CBP : public DAIAlgFG {
    private:
        /// Variable beliefs
        std::vector<Factor> _beliefsV;
        /// Factor beliefs
        std::vector<Factor> _beliefsF;
        /// Log-partition sum
        Real _logZ;

        /// Counts number of clampings at each leaf node
        Real _sum_level;

        /// Number of leaves of recursion tree
        size_t _num_leaves;

        /// Output stream where information about the clampings is written
        boost::shared_ptr<std::ofstream> _clamp_ofstream;

        /// Returns BBP cost function used
        BBPCostFunction BBP_cost_function() { return props.bbp_cfn; }

        /// Prints beliefs, variables and partition sum, in case of a debugging build
        void printDebugInfo();

        /// Called by 'run', and by itself. Implements the main algorithm.
        /** Chooses a variable to clamp, recurses, combines the logZ and
         *  beliefs estimates of the children, and returns the improved
         *  estimates in \a lz_out and \a beliefs_out to its parent
         */
        void runRecurse( InfAlg *bp, Real orig_logZ, std::vector<size_t> clamped_vars_list, size_t &num_leaves,
                         size_t &choose_count, Real &sum_level, Real &lz_out, std::vector<Factor> &beliefs_out );

        /// Choose the next variable to clamp
        /** Choose the next variable to clamp, given a converged InfAlg (\a bp),
         *  and a vector of variables that are already clamped (\a
         *  clamped_vars_list). Returns the chosen variable in \a i, and
         *  the set of states in \a xis. If \a maxVarOut is non-NULL and
         *  props.choose==CHOOSE_BBP then it is used to store the
         *  adjoint of the chosen variable
         */
        virtual bool chooseNextClampVar( InfAlg* bp, std::vector<size_t> &clamped_vars_list, size_t &i, std::vector<size_t> &xis, Real *maxVarOut );

        /// Return the InfAlg to use at each step of the recursion.
        /// \todo At present, only returns a BP instance
        InfAlg* getInfAlg();

        /// Numer of iterations needed
        size_t _iters;
        /// Maximum difference encountered so far
        Real _maxdiff;

        /// Sets variable beliefs, factor beliefs and logZ
        /** \param bs should be a concatenation of the variable beliefs followed by the factor beliefs
         */
        void setBeliefs( const std::vector<Factor> &bs, Real logZ );

        /// Constructor helper function
        void construct();

    public:
        /// Construct CBP object from FactorGraph fg and PropertySet opts
        CBP( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg) {
            props.set( opts );
            construct();
        }

        /// Name of this inference algorithm
        static const char *Name;

    /// \name General InfAlg interface
    //@{
        virtual CBP* clone() const { return new CBP(*this); }
        virtual std::string identify() const { return std::string(Name) + props.toString(); }
        virtual Factor belief (const Var &n) const { return _beliefsV[findVar(n)]; }
        virtual Factor belief (const VarSet &) const { DAI_THROW(NOT_IMPLEMENTED); }
        virtual Factor beliefV( size_t i ) const { return _beliefsV[i]; }
        virtual Factor beliefF( size_t I ) const { return _beliefsF[I]; }
        virtual std::vector<Factor> beliefs() const { return concat(_beliefsV, _beliefsF); }
        virtual Real logZ() const { return _logZ; }
        virtual void init() {};
        virtual void init( const VarSet & ) {};
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        virtual void setProperties( const PropertySet &opts ) { props.set( opts ); }
        virtual PropertySet getProperties() const { return props.get(); }
        virtual std::string printProperties() const { return props.toString(); }
    //@}

        //----------------------------------------------------------------

        /// Parameters of this inference algorithm
        /* PROPERTIES(props,CBP) {
            /// Enumeration of possible update schedules
            typedef BP::Properties::UpdateType UpdateType;
            /// Enumeration of possible methods for deciding when to stop recursing
            DAI_ENUM(RecurseType,REC_FIXED,REC_LOGZ,REC_BDIFF);
            /// Enumeration of possible heuristics for choosing clamping variable
            DAI_ENUM(ChooseMethodType,CHOOSE_RANDOM,CHOOSE_MAXENT,CHOOSE_BBP,CHOOSE_BP_L1,CHOOSE_BP_CFN);
            /// Enumeration of possible clampings: variables or factors
            DAI_ENUM(ClampType,CLAMP_VAR,CLAMP_FACTOR);

            /// Verbosity
            size_t verbose = 0;

            /// Tolerance to use in BP
            Real tol;
            /// Update style for BP
            UpdateType updates;
            /// Maximum number of iterations for BP
            size_t maxiter;

            /// Tolerance to use for controlling recursion depth (\a recurse is REC_LOGZ or REC_BDIFF)
            Real rec_tol;
            /// Maximum number of levels of recursion (\a recurse is REC_FIXED)
            size_t max_levels = 10;
            /// If choose==CHOOSE_BBP and maximum adjoint is less than this value, don't recurse
            Real min_max_adj;
            /// Heuristic for choosing clamping variable
            ChooseMethodType choose;
            /// Method for deciding when to stop recursing
            RecurseType recursion;
            /// Whether to clamp variables or factors
            ClampType clamp;
            /// Properties to pass to BBP
            PropertySet bbp_props;
            /// Cost function to use for BBP
            BBPCostFunction bbp_cfn;
            /// Random seed
            size_t rand_seed = 0;

            /// If non-empty, write clamping choices to this file
            std::string clamp_outfile = "";
        }
        */
/* {{{ GENERATED CODE: DO NOT EDIT. Created by
    ./scripts/regenerate-properties include/dai/cbp.h src/cbp.cpp
*/
        struct Properties {
            /// Enumeration of possible update schedules
            typedef BP::Properties::UpdateType UpdateType;
            /// Enumeration of possible methods for deciding when to stop recursing
            DAI_ENUM(RecurseType,REC_FIXED,REC_LOGZ,REC_BDIFF);
            /// Enumeration of possible heuristics for choosing clamping variable
            DAI_ENUM(ChooseMethodType,CHOOSE_RANDOM,CHOOSE_MAXENT,CHOOSE_BBP,CHOOSE_BP_L1,CHOOSE_BP_CFN);
            /// Enumeration of possible clampings: variables or factors
            DAI_ENUM(ClampType,CLAMP_VAR,CLAMP_FACTOR);
            /// Verbosity
            size_t verbose;
            /// Tolerance to use in BP
            Real tol;
            /// Update style for BP
            UpdateType updates;
            /// Maximum number of iterations for BP
            size_t maxiter;
            /// Tolerance to use for controlling recursion depth (\a recurse is REC_LOGZ or REC_BDIFF)
            Real rec_tol;
            /// Maximum number of levels of recursion (\a recurse is REC_FIXED)
            size_t max_levels;
            /// If choose==CHOOSE_BBP and maximum adjoint is less than this value, don't recurse
            Real min_max_adj;
            /// Heuristic for choosing clamping variable
            ChooseMethodType choose;
            /// Method for deciding when to stop recursing
            RecurseType recursion;
            /// Whether to clamp variables or factors
            ClampType clamp;
            /// Properties to pass to BBP
            PropertySet bbp_props;
            /// Cost function to use for BBP
            BBPCostFunction bbp_cfn;
            /// Random seed
            size_t rand_seed;
            /// If non-empty, write clamping choices to this file
            std::string clamp_outfile;

            /// Set members from PropertySet
            void set(const PropertySet &opts);
            /// Get members into PropertySet
            PropertySet get() const;
            /// Convert to a string which can be parsed as a PropertySet
            std::string toString() const;
        } props;
/* }}} END OF GENERATED CODE */

        /// Returns heuristic used for clamping variable
        Properties::ChooseMethodType ChooseMethod() { return props.choose; }
        /// Returns method used for deciding when to stop recursing
        Properties::RecurseType Recursion() { return props.recursion; }
        /// Returns clamping type used
        Properties::ClampType Clamping() { return props.clamp; }
        /// Returns maximum number of levels of recursion
        size_t maxClampLevel() { return props.max_levels; }
        /// Returns props.min_max_adj @see CBP::Properties::min_max_adj
        Real minMaxAdj() { return props.min_max_adj; }
        /// Returns tolerance used for controlling recursion depth
        Real recTol() { return props.rec_tol; }
};


/// Given a sorted vector of states \a xis and total state count \a n_states, return a vector of states not in \a xis
std::vector<size_t> complement( std::vector<size_t>& xis, size_t n_states );

/// Computes \f$\frac{\exp(a)}{\exp(a)+\exp(b)}\f$
Real unSoftMax( Real a, Real b );

/// Computes log of sum of exponents, i.e., \f$\log\left(\exp(a) + \exp(b)\right)\f$
Real logSumExp( Real a, Real b );

/// Compute sum of pairwise L-infinity distances of the first \a nv factors in each vector
Real dist( const std::vector<Factor>& b1, const std::vector<Factor>& b2, size_t nv );


} // end of namespace dai


#endif
