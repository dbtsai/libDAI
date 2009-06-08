/*  Copyright (C) 2009  Frederik Eaton [frederik at ofb dot net]

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
/// \brief Defines class CBP [\ref EaG09]
/// \todo Improve documentation
/// \todo Clean up


#ifndef __defined_libdai_cbp_h
#define __defined_libdai_cbp_h


#include <fstream>

#include <boost/shared_ptr.hpp>

#include <dai/daialg.h>

#include <dai/cbp.h>
#include <dai/bbp.h>


namespace dai {


/// Find a variable to clamp using BBP (goes with maximum adjoint) 
/// \see BBP
std::pair<size_t, size_t> bbpFindClampVar(const InfAlg &in_bp, bool clampingVar,
            const PropertySet &bbp_props, bbp_cfn_t cfn, 
            Real *maxVarOut);


/// Class for CBP (Clamped Belief Propagation)
/** This algorithm uses configurable heuristics to choose a variable
 *  x_i and a state x_i*. Inference is done with x_i "clamped" to x_i*
 *  (e.g. conditional on x_i==x_i*), and also with the negation of this
 *  condition. Clamping is done recursively up to a fixed number of
 *  levels (other stopping criteria are also implemented, see
 *  "recursion" property). The resulting approximate marginals are
 *  combined using logZ estimates.
 */
class CBP : public DAIAlgFG {
    std::vector<Factor> _beliefsV;
    std::vector<Factor> _beliefsF;
    double _logZ;
    double _est_logZ;
    double _old_est_logZ;

    double _sum_level;
    size_t _num_leaves;

    boost::shared_ptr<std::ofstream> _clamp_ofstream;

    bbp_cfn_t BBP_cost_function() {
      return props.bbp_cfn;
    }

    void printDebugInfo();

    /// Called by 'run', and by itself. Implements the main algorithm. 
    /** Chooses a variable to clamp, recurses, combines the logZ and
     *  beliefs estimates of the children, and returns the improved
     *  estimates in \a lz_out and \a beliefs_out to its parent
     */
    void runRecurse(InfAlg *bp,
                    double orig_logZ,
                    std::vector<size_t> clamped_vars_list,
                    size_t &num_leaves,
                    size_t &choose_count,
                    double &sum_level,
                    Real &lz_out,
                    std::vector<Factor>& beliefs_out);

    /// Choose the next variable to clamp
    /** Choose the next variable to clamp, given a converged InfAlg (\a bp),
     *  and a vector of variables that are already clamped (\a
     *  clamped_vars_list). Returns the chosen variable in \a i, and
     *  the set of states in \a xis. If \a maxVarOut is non-NULL and
     *  props.choose==CHOOSE_BBP then it is used to store the
     *  adjoint of the chosen variable
     */
    virtual bool chooseNextClampVar(InfAlg* bp,
                                    std::vector<size_t> &clamped_vars_list,
                                    size_t &i, std::vector<size_t> &xis, Real *maxVarOut);

    /// Return the InfAlg to use at each step of the recursion. 
    /// \todo At present, only returns a BP instance
    InfAlg* getInfAlg();

    size_t _iters;
    double _maxdiff;

    void setBeliefs(const std::vector<Factor>& bs, double logZ) {
        size_t i=0;
        _beliefsV.clear(); _beliefsV.reserve(nrVars());
        _beliefsF.clear(); _beliefsF.reserve(nrFactors());
        for(i=0; i<nrVars(); i++) _beliefsV.push_back(bs[i]);
        for(; i<nrVars()+nrFactors(); i++) _beliefsF.push_back(bs[i]);
        _logZ = logZ;
    }

    void construct();

  public:

    //----------------------------------------------------------------

    /// construct CBP object from FactorGraph
    CBP(const FactorGraph &fg, const PropertySet &opts) : DAIAlgFG(fg) {
        props.set(opts);
        construct();
    }

    static const char *Name;

    /// @name General InfAlg interface
    //@{
    virtual CBP* clone() const { return new CBP(*this); }
    virtual CBP* create() const { DAI_THROW(NOT_IMPLEMENTED); }
    virtual std::string identify() const { return std::string(Name) + props.toString(); }
    virtual Factor belief (const Var &n) const { return _beliefsV[findVar(n)]; }
    virtual Factor belief (const VarSet &) const { DAI_THROW(NOT_IMPLEMENTED); }
    virtual std::vector<Factor> beliefs() const { return concat(_beliefsV, _beliefsF); }
    virtual Real logZ() const { return _logZ; }
    virtual void init() {};
    virtual void init( const VarSet & ) {};
    virtual double run();
    virtual double maxDiff() const { return _maxdiff; }
    virtual size_t Iterations() const { return _iters; }
    //@}

    Factor beliefV (size_t i) const { return _beliefsV[i]; }
    Factor beliefF (size_t I) const { return _beliefsF[I]; }

    //----------------------------------------------------------------

 public:
/* PROPERTIES(props,CBP) {
        typedef BP::Properties::UpdateType UpdateType;
        DAI_ENUM(RecurseType,REC_FIXED,REC_LOGZ,REC_BDIFF);
        DAI_ENUM(ChooseMethodType,CHOOSE_RANDOM,CHOOSE_MAXENT,CHOOSE_BBP,CHOOSE_BP_L1,CHOOSE_BP_CFN);
        DAI_ENUM(ClampType,CLAMP_VAR,CLAMP_FACTOR);
        
        size_t verbose = 0;

        /// Tolerance to use in BP
        double tol;
        /// Update style for BP
        UpdateType updates;
        /// Maximum number of iterations for BP
        size_t maxiter;

        /// Tolerance to use for controlling recursion depth (\a recurse
        /// is REC_LOGZ or REC_BDIFF)
        double rec_tol;
        /// Maximum number of levels of recursion (\a recurse is REC_FIXED)
        size_t max_levels = 10;
        /// If choose=CHOOSE_BBP and maximum adjoint is less than this value, don't recurse
        double min_max_adj;
        /// Heuristic for choosing clamping variable
        ChooseMethodType choose;
        /// Method for deciding when to stop recursing
        RecurseType recursion;
        /// Whether to clamp variables or factors
        ClampType clamp;
        /// Properties to pass to BBP
        PropertySet bbp_props;
        /// Cost function to use for BBP
        bbp_cfn_t bbp_cfn;
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
        typedef BP::Properties::UpdateType UpdateType;
        DAI_ENUM(RecurseType,REC_FIXED,REC_LOGZ,REC_BDIFF);
        DAI_ENUM(ChooseMethodType,CHOOSE_RANDOM,CHOOSE_MAXENT,CHOOSE_BBP,CHOOSE_BP_L1,CHOOSE_BP_CFN);
        DAI_ENUM(ClampType,CLAMP_VAR,CLAMP_FACTOR);
        size_t verbose;
        /// Tolerance to use in BP
        double tol;
        /// Update style for BP
        UpdateType updates;
        /// Maximum number of iterations for BP
        size_t maxiter;
        /// Tolerance to use for controlling recursion depth (\a recurse
        /// is REC_LOGZ or REC_BDIFF)
        double rec_tol;
        /// Maximum number of levels of recursion (\a recurse is REC_FIXED)
        size_t max_levels;
        /// If choose=CHOOSE_BBP and maximum adjoint is less than this value, don't recurse
        double min_max_adj;
        /// Heuristic for choosing clamping variable
        ChooseMethodType choose;
        /// Method for deciding when to stop recursing
        RecurseType recursion;
        /// Whether to clamp variables or factors
        ClampType clamp;
        /// Properties to pass to BBP
        PropertySet bbp_props;
        /// Cost function to use for BBP
        bbp_cfn_t bbp_cfn;
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

    Properties::ChooseMethodType ChooseMethod() { return props.choose; }
    Properties::RecurseType Recursion() { return props.recursion; }
    Properties::ClampType Clamping() { return props.clamp; }
    size_t maxClampLevel() { return props.max_levels; }
    double minMaxAdj() { return props.min_max_adj; }
    double recTol() { return props.rec_tol; }
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
