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
/// \brief Defines class BBP [\ref EaG09]
/// \todo Improve documentation
/// \todo Clean up


#ifndef ___defined_libdai_bbp_h
#define ___defined_libdai_bbp_h


#include <vector>
#include <utility>

#include <dai/prob.h>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>

#include <dai/bp_dual.h>


namespace dai {


std::vector<Prob> get_zero_adj_F(const FactorGraph&);
std::vector<Prob> get_zero_adj_V(const FactorGraph&);


/// Implements BBP (Back-belief-propagation) [\ref EaG09]
class BBP {
protected:
    // ----------------------------------------------------------------
    // inputs
    BP_dual _bp_dual;
    const FactorGraph* _fg;
    const InfAlg *_ia;

    // iterations
    size_t _iters;

    // ----------------------------------------------------------------
    // Outputs
    std::vector<Prob> _adj_psi_V, _adj_psi_F;
    // The following vectors are indexed [i][_I]
    std::vector<std::vector<Prob> > _adj_n, _adj_m; 
    std::vector<Prob> _adj_b_V, _adj_b_F;

    // Helper quantities computed from the BP messages:
    // _T[i][_I]
    std::vector<std::vector<Prob > > _T;
    // _U[I][_i]
    std::vector<std::vector<Prob > > _U;
    // _S[i][_I][_j]
    std::vector<std::vector<std::vector<Prob > > > _S;
    // _R[I][_i][_J]
    std::vector<std::vector<std::vector<Prob > > > _R;

    std::vector<Prob> _adj_b_V_unnorm, _adj_b_F_unnorm;
    std::vector<Prob> _init_adj_psi_V;
    std::vector<Prob> _init_adj_psi_F;

    std::vector<std::vector<Prob> > _adj_n_unnorm, _adj_m_unnorm;
    std::vector<std::vector<Prob> > _new_adj_n, _new_adj_m;

    // ----------------------------------------------------------------
    // Indexing for performance
  
    /// Calculates _indices, which is a cache of IndexFor (see bp.cpp)
    void RegenerateInds();
    
    typedef std::vector<size_t>  _ind_t;
    std::vector<std::vector<_ind_t> >  _indices; 
    const _ind_t& _index(size_t i, size_t _I) const { return _indices[i][_I]; }

    // ----------------------------------------------------------------
    // Initialization

    /// Calculate T values (see paper)
    void RegenerateT();
    /// Calculate U values (see paper)
    void RegenerateU();
    /// Calculate S values (see paper)
    void RegenerateS();
    /// Calculate R values (see paper)
    void RegenerateR();
    /// Calculate _adj_b_V_unnorm and _adj_b_F_unnorm from _adj_b_V and _adj_b_F
    void RegenerateInputs();
    /// Initialise members for factor adjoints (call after RegenerateInputs)
    void RegeneratePsiAdjoints();
    /// Initialise members for messages adjoints (call after RegenerateInputs)
    void RegenerateParMessageAdjoints();
    /** Same as RegenerateMessageAdjoints, but calls sendSeqMsgN rather
     *  than updating _adj_n (and friends) which are unused in sequential algorithm
     */
    void RegenerateSeqMessageAdjoints(); 

    DAI_ACCMUT(Prob & T(size_t i, size_t _I), { return _T[i][_I]; });
    DAI_ACCMUT(Prob & U(size_t I, size_t _i), { return _U[I][_i]; });
    DAI_ACCMUT(Prob & S(size_t i, size_t _I, size_t _j), { return _S[i][_I][_j]; });
    DAI_ACCMUT(Prob & R(size_t I, size_t _i, size_t _J), { return _R[I][_i][_J]; });

    void calcNewN(size_t i, size_t _I);
    void calcNewM(size_t i, size_t _I);
    void calcUnnormMsgM(size_t i, size_t _I);
    void calcUnnormMsgN(size_t i, size_t _I);
    void upMsgM(size_t i, size_t _I);
    void upMsgN(size_t i, size_t _I);
    void doParUpdate();
    Real getUnMsgMag();
    void getMsgMags(Real &s, Real &new_s);

    void zero_adj_b_F() {
        _adj_b_F.clear();
        _adj_b_F.reserve(_fg->nrFactors());
        for(size_t I=0; I<_fg->nrFactors(); I++) {
            _adj_b_F.push_back(Prob(_fg->factor(I).states(),Real(0.0)));
        }
    }

    //----------------------------------------------------------------
    // new interface

    void incrSeqMsgM(size_t i, size_t _I, const Prob& p);
    void updateSeqMsgM(size_t i, size_t _I);
    void sendSeqMsgN(size_t i, size_t _I, const Prob &f);
    void sendSeqMsgM(size_t i, size_t _I);
    /// used instead of upMsgM / calcNewM, calculates adj_m_unnorm as well
    void setSeqMsgM(size_t i, size_t _I, const Prob &p); 

    Real getMaxMsgM();
    Real getTotalMsgM();
    Real getTotalNewMsgM();
    Real getTotalMsgN();

    void getArgmaxMsgM(size_t &i, size_t &_I, Real &mag);

public:
    /// Called by 'init', recalculates intermediate values
    void Regenerate();

    BBP(const InfAlg *ia, const PropertySet &opts) :
        _bp_dual(ia), _fg(&(ia->fg())), _ia(ia)
    {
        props.set(opts);
    }

    void init(const std::vector<Prob> &adj_b_V, const std::vector<Prob> &adj_b_F,
              const std::vector<Prob> &adj_psi_V, const std::vector<Prob> &adj_psi_F) {
        _adj_b_V = adj_b_V;
        _adj_b_F = adj_b_F;
        _init_adj_psi_V = adj_psi_V;
        _init_adj_psi_F = adj_psi_F;
        Regenerate(); 
    }
    void init(const std::vector<Prob> &adj_b_V, const std::vector<Prob> &adj_b_F) {
        init(adj_b_V, adj_b_F, get_zero_adj_V(*_fg), get_zero_adj_F(*_fg));
    }
    void init(const std::vector<Prob> &adj_b_V) {
        init(adj_b_V, get_zero_adj_F(*_fg));
    }

    /// run until change is less than given tolerance
    void run();

    size_t doneIters() { return _iters; }

    DAI_ACCMUT(Prob& adj_psi_V(size_t i), { return _adj_psi_V[i]; });
    DAI_ACCMUT(Prob& adj_psi_F(size_t I), { return _adj_psi_F[I]; });
    DAI_ACCMUT(Prob& adj_b_V(size_t i), { return _adj_b_V[i]; });
    DAI_ACCMUT(Prob& adj_b_F(size_t I), { return _adj_b_F[I]; });
 protected:
    DAI_ACCMUT(Prob& adj_n(size_t i, size_t _I), { return _adj_n[i][_I]; });
    DAI_ACCMUT(Prob& adj_m(size_t i, size_t _I), { return _adj_m[i][_I]; });
 public: 

    /// Parameters of this inference algorithm
/* PROPERTIES(props,BBP) {
       DAI_ENUM(UpdateType,SEQ_FIX,SEQ_MAX,SEQ_BP_REV,SEQ_BP_FWD,PAR);
       size_t verbose;
       /// tolerance (not used for updates=SEQ_BP_{REV,FWD})
       double tol;
       size_t maxiter;
       /// damping (0 for none)
       double damping;
       UpdateType updates;
       bool clean_updates;
    }
*/
/* {{{ GENERATED CODE: DO NOT EDIT. Created by 
    ./scripts/regenerate-properties include/dai/bbp.h src/bbp.cpp 
*/
    struct Properties {
        DAI_ENUM(UpdateType,SEQ_FIX,SEQ_MAX,SEQ_BP_REV,SEQ_BP_FWD,PAR);
        size_t verbose;
        /// tolerance (not used for updates=SEQ_BP_{REV,FWD})
        double tol;
        size_t maxiter;
        /// damping (0 for none)
        double damping;
        UpdateType updates;
        bool clean_updates;

        /// Set members from PropertySet
        void set(const PropertySet &opts);
        /// Get members into PropertySet
        PropertySet get() const;
        /// Convert to a string which can be parsed as a PropertySet
        std::string toString() const;
    } props;
/* }}} END OF GENERATED CODE */
};

/// Cost functions. Not used by BBP class, only used by following functions.
DAI_ENUM(bbp_cfn_t,cfn_gibbs_b,cfn_gibbs_b2,cfn_gibbs_exp,cfn_gibbs_b_factor,cfn_gibbs_b2_factor,cfn_gibbs_exp_factor,cfn_var_ent,cfn_factor_ent,cfn_bethe_ent);

/// Initialise BBP using InfAlg, cost function, and stateP
/** Calls bbp.init with adjoints calculated from ia.beliefV and
 *  ia.beliefF. stateP is a Gibbs state and can be NULL, it will be
 *  initialised using a Gibbs run of 2*fg.Iterations() iterations.
 */
void initBBPCostFnAdj(BBP& bbp, const InfAlg& ia, bbp_cfn_t cfn_type, const std::vector<size_t>* stateP);

/// Answers question: Does given cost function depend on having a Gibbs state?
bool needGibbsState(bbp_cfn_t cfn);

/// Calculate actual value of cost function (cfn_type, stateP)
/** This function returns the actual value of the cost function whose
 *  gradient with respect to singleton beliefs is given by
 *  gibbsToB1Adj on the same arguments
 */
Real getCostFn(const InfAlg& fg, bbp_cfn_t cfn_type, const std::vector<size_t> *stateP);

/// Function to test the validity of adjoints computed by BBP
/** given a state for each variable, use numerical derivatives
 *  (multiplying a factor containing a variable by psi_1 adjustments)
 *  to verify accuracy of _adj_psi_V.
 *  'h' controls size of perturbation.
 *  'bbpTol' controls tolerance of BBP run.
 */
double numericBBPTest(const InfAlg& bp, const std::vector<size_t> *state, const PropertySet& bbp_props, bbp_cfn_t cfn, double h);

// ----------------------------------------------------------------
// Utility functions, some of which are used elsewhere

/// Subtract 1 from a size_t, or return 0 if the argument is 0
inline size_t oneLess(size_t v) { return v==0?v:v-1; }

/// function to compute adj_w_unnorm from w, Z_w, adj_w
Prob unnormAdjoint(const Prob &w, Real Z_w, const Prob &adj_w);

/// Runs Gibbs sampling for 'iters' iterations on ia.fg(), and returns state
std::vector<size_t> getGibbsState(const InfAlg& ia, size_t iters);


} // end of namespace dai


#endif
