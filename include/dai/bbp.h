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


/// Computes the adjoint of the unnormed probability vector from the normalizer and the adjoint of the normalized probability vector @see eqn. (13) in [\ref EaG09]
Prob unnormAdjoint( const Prob &w, Real Z_w, const Prob &adj_w );

/// Runs Gibbs sampling for \a iters iterations on ia.fg(), and returns state
std::vector<size_t> getGibbsState( const InfAlg &ia, size_t iters );


/// Implements BBP (Back-Belief-Propagation) [\ref EaG09]
class BBP {
    protected:
        /// @name Inputs
        //@{
        BP_dual _bp_dual;
        const FactorGraph *_fg;
        const InfAlg *_ia;
        //@}

        /// Number of iterations done
        size_t _iters;

        /// @name Outputs
        //@{
        /// Variable factor adjoints
        std::vector<Prob> _adj_psi_V;
        /// Factor adjoints
        std::vector<Prob> _adj_psi_F;
        /// Variable->factor message adjoints (indexed [i][_I])
        std::vector<std::vector<Prob> > _adj_n;
        /// Factor->variable message adjoints (indexed [i][_I])
        std::vector<std::vector<Prob> > _adj_m;
        /// Normalized variable belief adjoints
        std::vector<Prob> _adj_b_V;
        /// Normalized factor belief adjoints
        std::vector<Prob> _adj_b_F;
        //@}

        /// @name Helper quantities computed from the BP messages
        //@{
        /// _T[i][_I] (see eqn. (41) in [\ref EaG09])
        std::vector<std::vector<Prob > > _T;
        /// _U[I][_i] (see eqn. (42) in [\ref EaG09])
        std::vector<std::vector<Prob > > _U;
        /// _S[i][_I][_j] (see eqn. (43) in [\ref EaG09])
        std::vector<std::vector<std::vector<Prob > > > _S;
        /// _R[I][_i][_J] (see eqn. (44) in [\ref EaG09])
        std::vector<std::vector<std::vector<Prob > > > _R;
        //@}

        /// Unnormalized variable belief adjoints
        std::vector<Prob> _adj_b_V_unnorm;
        /// Unnormalized factor belief adjoints
        std::vector<Prob> _adj_b_F_unnorm;

        /// Initial variable factor adjoints
        std::vector<Prob> _init_adj_psi_V;
        /// Initial factor adjoints
        std::vector<Prob> _init_adj_psi_F;

        /// Unnormalized variable->factor message adjoint (indexed [i][_I])
        std::vector<std::vector<Prob> > _adj_n_unnorm;
        /// Unnormalized factor->variable message adjoint (indexed [i][_I])
        std::vector<std::vector<Prob> > _adj_m_unnorm;
        /// Updated normalized variable->factor message adjoint (indexed [i][_I])
        std::vector<std::vector<Prob> > _new_adj_n;
        /// Updated normalized factor->variable message adjoint (indexed [i][_I])
        std::vector<std::vector<Prob> > _new_adj_m;

        /// @name Optimized indexing (for performance)
        //@{
        /// Calculates _indices, which is a cache of IndexFor @see bp.cpp
        void RegenerateInds();
        
        /// Index type
        typedef std::vector<size_t>  _ind_t;
        /// Cached indices (indexed [i][_I])
        std::vector<std::vector<_ind_t> >  _indices; 
        /// Returns an index from the cache
        const _ind_t& _index(size_t i, size_t _I) const { return _indices[i][_I]; }
        //@}

        /// @name Initialization
        //@{
        /// Calculate T values; see eqn. (41) in [\ref EaG09]
        void RegenerateT();
        /// Calculate U values; see eqn. (42) in [\ref EaG09]
        void RegenerateU();
        /// Calculate S values; see eqn. (43) in [\ref EaG09]
        void RegenerateS();
        /// Calculate R values; see eqn. (44) in [\ref EaG09]
        void RegenerateR();
        /// Calculate _adj_b_V_unnorm and _adj_b_F_unnorm from _adj_b_V and _adj_b_F
        void RegenerateInputs();
        /// Initialise members for factor adjoints (call after RegenerateInputs)
        void RegeneratePsiAdjoints();
        /// Initialise members for message adjoints (call after RegenerateInputs) for parallel algorithm
        void RegenerateParMessageAdjoints();
        /// Initialise members for message adjoints (call after RegenerateInputs) for sequential algorithm
        /** Same as RegenerateMessageAdjoints, but calls sendSeqMsgN rather
         *  than updating _adj_n (and friends) which are unused in the sequential algorithm.
         */
        void RegenerateSeqMessageAdjoints();
        //@}

        /// Returns T value; see eqn. (41) in [\ref EaG09]
        DAI_ACCMUT(Prob & T(size_t i, size_t _I), { return _T[i][_I]; });
        /// Retunrs U value; see eqn. (42) in [\ref EaG09]
        DAI_ACCMUT(Prob & U(size_t I, size_t _i), { return _U[I][_i]; });
        /// Returns S value; see eqn. (43) in [\ref EaG09]
        DAI_ACCMUT(Prob & S(size_t i, size_t _I, size_t _j), { return _S[i][_I][_j]; });
        /// Returns R value; see eqn. (44) in [\ref EaG09]
        DAI_ACCMUT(Prob & R(size_t I, size_t _i, size_t _J), { return _R[I][_i][_J]; });

        /// @name Parallel algorithm
        //@{
        /// Calculates new variable->factor message adjoint
        /** Increases variable factor adjoint according to eqn. (27) in [\ref EaG09] and
         *  calculates the new variable->factor message adjoint according to eqn. (29) in [\ref EaG09].
         */
        void calcNewN( size_t i, size_t _I );
        /// Calculates new factor->variable message adjoint
        /** Increases factor adjoint according to eqn. (28) in [\ref EaG09] and
         *  calculates the new factor->variable message adjoint according to the r.h.s. of eqn. (30) in [\ref EaG09].
         */
        void calcNewM( size_t i, size_t _I );
        /// Calculates unnormalized variable->factor message adjoint from the normalized one
        void calcUnnormMsgN( size_t i, size_t _I );
        /// Calculates unnormalized factor->variable message adjoint from the normalized one
        void calcUnnormMsgM( size_t i, size_t _I );
        /// Updates (un)normalized variable->factor message adjoints
        void upMsgN( size_t i, size_t _I );
        /// Updates (un)normalized factor->variable message adjoints
        void upMsgM( size_t i, size_t _I );
        /// Do one parallel update of all message adjoints
        void doParUpdate();
        //@}

        /// @name Sequential algorithm
        //@{
        /// Helper function for sendSeqMsgM: increases factor->variable message adjoint by p and calculates the corresponding unnormalized adjoint
        void incrSeqMsgM( size_t i, size_t _I, const Prob& p );
        //  DISABLED BECAUSE IT IS BUGGY:
        //  void updateSeqMsgM( size_t i, size_t _I );
        /// Sets normalized factor->variable message adjoint and calculates the corresponding unnormalized adjoint
        void setSeqMsgM( size_t i, size_t _I, const Prob &p ); 
        /// Implements routine Send-n in Figure 5 in [\ref EaG09]
        void sendSeqMsgN( size_t i, size_t _I, const Prob &f );
        /// Implements routine Send-m in Figure 5 in [\ref EaG09]
        void sendSeqMsgM( size_t i, size_t _I );
        //@}

        /// Calculates averaged L-1 norm of unnormalized message adjoints
        Real getUnMsgMag();
        /// Calculates averaged L-1 norms of current and new normalized message adjoints
        void getMsgMags( Real &s, Real &new_s );

        /// Sets all vectors _adj_b_F to zero
        void zero_adj_b_F() {
            _adj_b_F.clear();
            _adj_b_F.reserve( _fg->nrFactors() );
            for( size_t I = 0; I < _fg->nrFactors(); I++ )
                _adj_b_F.push_back( Prob( _fg->factor(I).states(), Real( 0.0 ) ) );
        }

        /// Returns indices and magnitude of the largest normalized factor->variable message adjoint
        void getArgmaxMsgM( size_t &i, size_t &_I, Real &mag );
        /// Returns magnitude of the largest (in L1-norm) normalized factor->variable message adjoint
        Real getMaxMsgM();
        /// Calculates sum of L1 norms of all normalized factor->variable message adjoints
        Real getTotalMsgM();
        /// Calculates sum of L1 norms of all updated normalized factor->variable message adjoints
        Real getTotalNewMsgM();
        /// Calculates sum of L1 norms of all normalized variable->factor message adjoints
        Real getTotalMsgN();

    public:
        /// Called by \a init, recalculates intermediate values
        void Regenerate();

        /// Constructor
        BBP( const InfAlg *ia, const PropertySet &opts ) : _bp_dual(ia), _fg(&(ia->fg())), _ia(ia) {
            props.set(opts);
        }

        /// Returns a vector of Probs (filled with zeroes) with state spaces corresponding to the factors in the factor graph fg
        std::vector<Prob> getZeroAdjF( const FactorGraph &fg );
        /// Returns a vector of Probs (filled with zeroes) with state spaces corresponding to the variables in the factor graph fg
        std::vector<Prob> getZeroAdjV( const FactorGraph &fg );

        /// Initializes belief adjoints and initial factor adjoints and regenerates
        void init( const std::vector<Prob> &adj_b_V, const std::vector<Prob> &adj_b_F, const std::vector<Prob> &adj_psi_V, const std::vector<Prob> &adj_psi_F ) {
            _adj_b_V = adj_b_V;
            _adj_b_F = adj_b_F;
            _init_adj_psi_V = adj_psi_V;
            _init_adj_psi_F = adj_psi_F;
            Regenerate(); 
        }

        /// Initializes belief adjoints and with zero initial factor adjoints and regenerates
        void init( const std::vector<Prob> &adj_b_V, const std::vector<Prob> &adj_b_F ) {
            init( adj_b_V, adj_b_F, getZeroAdjV(*_fg), getZeroAdjF(*_fg) );
        }

        /// Initializes variable belief adjoints (and sets factor belief adjoints to zero) and with zero initial factor adjoints and regenerates
        void init( const std::vector<Prob> &adj_b_V ) {
            init(adj_b_V, getZeroAdjF(*_fg));
        }

        /// Run until change is less than given tolerance
        void run();

        /// Return number of iterations done so far
        size_t doneIters() { return _iters; }

        /// Returns variable factor adjoint
        DAI_ACCMUT(Prob& adj_psi_V(size_t i), { return _adj_psi_V[i]; });
        /// Returns factor adjoint
        DAI_ACCMUT(Prob& adj_psi_F(size_t I), { return _adj_psi_F[I]; });
        /// Returns variable belief adjoint
        DAI_ACCMUT(Prob& adj_b_V(size_t i), { return _adj_b_V[i]; });
        /// Returns factor belief adjoint
        DAI_ACCMUT(Prob& adj_b_F(size_t I), { return _adj_b_F[I]; });

     protected:
        /// Returns variable->factor message adjoint
        DAI_ACCMUT(Prob& adj_n(size_t i, size_t _I), { return _adj_n[i][_I]; });
        /// Returns factor->variable message adjoint
        DAI_ACCMUT(Prob& adj_m(size_t i, size_t _I), { return _adj_m[i][_I]; });

     public: 
        /// Parameters of this algorithm
        /* PROPERTIES(props,BBP) {
           /// Enumeration of possible update schedules
           DAI_ENUM(UpdateType,SEQ_FIX,SEQ_MAX,SEQ_BP_REV,SEQ_BP_FWD,PAR);

           /// Verbosity
           size_t verbose;

           /// Maximum number of iterations
           size_t maxiter;

           /// Tolerance (not used for updates = SEQ_BP_REV, SEQ_BP_FWD)
           double tol;

           /// Damping constant (0 for none); damping = 1 - lambda where lambda is the damping constant used in [\ref EaG09]
           double damping;

           /// Update schedule
           UpdateType updates;

           // DISABLED BECAUSE IT IS BUGGY:
           // bool clean_updates;
        } 
        */
/* {{{ GENERATED CODE: DO NOT EDIT. Created by 
    ./scripts/regenerate-properties include/dai/bbp.h src/bbp.cpp 
*/
        struct Properties {
            /// Enumeration of possible update schedules
            DAI_ENUM(UpdateType,SEQ_FIX,SEQ_MAX,SEQ_BP_REV,SEQ_BP_FWD,PAR);
            /// Verbosity
            size_t verbose;
            /// Maximum number of iterations
            size_t maxiter;
            /// Tolerance (not used for updates = SEQ_BP_REV, SEQ_BP_FWD)
            double tol;
            /// Damping constant (0 for none); damping = 1 - lambda where lambda is the damping constant used in [\ref EaG09]
            double damping;
            /// Update schedule
            UpdateType updates;

            /// Set members from PropertySet
            void set(const PropertySet &opts);
            /// Get members into PropertySet
            PropertySet get() const;
            /// Convert to a string which can be parsed as a PropertySet
            std::string toString() const;
        } props;
/* }}} END OF GENERATED CODE */
};


/// Enumeration of several cost functions that can be used with BBP.
DAI_ENUM(bbp_cfn_t,CFN_GIBBS_B,CFN_GIBBS_B2,CFN_GIBBS_EXP,CFN_GIBBS_B_FACTOR,CFN_GIBBS_B2_FACTOR,CFN_GIBBS_EXP_FACTOR,CFN_VAR_ENT,CFN_FACTOR_ENT,CFN_BETHE_ENT);

/// Initialise BBP using InfAlg, cost function, and stateP
/** Calls bbp.init with adjoints calculated from ia.beliefV and
 *  ia.beliefF. stateP is a Gibbs state and can be NULL, it will be
 *  initialised using a Gibbs run of 2*fg.Iterations() iterations.
 */
void initBBPCostFnAdj( BBP &bbp, const InfAlg &ia, bbp_cfn_t cfn_type, const std::vector<size_t> *stateP );

/// Answers question: does the given cost function depend on having a Gibbs state?
bool needGibbsState( bbp_cfn_t cfn );

/// Calculate actual value of cost function (cfn_type, stateP)
/** This function returns the actual value of the cost function whose
 *  gradient with respect to singleton beliefs is given by
 *  gibbsToB1Adj on the same arguments
 */
Real getCostFn( const InfAlg &fg, bbp_cfn_t cfn_type, const std::vector<size_t> *stateP );

/// Function to test the validity of adjoints computed by BBP given a state for each variable using numerical derivatives.
/** Factors containing a variable are multiplied by psi_1 adjustments to verify accuracy of _adj_psi_V.
 *  \param bp BP object.
 *  \param state Global state of all variables.
 *  \param bbp_props BBP Properties.
 *  \param cfn Cost function to be used.
 *  \param h controls size of perturbation.
 */
double numericBBPTest( const InfAlg &bp, const std::vector<size_t> *state, const PropertySet &bbp_props, bbp_cfn_t cfn, double h );


} // end of namespace dai


#endif
