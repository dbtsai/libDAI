/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010 Joris Mooij
 */


/// \file
/// \brief Defines class TRWBP, which implements Tree-Reweighted Belief Propagation


#ifndef __defined_libdai_trwbp_h
#define __defined_libdai_trwbp_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>
#include <dai/bp.h>


namespace dai {


/// Approximate inference algorithm "Tree-Reweighted Belief Propagation" [\ref WJW03]
/** The Tree-Reweighted Belief Propagation algorithm is like Belief
 *  Propagation, but associates each factor with a scale parameter.
 *  which controls the divergence measure being minimized.
 *
 *  The messages \f$m_{I\to i}(x_i)\f$ are passed from factors \f$I\f$ to variables \f$i\f$. 
 *  The update equation is given by:
 *    \f[ m_{I\to i}(x_i) \propto \sum_{x_{N_I\setminus\{i\}}} f_I(x_I)^{1/c_I} \prod_{j\in N_I\setminus\{i\}} m_{I\to j}^{c_I-1} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}^{c_J} \f]
 *  After convergence, the variable beliefs are calculated by:
 *    \f[ b_i(x_i) \propto \prod_{I\in N_i} m_{I\to i}^{c_I} \f]
 *  and the factor beliefs are calculated by:
 *    \f[ b_I(x_I) \propto f_I(x_I)^{1/c_I} \prod_{j \in N_I} m_{I\to j}^{c_I-1} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}^{c_J} \f]
 *
 *  \todo Fix documentation
 */
class TRWBP : public BP {
    protected:
        /// Factor scale parameters (indexed by factor ID)
        std::vector<Real> _edge_weight;

    public:
        /// Name of this inference algorithm
        static const char *Name;

    public:
    /// \name Constructors/destructors
    //@{
        /// Default constructor
        TRWBP() : BP(), _edge_weight() {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param opts Parameters @see BP::Properties
         */
        TRWBP( const FactorGraph &fg, const PropertySet &opts ) : BP(fg, opts), _edge_weight() {
            setProperties( opts );
            construct();
        }
    //@}

    /// \name General InfAlg interface
    //@{
        virtual TRWBP* clone() const { return new TRWBP(*this); }
        virtual std::string identify() const;
        virtual Real logZ() const;
    //@}

    /// \name TRWBP accessors/mutators for scale parameters
    //@{
        /// Returns scale parameter of edge corresponding to the \a I 'th factor
        Real edgeWeight( size_t I ) const { return _edge_weight[I]; }

        /// Returns constant reference to vector of all factor scale parameters
        const std::vector<Real>& edgeWeights() const { return _edge_weight; }

        /// Sets the scale parameter of the \a I 'th factor to \a c
        void setEdgeWeight( size_t I, Real c ) { _edge_weight[I] = c; }

        /// Sets the scale parameters of all factors simultaenously
        /** \note Faster than calling setScaleF(size_t,Real) for each factor
         */
        void setEdgeWeights( const std::vector<Real> &c ) { _edge_weight = c; }

    protected:
        // Calculate the updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        virtual void calcNewMessage( size_t i, size_t _I );

        /// Calculates unnormalized belief of variable \a i
        virtual void calcBeliefV( size_t i, Prob &p ) const;

        // Calculates unnormalized belief of factor \a I
        virtual void calcBeliefF( size_t I, Prob &p ) const;

        // Helper function for constructors
        virtual void construct();
};


} // end of namespace dai


#endif
