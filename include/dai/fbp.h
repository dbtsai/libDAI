/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2009 Frederik Eaton
 */


/// \file
/// \brief Defines class FBP, which implements Fractional Belief Propagation


#ifndef __defined_libdai_fbp_h
#define __defined_libdai_fbp_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>
#include <dai/bp.h>


namespace dai {


/// Approximate inference algorithm "Fractional Belief Propagation" [\ref WiH03]
/** The Fractional Belief Propagation algorithm is like Belief
 *  Propagation, but associates each factor with a scale parameter
 *  which controls the divergence measure being minimized. Standard
 *  Belief Propagation corresponds to the case of FBP where each scale
 *  parameter is 1. When cast as an EP algorithm, BP (and EP) minimize
 *  the inclusive KL-divergence, i.e. \f$\min_q KL(p||q)\f$ (note that the
 *  Bethe free energy is typically derived from \f$ KL(q||p) \f$). If each
 *  factor \a I has scale parameter \f$ c_I \f$, then FBP minimizes the
 *  alpha-divergence with \f$ \alpha=1/c_I \f$ for that factor, which also
 *  corresponds to Power EP [\ref Min05].
 *
 *  The messages \f$m_{I\to i}(x_i)\f$ are passed from factors \f$I\f$ to variables \f$i\f$. 
 *  The update equation is given by:
 *    \f[ m_{I\to i}(x_i) \propto \left( \sum_{x_{N_I\setminus\{i\}}} f_I(x_I)^{1/c_I} \prod_{j\in N_I\setminus\{i\}} m_{I\to j}^{1-1/c_I}\right)^{c_I} \prod_{J\in N_j\setminus\{I\}} m_{J\to j} \f]
 *  After convergence, the variable beliefs are calculated by:
 *    \f[ b_i(x_i) \propto \prod_{I\in N_i} m_{I\to i} \f]
 *  and the factor beliefs are calculated by:
 *    \f[ b_I(x_I) \propto f_I(x_I)^{1/c_I} \prod_{j \in N_I} m_{I\to j}^{1-1/c_I} \prod_{J\in N_j\setminus\{I\}} m_{J\to j} \f]
 *
 *  \todo Implement logZ
 *  \todo Why are the _scale_var necessary?
 *  \todo Add nice way to set scale parameters
 *
 *  \author Frederik Eaton
 */
class FBP : public BP {
    protected:
        /// Factor scale parameters (indexed by factor ID)
        std::vector<Real> _scale_factor;
        /// Variable scale parameters (indexed by variable ID)
        /** \note Equal to sum of scale parameters of neighboring factors
         */
        std::vector<Real> _scale_var;

    public:
        /// Name of this inference algorithm
        static const char *Name;

    public:
    /// \name Constructors/destructors
    //@{
        /// Default constructor
        FBP() : BP(), _scale_factor(), _scale_var() {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param opts Parameters @see BP::Properties
         */
        FBP( const FactorGraph &fg, const PropertySet &opts ) : BP(fg, opts), _scale_factor(), _scale_var() {
            setProperties( opts );
            construct();
        }
    //@}

    /// \name General InfAlg interface
    //@{
        virtual FBP* clone() const { return new FBP(*this); }
        virtual std::string identify() const;
    //@}

    /// \name FBP accessors/mutators for scale parameters
    //@{
        /// Returns scale parameter of the \a I 'th factor
        Real scaleF( size_t I ) const { return _scale_factor[I]; }

        /// Returns constant reference to vector of all factor scale parameters
        const std::vector<Real>& scaleFs() const { return _scale_factor; }

        /// Returns scale parameter of the \a i 'th variable
        Real scaleV( size_t i ) const { return _scale_var[i]; }

        /// Returns constant reference to vector of all variable scale parameters
        const std::vector<Real>& scaleVs() const { return _scale_var; }

        /// Sets the scale parameter of the \a I 'th factor to \a c
        void setScaleF( size_t I, Real c ) {
            _scale_factor[I] = c;
            foreach( const Neighbor &i, nbF(I) )
                recalcScaleV(i);
        }

        /// Sets the scale parameters of all factors simultaenously
        /** \note Faster than calling setScaleF(size_t,Real) for each factor
         */
        void setScaleFs( const std::vector<Real> &c ) {
            _scale_factor = c;
            recalcScaleVs();
        }

        /// Recalculates all variable scale parameters
        /** \note For each variable, its scale parameter is set to 
         *  the sum of the scale parameters of its neighboring factors.
         */
        void recalcScaleVs() {
            for( size_t i = 0; i < nrVars(); i++ )
                recalcScaleV(i);
        }

        /// Recalculates the scale parameter of the \a i 'th variable
        /** \note The scale parameter is set to the sum of the scale parameters of its neighboring factors.
         */
        void recalcScaleV( size_t i ) {
            // Set _scale_var[i] to the sum of its neighbors
            Real c_i = 0.0;
            foreach( const Neighbor &I, nbV(i) )
                c_i += scaleF(I);
            _scale_var[i] = c_i;
        }

    protected:
        // Calculate the updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        virtual void calcNewMessage( size_t i, size_t _I );

        // Calculates unnormalized belief of factor \a I
        virtual void calcBeliefF( size_t I, Prob &p ) const;

        // Helper function for constructors
        virtual void construct();

        /// (Re)constructs the scale parameters data structures
        void constructScaleParams() {
            _scale_factor.resize( nrFactors(), 1.0 );
            _scale_var.resize( nrVars() );
            recalcScaleVs();
        }
};


} // end of namespace dai


#endif
