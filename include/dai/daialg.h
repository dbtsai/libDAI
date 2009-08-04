/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

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
/// \brief Defines abstract base class InfAlg, its descendants DAIAlg<T>, the specializations DAIAlgFG and DAIAlgRG and some generic inference methods.
/// \todo Improve documentation


#ifndef __defined_libdai_daialg_h
#define __defined_libdai_daialg_h


#include <string>
#include <iostream>
#include <vector>
#include <dai/factorgraph.h>
#include <dai/regiongraph.h>


namespace dai {


/// InfAlg is an abstract base class, defining the common interface of all inference algorithms in libDAI.
/** \todo General marginalization functions like calcMarginal now copy a complete InfAlg object. Instead, 
 *  it would make more sense that they construct a new object without copying the FactorGraph or RegionGraph. 
 *  Or they can simply be made methods of the general InfAlg class.
 *  \idea Use a PropertySet as output of an InfAlg, instead of functions like maxDiff() and Iterations().
 */
class InfAlg {
    public:
        /// Virtual desctructor (needed because this class contains virtual functions)
        virtual ~InfAlg() {}

    public:
        /// Returns a pointer to a new, cloned copy of *this (i.e., virtual copy constructor)
        virtual InfAlg* clone() const = 0;

        /// Identifies itself for logging purposes
        virtual std::string identify() const = 0;

        /// Returns the "belief" (i.e., approximate marginal probability distribution) of a variable
        virtual Factor belief( const Var &n ) const = 0;

        /// Returns the "belief" (i.e., approximate marginal probability distribution) of a set of variables
        virtual Factor belief( const VarSet &n ) const = 0;

        /// Returns marginal for a variable.
        /** Sometimes preferred to belief() for performance reasons.
          * Faster implementations exist in e.g. BP.
          */
        virtual Factor beliefV( size_t i ) const { return belief( fg().var(i) ); }

        /// Returns marginal for a factor.
        /** Sometimes preferred to belief() for performance reasons.
          * Faster implementations exist in e.g. BP.
          */
        virtual Factor beliefF( size_t I ) const { return belief( fg().factor(I).vars() ); }

        /// Returns all "beliefs" (i.e., approximate marginal probability distribution) calculated by the algorithm
        virtual std::vector<Factor> beliefs() const = 0;

        /// Returns the logarithm of the (approximated) partition sum (normalizing constant of the factor graph)
        virtual Real logZ() const = 0;

        /// Initializes all data structures of the approximate inference algorithm
        /** This method should be called at least once before run() is called
         */
        virtual void init() = 0;

        /// Initializes all data structures corresponding to some set of variables
        /** This method can be used to do a partial initialization after a part of the factor graph has changed.
         *  Instead of initializing all data structures, it only initializes those involving the variables in ns.
         */
        virtual void init( const VarSet &ns ) = 0;

        /// Runs the approximate inference algorithm
        /*  Before run() is called the first time, init() should be called.
         *  If run() returns successfully, the results can be queried using the methods belief(), beliefs() and logZ().
         */
        virtual double run() = 0;

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$)
        virtual void clamp( const Var & n, size_t i, bool backup = false ) = 0;

        /// Set all factors interacting with var(i) to 1
        virtual void makeCavity( size_t i, bool backup = false ) = 0;

        /// Return maximum difference between single node beliefs in the last pass
        /// \throw Exception if not implemented/supported
        virtual double maxDiff() const = 0;

        /// Return number of passes over the factorgraph
        /// \throw Exception if not implemented/supported
        virtual size_t Iterations() const = 0;


        /// Get reference to underlying FactorGraph
        virtual FactorGraph &fg() = 0;

        /// Get const reference to underlying FactorGraph
        virtual const FactorGraph &fg() const = 0;

        /// Save factor I
        virtual void backupFactor( size_t I ) = 0;
        /// Save Factors involving ns
        virtual void backupFactors( const VarSet &ns ) = 0;

        /// Restore factor I
        virtual void restoreFactor( size_t I ) = 0;
        /// Restore Factors involving ns
        virtual void restoreFactors( const VarSet &ns ) = 0;
};


/// Combines an InfAlg and a graphical model, e.g., a FactorGraph or RegionGraph
/** \tparam GRM Should be castable to FactorGraph
 *  \todo A DAIAlg should not inherit from a FactorGraph or RegionGraph, but should
 *  store a reference to the graphical model object. This prevents needless copying 
 *  of (possibly large) data structures. Disadvantage: the caller must not change 
 *  the graphical model between calls to the inference algorithm (maybe a smart_ptr 
 *  or some locking mechanism would help here?). 
 */
template <class GRM>
class DAIAlg : public InfAlg, public GRM {
    public:
        /// Default constructor
        DAIAlg() : InfAlg(), GRM() {}
        
        /// Construct from GRM 
        DAIAlg( const GRM &grm ) : InfAlg(), GRM(grm) {}

        /// Save factor I
        void backupFactor( size_t I ) { GRM::backupFactor( I ); }
        /// Save Factors involving ns
        void backupFactors( const VarSet &ns ) { GRM::backupFactors( ns ); }

        /// Restore factor I
        void restoreFactor( size_t I ) { GRM::restoreFactor( I ); }
        /// Restore Factors involving ns
        void restoreFactors( const VarSet &ns ) { GRM::restoreFactors( ns ); }

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$)
        void clamp( const Var & n, size_t i, bool backup = false ) { GRM::clamp( n, i, backup ); }

        /// Set all factors interacting with var(i) to 1
        void makeCavity( size_t i, bool backup = false ) { GRM::makeCavity( i, backup ); }

        /// Get reference to underlying FactorGraph
        FactorGraph &fg() { return (FactorGraph &)(*this); }

        /// Get const reference to underlying FactorGraph
        const FactorGraph &fg() const { return (const FactorGraph &)(*this); }
};


/// Base class for inference algorithms that operate on a FactorGraph
typedef DAIAlg<FactorGraph> DAIAlgFG;

/// Base class for inference algorithms that operate on a RegionGraph
typedef DAIAlg<RegionGraph> DAIAlgRG;


Factor calcMarginal( const InfAlg & obj, const VarSet & ns, bool reInit );
std::vector<Factor> calcPairBeliefs( const InfAlg & obj, const VarSet& ns, bool reInit );
std::vector<Factor> calcPairBeliefsNew( const InfAlg & obj, const VarSet& ns, bool reInit );
Factor calcMarginal2ndO( const InfAlg & obj, const VarSet& ns, bool reInit );


} // end of namespace dai


#endif
