/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
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


#ifndef __defined_libdai_daialg_h
#define __defined_libdai_daialg_h


#include <string>
#include <iostream>
#include <vector>
#include <dai/factorgraph.h>
#include <dai/regiongraph.h>


namespace dai {


/// The InfAlg class is the common denominator of the various approximate inference algorithms.
/// A InfAlg object represents a discrete factorized probability distribution over multiple variables 
/// together with an inference algorithm.
class InfAlg {
    public:
        /// Clone (virtual copy constructor)
        virtual InfAlg* clone() const = 0;

        /// Virtual desctructor
        // (this is needed because this class contains virtual functions)
        virtual ~InfAlg() {}
        
        /// Identifies itself for logging purposes
        virtual std::string identify() const = 0;

        /// Get single node belief
        virtual Factor belief( const Var &n ) const = 0;

        /// Get general belief
        virtual Factor belief( const VarSet &n ) const = 0;

        /// Get all beliefs
        virtual std::vector<Factor> beliefs() const = 0;

        /// Get log partition sum
        virtual Complex logZ() const = 0;

        /// Clear messages and beliefs
        virtual void init() = 0;

        /// The actual approximate inference algorithm
        virtual double run() = 0;

        /// Save factor I
        virtual void saveProb( size_t I ) = 0;
        /// Save Factors involving ns
        virtual void saveProbs( const VarSet &ns ) = 0;

        /// Restore factor I
        virtual void undoProb( size_t I ) = 0;
        /// Restore Factors involving ns
        virtual void undoProbs( const VarSet &ns ) = 0;

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$)
        virtual void clamp( const Var & n, size_t i ) = 0;

        /// Return all variables that interact with var(i)
        virtual VarSet delta( size_t i ) const = 0;

        /// Set all factors interacting with var(i) to 1
        virtual void makeCavity( size_t i ) = 0;

        /// Get index of variable n
        virtual size_t findVar( const Var & n ) const = 0;

        /// Get index of first factor involving ns
        virtual size_t findFactor( const VarSet &ns ) const = 0;

        /// Get number of variables
        virtual size_t nrVars() const = 0;

        /// Get number of factors
        virtual size_t nrFactors() const = 0;

        /// Get const reference to variable i
        virtual const Var & var(size_t i) const = 0;

        /// Get reference to variable i
        virtual Var & var(size_t i) = 0;

        /// Get const reference to factor I
        virtual const Factor & factor( size_t I ) const = 0;

        /// Get reference to factor I
        virtual Factor & factor( size_t I ) = 0;

        /// Factor I has been updated
        virtual void updatedFactor( size_t I ) = 0;

        /// Return maximum difference between beliefs in the last pass
        virtual double maxDiff() const = 0;
};


template <class T>
class DAIAlg : public InfAlg, public T {
    public:
        /// Default constructor
        DAIAlg() : InfAlg(), T() {}
        
        /// Construct from T
        DAIAlg( const T &t ) : InfAlg(), T(t) {}

        /// Copy constructor
        DAIAlg( const DAIAlg & x ) : InfAlg(x), T(x) {}

        /// Save factor I (using T::saveProb)
        void saveProb( size_t I ) { T::saveProb( I ); }
        /// Save Factors involving ns (using T::saveProbs)
        void saveProbs( const VarSet &ns ) { T::saveProbs( ns ); }

        /// Restore factor I (using T::undoProb)
        void undoProb( size_t I ) { T::undoProb( I ); }
        /// Restore Factors involving ns (using T::undoProbs)
        void undoProbs( const VarSet &ns ) { T::undoProbs( ns ); }

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$) (using T::clamp)
        void clamp( const Var & n, size_t i ) { T::clamp( n, i ); }

        /// Return all variables that interact with var(i) (using T::delta)
        VarSet delta( size_t i ) const { return T::delta( i ); }

        /// Set all factors interacting with var(i) to 1 (using T::makeCavity)
        void makeCavity( size_t i ) { T::makeCavity( i ); }

        /// Get index of variable n (using T::findVar)
        size_t findVar( const Var & n ) const { return T::findVar(n); }

        /// Get index of first factor involving ns (using T::findFactor)
        size_t findFactor( const VarSet &ns ) const { return T::findFactor(ns); }

        /// Get number of variables (using T::nrFactors)
        size_t nrVars() const { return T::nrVars(); }

        /// Get number of factors (using T::nrFactors)
        size_t nrFactors() const { return T::nrFactors(); }

        /// Get const reference to variable i (using T::var)
        const Var & var( size_t i ) const { return T::var(i); }

        /// Get reference to variable i (using T::var)
        Var & var(size_t i) { return T::var(i); }

        /// Get const reference to factor I (using T::factor)
        const Factor & factor( size_t I ) const { return T::factor(I); }

        /// Get reference to factor I (using T::factor)
        Factor & factor( size_t I ) { return T::factor(I); }

        /// Factor I has been updated (using T::updatedFactor)
        void updatedFactor( size_t I ) { T::updatedFactor(I); }
};


typedef DAIAlg<FactorGraph> DAIAlgFG;
typedef DAIAlg<RegionGraph> DAIAlgRG;


/// Calculate the marginal of obj on ns by clamping 
/// all variables in ns and calculating logZ for each joined state
Factor calcMarginal( const InfAlg & obj, const VarSet & ns, bool reInit );


/// Calculate beliefs of all pairs in ns (by clamping
/// nodes in ns and calculating logZ and the beliefs for each state)
std::vector<Factor> calcPairBeliefs( const InfAlg & obj, const VarSet& ns, bool reInit );


/// Calculate beliefs of all pairs in ns (by clamping
/// pairs in ns and calculating logZ for each joined state)
std::vector<Factor> calcPairBeliefsNew( const InfAlg & obj, const VarSet& ns, bool reInit );


/// Calculate 2nd order interactions of the marginal of obj on ns
Factor calcMarginal2ndO( const InfAlg & obj, const VarSet& ns, bool reInit );


} // end of namespace dai


#endif
