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
        /// Clone *this (virtual copy constructor)
        virtual InfAlg* clone() const = 0;

        /// Create (virtual default constructor)
        virtual InfAlg* create() const = 0;
        
        /// Virtual desctructor (needed because this class contains virtual functions)
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
        virtual Real logZ() const = 0;

        /// Clear messages and beliefs
        virtual void init() = 0;

        /// Clear messages and beliefs corresponding to the nodes in ns
        virtual void init( const VarSet &ns ) = 0;

        /// The actual approximate inference algorithm
        virtual double run() = 0;

        /// Save factor I
        virtual void backupFactor( size_t I ) = 0;
        /// Save Factors involving ns
        virtual void backupFactors( const VarSet &ns ) = 0;

        /// Restore factor I
        virtual void restoreFactor( size_t I ) = 0;
        /// Restore Factors involving ns
        virtual void restoreFactors( const VarSet &ns ) = 0;

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$)
        virtual void clamp( const Var & n, size_t i, bool backup = false ) = 0;

        /// Set all factors interacting with var(i) to 1
        virtual void makeCavity( size_t i, bool backup = false ) = 0;

        /// Get reference to underlying FactorGraph
        virtual FactorGraph &fg() = 0;

        /// Get const reference to underlying FactorGraph
        virtual const FactorGraph &fg() const = 0;

        /// Return maximum difference between single node beliefs in the last pass
        virtual double maxDiff() const = 0;

        /// Return number of passes over the factorgraph
        virtual size_t Iterations() const = 0;
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

        /// Assignment operator
        DAIAlg & operator=( const DAIAlg &x ) {
            if( this != &x ) {
                InfAlg::operator=(x);
                T::operator=(x);
            }
            return *this;
        }

        /// Save factor I (using T::backupFactor)
        void backupFactor( size_t I ) { T::backupFactor( I ); }
        /// Save Factors involving ns (using T::backupFactors)
        void backupFactors( const VarSet &ns ) { T::backupFactors( ns ); }

        /// Restore factor I (using T::restoreFactor)
        void restoreFactor( size_t I ) { T::restoreFactor( I ); }
        /// Restore Factors involving ns (using T::restoreFactors)
        void restoreFactors( const VarSet &ns ) { T::restoreFactors( ns ); }

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$) (using T::clamp)
        void clamp( const Var & n, size_t i, bool backup = false ) { T::clamp( n, i, backup ); }

        /// Set all factors interacting with var(i) to 1 (using T::makeCavity)
        void makeCavity( size_t i, bool backup = false ) { T::makeCavity( i, backup ); }

        /// Get reference to underlying FactorGraph
        FactorGraph &fg() { return (FactorGraph &)(*this); }

        /// Get const reference to underlying FactorGraph
        const FactorGraph &fg() const { return (const FactorGraph &)(*this); }
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
