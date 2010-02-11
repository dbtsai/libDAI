/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines classes Region, FRegion and RegionGraph, which implement a particular subclass of region graphs.


#ifndef __defined_libdai_regiongraph_h
#define __defined_libdai_regiongraph_h


#include <iostream>
#include <dai/bipgraph.h>
#include <dai/factorgraph.h>
#include <dai/weightedgraph.h>


namespace dai {


/// A Region is a set of variables with a counting number
class Region : public VarSet {
    private:
        /// Counting number
        Real          _c;

    public:
        /// Default constructor
        Region() : VarSet(), _c(1.0) {}

        /// Construct from a set of variables and a counting number
        Region( const VarSet &x, Real c ) : VarSet(x), _c(c) {}

        /// Returns constant reference to counting number
        const Real & c() const { return _c; }

        /// Returns reference to counting number
        Real & c() { return _c; }
};


/// An FRegion is a factor with a counting number
class FRegion : public Factor {
    private:
        /// Counting number
        Real _c;

    public:
        /// Default constructor
        FRegion() : Factor(), _c(1.0) {}

        /// Constructs from a factor and a counting number
        FRegion( const Factor & x, Real c ) : Factor(x), _c(c) {}

        /// Returns constant reference to counting number
        const Real & c() const { return _c; }

        /// Returns reference to counting number
        Real & c() { return _c; }
};


/// A RegionGraph combines a bipartite graph consisting of outer regions (type FRegion) and inner regions (type Region) with a FactorGraph
/** A RegionGraph inherits from a FactorGraph and adds additional structure in the form of a "region graph". Our definition of region graph
 *  is inspired by [\ref HAK03], which is less general than the definition given in [\ref YFW05].
 *
 *  The extra structure described by a RegionGraph over that described by a FactorGraph is:
 *  - a set of outer regions (indexed by \f$\alpha\f$), where each outer region consists of
 *    - a factor defined on a subset of variables
 *    - a counting number
 *  - a set of inner regions (indexed by \f$\beta\f$), where each inner region consists of
 *    - a subset of variables
 *    - a counting number
 *  - edges between inner and outer regions
 *
 *  Each factor in the factor graph belongs to an outer region; normally, the factor contents
 *  of an outer region would be the product of all the factors that belong to that region.
 */
class RegionGraph : public FactorGraph {
    public:
        /// Stores the neighborhood structure
        BipartiteGraph          G;

        /// The outer regions (corresponding to nodes of type 1)
        std::vector<FRegion>    ORs;

        /// The inner regions (corresponding to nodes of type 2)
        std::vector<Region>     IRs;

        /// Stores for each factor index the index of the outer region it belongs to
        std::vector<size_t>     fac2OR;


    public:
    /// \name Constructors and destructors
    //@{
        /// Default constructor
        RegionGraph() : FactorGraph(), G(), ORs(), IRs(), fac2OR() {}

        /// Partially constructs a region graph from a factor graph
        RegionGraph( const FactorGraph &fg ) : FactorGraph(fg), G(), ORs(), IRs(), fac2OR() {}

        /// Constructs a region graph from a factor graph, a vector of outer regions, a vector of inner regions and a vector of edges
        /** The counting numbers for the outer regions are set to 1.
         */
        RegionGraph( const FactorGraph &fg, const std::vector<VarSet> &ors, const std::vector<Region> &irs, const std::vector<std::pair<size_t,size_t> > &edges ) : FactorGraph(), G(), ORs(), IRs(), fac2OR() {
            construct( fg, ors, irs, edges );

            // Check counting numbers
#ifdef DAI_DEBUG
            checkCountingNumbers();
#endif
        }

        /// Constructs a region graph from a factor graph and a vector of outer clusters (CVM style)
        /** The region graph is constructed as in the Cluster Variation Method. 
         *
         *  The outer regions have as variable subsets the clusters specified in \a cl. 
         *  Each factor in the factor graph \a fg is assigned to one of the outer regions. 
         *  Each outer region gets counting number 1. 
         *
         *  The inner regions are (repeated) intersections of outer regions.
         *  An inner and an outer region are connected if the variables in the inner region form a
         *  subset of the variables in the outer region. The counting numbers for the inner
         *  regions are calculated by calcCountingNumbers() and satisfy the Moebius formula.
         */
        RegionGraph( const FactorGraph &fg, const std::vector<VarSet> &cl ) : FactorGraph(), G(), ORs(), IRs(), fac2OR() {
            constructCVM( fg, cl );

            // Check counting numbers
#ifdef DAI_DEBUG
            checkCountingNumbers();
#endif
        }

        /// Clone \c *this (virtual copy constructor)
        virtual RegionGraph* clone() const { return new RegionGraph(*this); }
    //@}

    /// \name Queries
    //@{
        /// Returns number of outer regions
        size_t nrORs() const { return ORs.size(); }
        /// Returns number of inner regions
        size_t nrIRs() const { return IRs.size(); }

        /// Returns constant reference to outer region \a alpha
        const FRegion & OR(size_t alpha) const { return ORs[alpha]; }
        /// Returns reference to outer region \a alpha
        FRegion & OR(size_t alpha) { return ORs[alpha]; }

        /// Returns constant reference to inner region \a beta
        const Region & IR(size_t beta) const { return IRs[beta]; }
        /// Returns reference to inner region \a beta
        Region & IR(size_t beta) { return IRs[beta]; }

        /// Returns constant reference to the neighbors of outer region \a alpha
        const Neighbors & nbOR( size_t alpha ) const { return G.nb1(alpha); }
        /// Returns constant reference to the neighbors of inner region \a beta
        const Neighbors & nbIR( size_t beta ) const { return G.nb2(beta); }

        /// Check whether the counting numbers are valid
        /** Counting numbers are said to be (variable) valid if for each variable \f$x\f$,
         *    \f[\sum_{\alpha \ni x} c_\alpha + \sum_{\beta \ni x} c_\beta = 1\f]
         *  or in words, if the sum of the counting numbers of the regions
         *  that contain the variable equals one.
         */
        bool checkCountingNumbers() const;
    //@}

    /// \name Operations
    //@{
        /// Set the content of the \a I 'th factor and make a backup of its old content if \a backup == \c true
        virtual void setFactor( size_t I, const Factor &newFactor, bool backup = false ) {
            FactorGraph::setFactor( I, newFactor, backup );
            RecomputeOR( I );
        }

        /// Set the contents of all factors as specified by \a facs and make a backup of the old contents if \a backup == \c true
        virtual void setFactors( const std::map<size_t, Factor> & facs, bool backup = false ) {
            FactorGraph::setFactors( facs, backup );
            VarSet ns;
            for( std::map<size_t, Factor>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ )
                ns |= fac->second.vars();
            RecomputeORs( ns );
        }

        /// Recompute all outer regions
        /** The factor contents of each outer region is set to the product of the factors belonging to that region.
         */
        void RecomputeORs();

        /// Recompute all outer regions involving the variables in \a vs
        /** The factor contents of each outer region involving at least one of the variables in \a vs is set to the product of the factors belonging to that region.
         */
        void RecomputeORs( const VarSet &vs );

        /// Recompute all outer regions involving factor \a I
        /** The factor contents of each outer region involving the \a I 'th factor is set to the product of the factors belonging to that region.
         */
        void RecomputeOR( size_t I );

        /// Calculates counting numbers of inner regions based upon counting numbers of outer regions
        /** The counting numbers of the inner regions are set using the Moebius inversion formula:
         *    \f[ c_\beta := 1 - \sum_{\gamma \in \mathrm{an}(\beta)} c_\gamma \f]
         *  where \f$\mathrm{an}(\beta)\f$ are the ancestors of inner region \f$\beta\f$ according to
         *  the partial ordering induced by the subset relation (i.e., a region is a child of another
         *  region if its variables are a subset of the variables of its parent region).
         */
        void calcCountingNumbers();
    //@}

    /// \name Input/output
    //@{
        /// Writes a RegionGraph to an output stream
        friend std::ostream & operator << ( std::ostream & os, const RegionGraph & rg );
    //@}

    protected:
        /// Helper function for constructors
        void construct( const FactorGraph &fg, const std::vector<VarSet> &ors, const std::vector<Region> &irs, const std::vector<std::pair<size_t,size_t> > &edges );

        /// Helper function for constructors (CVM style)
        void constructCVM( const FactorGraph &fg, const std::vector<VarSet> &cl );
};


} // end of namespace dai


#endif
