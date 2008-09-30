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
/// \brief Defines classes Region, FRegion and RegionGraph


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
        double          _c;

    public:
        /// Default constructor
        Region() : VarSet(), _c(1.0) {}

        /// Construct Region from a VarSet and a counting number
        Region(const VarSet & x, double c) : VarSet(x), _c(c) {}
        
        /// Copy constructor
        Region(const Region & x) : VarSet(x), _c(x._c) {}

        /// Assignment operator
        Region & operator=(const Region & x) {
            if( this != &x ) {
                VarSet::operator=(x);
                _c          = x._c;
            }
            return *this;
        }

        /// Provide read access to counting number
        const double & c() const { return _c; }
        /// Provide full access to counting number
        double & c() { return _c; }
};


/// A FRegion is a factor with a counting number
class FRegion : public Factor {
    private:
        /// Counting number
        double _c;

    public:
        /// Default constructor
        FRegion() : Factor(), _c(1.0) {}

        /// Constructs FRegion from a Factor and a counting number
        FRegion( const Factor & x, double c ) : Factor(x), _c(c) {}
        
        /// Copy constructor
        FRegion( const FRegion & x ) : Factor(x), _c(x._c) {}

        /// Assignment operator
        FRegion & operator=(const FRegion & x) {
            if( this != &x ) {
                Factor::operator=(x);
                _c = x._c;
            }
            return *this;
        }

        /// Provide read access to counting number
        const double & c() const { return _c; }
        /// Provide full access to counting number
        double & c() { return _c; }
};


/// A RegionGraph is a bipartite graph consisting of outer regions (type FRegion) and inner regions (type Region)
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
        /// Default constructor
        RegionGraph() : FactorGraph(), G(), ORs(), IRs(), fac2OR() {}

        /// Constructs a RegionGraph from a FactorGraph
        RegionGraph( const FactorGraph &fg ) : FactorGraph(fg), G(), ORs(), IRs(), fac2OR() {}

        /// Constructs a RegionGraph from a FactorGraph, a vector of outer regions, a vector of inner regions and a vector of edges
        RegionGraph( const FactorGraph &fg, const std::vector<Region> &ors, const std::vector<Region> &irs, const std::vector<std::pair<size_t,size_t> > &edges );
        
        /// Constructs a RegionGraph from a FactorGraph and a vector of outer VarSets (CVM style)
        RegionGraph( const FactorGraph &fg, const std::vector<VarSet> &cl );

        /// Copy constructor
        RegionGraph( const RegionGraph &x ) : FactorGraph(x), G(x.G), ORs(x.ORs), IRs(x.IRs), fac2OR(x.fac2OR) {}

        /// Assignment operator
        RegionGraph & operator=( const RegionGraph &x ) {
            if( this != &x ) {
                FactorGraph::operator=( x );
                G = x.G;
                ORs = x.ORs;
                IRs = x.IRs;
                fac2OR = x.fac2OR;
            }
            return *this;
        }

        /// Clone *this (virtual copy constructor)
        virtual RegionGraph* clone() const { return new RegionGraph(*this); }

        /// Create (virtual default constructor)
        virtual RegionGraph* create() const { return new RegionGraph(); }

        /// Set the content of the I'th factor and make a backup of its old content if backup == true
        virtual void setFactor( size_t I, const Factor &newFactor, bool backup = false ) {
            FactorGraph::setFactor( I, newFactor, backup ); 
            RecomputeOR( I ); 
        }

        /// Set the contents of all factors as specified by facs and make a backup of the old contents if backup == true
        virtual void setFactors( const std::map<size_t, Factor> & facs, bool backup = false ) {
            FactorGraph::setFactors( facs, backup );
            VarSet ns;
            for( std::map<size_t, Factor>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ )
                ns |= fac->second.vars();
            RecomputeORs( ns ); 
        }


        /// Provides read access to outer region
        const FRegion & OR(size_t alpha) const { return ORs[alpha]; }
        /// Provides access to outer region
        FRegion & OR(size_t alpha) { return ORs[alpha]; }

        /// Provides read access to inner region
        const Region & IR(size_t beta) const { return IRs[beta]; }
        /// Provides access to inner region
        Region & IR(size_t beta) { return IRs[beta]; }

        /// Returns number of outer regions
        size_t nrORs() const { return ORs.size(); }
        /// Returns number of inner regions
        size_t nrIRs() const { return IRs.size(); }


        /// Provides read access to neighbors of outer region
        const Neighbors & nbOR( size_t alpha ) const { return G.nb1(alpha); }
        /// Provides full access to neighbors of outer region
        Neighbors & nbOR( size_t alpha ) { return G.nb1(alpha); }

        /// Provides read access to neighbors of inner region
        const Neighbors & nbIR( size_t beta ) const { return G.nb2(beta); }
        /// Provides full access to neighbors of inner region
        Neighbors & nbIR( size_t beta ) { return G.nb2(beta); }


        /// Calculates counting numbers of inner regions based upon counting numbers of outer regions
        void Calc_Counting_Numbers();
        /// Check whether the counting numbers are valid
        bool Check_Counting_Numbers();

        /// Recompute all outer regions
        void RecomputeORs();

        /// Recompute all outer regions involving the variables in ns
        void RecomputeORs( const VarSet & ns );

        /// Recompute all outer regions involving factor I
        void RecomputeOR( size_t I );

        // Friends
        friend std::ostream & operator << ( std::ostream & os, const RegionGraph & rg );
};


} // end of namespace dai


#endif
