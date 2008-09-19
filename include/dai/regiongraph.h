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


#ifndef __defined_libdai_regiongraph_h
#define __defined_libdai_regiongraph_h


#include <iostream>
#include <dai/bipgraph.h>
#include <dai/factorgraph.h>
#include <dai/weightedgraph.h>


namespace dai {


/// A Region is a set of variables with a counting number
class Region : public VarSet {
    protected:
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
    protected:
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
        BipartiteGraph          G;
        std::vector<FRegion>    ORs;
        std::vector<Region>     IRs;

    protected:
        /// Give back the OR index that corresponds to a factor index
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

        /// We have to overload FactorGraph::clamp because the corresponding outer regions have to be recomputed
        void clamp( const Var &n, size_t i ) { FactorGraph::clamp( n, i ); RecomputeORs( n ); }

        /// We have to overload FactorGraph::makeCavity because the corresponding outer regions have to be recomputed
        void makeCavity( size_t i ) { FactorGraph::makeCavity( i ); RecomputeORs( var(i) ); }

        /// We have to overload FactorGraph::undoProbs because the corresponding outer regions have to be recomputed
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs( ns ); RecomputeORs( ns ); }

        /// We have to overload FactorGraph::undoProb because the corresponding outer regions have to be recomputed
        void undoProb( size_t I ) { FactorGraph::undoProb( I ); RecomputeOR( I ); }

        /// Send RegionGraph to output stream
        friend std::ostream & operator << ( std::ostream & os, const RegionGraph & rg );
};


} // end of namespace dai


#endif
