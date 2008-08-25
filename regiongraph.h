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


#ifndef __REGIONGRAPH_H__
#define __REGIONGRAPH_H__


#include <iostream>
#include "bipgraph.h"
#include "factorgraph.h"
#include "weightedgraph.h"


namespace dai {


using namespace std;


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


typedef BipartiteGraph<FRegion,Region> BipRegGraph;


/// A RegionGraph is a bipartite graph consisting of outer regions (type FRegion) and inner regions (type Region)
class RegionGraph : public FactorGraph, BipRegGraph {
    public:
        typedef BipRegGraph::_nb_t      R_nb_t;
        typedef R_nb_t::const_iterator  R_nb_cit;
        typedef BipRegGraph::_edge_t    R_edge_t;

        
    protected:
        /// Give back the OR index that corresponds to a factor index
        typedef map<size_t,size_t>::const_iterator fac2OR_cit;
        map<size_t,size_t>              _fac2OR;


    public:
        /// Default constructor
        RegionGraph() : FactorGraph(), BipRegGraph(), _fac2OR() {}

        /// Constructs a RegionGraph from a FactorGraph
        RegionGraph(const FactorGraph & fg) : FactorGraph(fg), BipRegGraph(), _fac2OR() {}

        /// Constructs a RegionGraph from a FactorGraph, a vector of outer regions, a vector of inner regions and a vector of edges
        RegionGraph(const FactorGraph & fg, const vector<Region> & ors, const vector<Region> & irs, const vector<R_edge_t> & edges);
        
        /// Constructs a RegionGraph from a FactorGraph and a vector of outer VarSets (CVM style)
        RegionGraph(const FactorGraph & fg, const vector<VarSet> & cl);

        /// Copy constructor
        RegionGraph(const RegionGraph & x) : FactorGraph(x), BipRegGraph(x), _fac2OR(x._fac2OR) {}

        /// Assignment operator
        RegionGraph & operator=(const RegionGraph & x) {
            if( this != &x ) {
                FactorGraph::operator=(x);
                BipRegGraph::operator=(x);
                _fac2OR = x._fac2OR;
            }
            return *this;
        }

        /// Provides read access to outer region
        const FRegion & OR(long alpha) const { return BipRegGraph::V1(alpha); }
        /// Provides access to outer region
        FRegion & OR(long alpha) { return BipRegGraph::V1(alpha); }
        /// Provides read access to all outer regions
        const vector<FRegion> & ORs() const { return BipRegGraph::V1s(); }
        /// Provides access to all outer regions
        vector<FRegion> &ORs() { return BipRegGraph::V1s(); }
        /// Returns number of outer regions
        size_t nr_ORs() const { return BipRegGraph::V1s().size(); }

        /// Provides read access to inner region
        const Region & IR(long beta) const { return BipRegGraph::V2(beta); }
        /// Provides access to inner region
        Region & IR(long beta) { return BipRegGraph::V2(beta); }
        /// Provides read access to all inner regions
        const vector<Region> & IRs() const { return BipRegGraph::V2s(); }
        /// Provides access to all inner regions
        vector<Region> & IRs() { return BipRegGraph::V2s(); }
        /// Returns number of inner regions
        size_t nr_IRs() const { return BipRegGraph::V2s().size(); }

        /// Provides read access to edge
        const R_edge_t & Redge(size_t ind) const { return BipRegGraph::edge(ind); }
        /// Provides full access to edge
        R_edge_t & Redge(size_t ind) { return BipRegGraph::edge(ind); }
        /// Provides read access to all edges
        const std::vector<R_edge_t> & Redges() const { return BipRegGraph::edges(); }
        /// Provides full access to all edges
        std::vector<R_edge_t> & Redges() { return BipRegGraph::edges(); }
        /// Returns number of edges
        size_t nr_Redges() const { return BipRegGraph::edges().size(); }

        /// Provides read access to neighbours of outer region
        const R_nb_t & nbOR( size_t i1 ) const { return BipRegGraph::nb1(i1); }
        /// Provides full access to neighbours of outer region
        R_nb_t & nbOR( size_t i1 ) { return BipRegGraph::nb1(i1); }

        /// Provides read access to neighbours of inner region
        const R_nb_t & nbIR( size_t i2 ) const { return BipRegGraph::nb2(i2); }
        /// Provides full access to neighbours of inner region
        R_nb_t & nbIR( size_t i2 ) { return BipRegGraph::nb2(i2); }

        /// Converts the pair of outer/inner region indices (i1,i2) to the corresponding edge index
        size_t ORIR2E( const size_t i1, const size_t i2 ) const { return BipRegGraph::VV2E(i1, i2); }

        void Regenerate() { BipRegGraph::Regenerate(); }
        
        
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
        void makeCavity( const Var &n ) { FactorGraph::makeCavity( n ); RecomputeORs( n ); }

        /// We have to overload FactorGraph::makeFactorCavity because the corresponding outer regions have to be recomputed
        void makeFactorCavity( size_t I ) { FactorGraph::makeFactorCavity( I ); RecomputeOR( I ); }

        /// We have to overload FactorGraph::undoProbs because the corresponding outer regions have to be recomputed
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs( ns ); RecomputeORs( ns ); }

        /// We have to overload FactorGraph::undoProb because the corresponding outer regions have to be recomputed
        void undoProb( size_t I ) { FactorGraph::undoProb( I ); RecomputeOR( I ); }

        /// If updateFactor is called, we know that factor I has been changed and we should recompute the outer regions involving I
        void updatedFactor( size_t I ) { RecomputeOR( I ); }

        /// Send RegionGraph to output stream
        friend ostream & operator << ( ostream & os, const RegionGraph & rg );
};


}


#endif
