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


#ifndef __defined_libdai_factorgraph_h
#define __defined_libdai_factorgraph_h


#include <iostream>
#include <map>
#include <tr1/unordered_map>
#include <dai/bipgraph.h>
#include <dai/factor.h>


namespace dai {


bool hasShortLoops( const std::vector<Factor> &P );
void RemoveShortLoops( std::vector<Factor> &P );


class FactorGraph {
    public:
        BipartiteGraph         G;
        std::vector<Var>       vars;
        std::vector<Factor>    factors;
        typedef BipartiteGraph::Neighbor  Neighbor;
        typedef BipartiteGraph::Neighbors Neighbors;
        typedef BipartiteGraph::Edge      Edge;

    protected:
        std::map<size_t,Prob>  _undoProbs;
        Prob::NormType         _normtype;

    public:
        /// Default constructor
        FactorGraph() : G(), vars(), factors(), _undoProbs(), _normtype(Prob::NORMPROB) {};
        /// Copy constructor
        FactorGraph(const FactorGraph & x) : G(x.G), vars(x.vars), factors(x.factors), _undoProbs(x._undoProbs), _normtype(x._normtype) {};
        /// Construct FactorGraph from vector of Factors
        FactorGraph(const std::vector<Factor> &P);
        // Construct a FactorGraph from given factor and variable iterators
        template<typename FactorInputIterator, typename VarInputIterator>
        FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint = 0, size_t nr_var_hint = 0 );
        
        /// Assignment operator
        FactorGraph & operator=(const FactorGraph & x) {
            if( this != &x ) {
                G          = x.G;
                vars       = x.vars;
                factors    = x.factors;
                _undoProbs = x._undoProbs;
                _normtype  = x._normtype;
            }
            return *this;
        }
        virtual ~FactorGraph() {}

        // aliases
        Var & var(size_t i) { return vars[i]; }
        const Var & var(size_t i) const { return vars[i]; }
        Factor & factor(size_t I) { return factors[I]; }
        const Factor & factor(size_t I) const { return factors[I]; }

        size_t nrVars() const { return vars.size(); }
        size_t nrFactors() const { return factors.size(); }
        size_t nrEdges() const { return G.nrEdges(); }

        /// Provides read access to neighbors of variable
        const Neighbors & nbV( size_t i ) const { return G.nb1(i); }
        /// Provides full access to neighbors of variable
        Neighbors & nbV( size_t i ) { return G.nb1(i); }
        /// Provides read access to neighbors of factor
        const Neighbors & nbF( size_t I ) const { return G.nb2(I); }
        /// Provides full access to neighbors of factor
        Neighbors & nbF( size_t I ) { return G.nb2(I); }
        /// Provides read access to neighbor of variable
        const Neighbor & nbV( size_t i, size_t _I ) const { return G.nb1(i)[_I]; }
        /// Provides full access to neighbor of variable
        Neighbor & nbV( size_t i, size_t _I ) { return G.nb1(i)[_I]; }
        /// Provides read access to neighbor of factor
        const Neighbor & nbF( size_t I, size_t _i ) const { return G.nb2(I)[_i]; }
        /// Provides full access to neighbor of factor
        Neighbor & nbF( size_t I, size_t _i ) { return G.nb2(I)[_i]; }

        size_t findVar(const Var & n) const {
            size_t i = find( vars.begin(), vars.end(), n ) - vars.begin();
            assert( i != nrVars() );
            return i;
        }
        size_t findFactor(const VarSet &ns) const {
            size_t I;
            for( I = 0; I < nrFactors(); I++ )
                if( factor(I).vars() == ns )
                    break;
            assert( I != nrFactors() );
            return I;
        }

        friend std::ostream& operator << (std::ostream& os, const FactorGraph& fg);
        friend std::istream& operator >> (std::istream& is, FactorGraph& fg);

        VarSet delta( unsigned i ) const;
        VarSet Delta( unsigned i ) const;
        virtual void makeCavity( unsigned i );

        long ReadFromFile(const char *filename);
        long WriteToFile(const char *filename) const;
        long WriteToDotFile(const char *filename) const;

        virtual void clamp( const Var & n, size_t i );
        
        bool hasNegatives() const;
        Prob::NormType NormType() const { return _normtype; }
        
        std::vector<VarSet> Cliques() const;

        virtual void undoProbs( const VarSet &ns );
        void saveProbs( const VarSet &ns );
        virtual void undoProb( size_t I );
        void saveProb( size_t I );

        virtual void updatedFactor( size_t /*I*/ ) {};

    private:
        /// Part of constructors (creates edges, neighbors and adjacency matrix)
        void createGraph( size_t nrEdges );
};


// assumes that the set of variables in [var_begin,var_end) is the union of the variables in the factors in [fact_begin, fact_end)
template<typename FactorInputIterator, typename VarInputIterator>
FactorGraph::FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint, size_t nr_var_hint ) : G(), _undoProbs(), _normtype(Prob::NORMPROB) {
    // add factors
    size_t nrEdges = 0;
    factors.reserve( nr_fact_hint );
    for( FactorInputIterator p2 = fact_begin; p2 != fact_end; ++p2 ) {
        factors.push_back( *p2 );
        nrEdges += p2->vars().size();
    }
 
    // add variables
    vars.reserve( nr_var_hint );
    for( VarInputIterator p1 = var_begin; p1 != var_end; ++p1 )
        vars.push_back( *p1 );

    // create graph structure
    createGraph( nrEdges );
}


} // end of namespace dai


#endif
