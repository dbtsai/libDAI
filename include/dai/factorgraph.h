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


class FactorGraph : public BipartiteGraph<Var,Factor> {
    protected:
        std::map<size_t,Prob>    _undoProbs;
        Prob::NormType           _normtype;

    public:
        /// Default constructor
        FactorGraph() : BipartiteGraph<Var,Factor>(), _undoProbs(), _normtype(Prob::NORMPROB) {};
        /// Copy constructor
        FactorGraph(const FactorGraph & x) : BipartiteGraph<Var,Factor>(x), _undoProbs(), _normtype(x._normtype) {};
        /// Construct FactorGraph from vector of Factors
        FactorGraph(const std::vector<Factor> &P);
        // Construct a FactorGraph from given factor and variable iterators
        template<typename FactorInputIterator, typename VarInputIterator>
        FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint = 0, size_t nr_var_hint = 0 );
        
        /// Assignment operator
        FactorGraph & operator=(const FactorGraph & x) {
            if(this!=&x) {
                BipartiteGraph<Var,Factor>::operator=(x);
                _undoProbs      = x._undoProbs;
                _normtype       = x._normtype;
            }
            return *this;
        }
        virtual ~FactorGraph() {}

        // aliases
        Var & var(size_t i) { return V1(i); }
        const Var & var(size_t i) const { return V1(i); }
        const std::vector<Var> & vars() const { return V1s(); }
        std::vector<Var> & vars() { return V1s(); }
        size_t nrVars() const { return V1s().size(); }
        Factor & factor(size_t I) { return V2(I); }
        const Factor & factor(size_t I) const { return V2(I); }
        const std::vector<Factor> & factors() const { return V2s(); }
        std::vector<Factor> & factors() { return V2s(); }
        size_t nrFactors() const { return V2s().size(); }

        /// Provides read access to neighbours of variable
        const _nb_t & nbV( size_t i1 ) const { return nb1(i1); }
        /// Provides full access to neighbours of variable
        _nb_t & nbV( size_t i1 ) { return nb1(i1); }
        /// Provides read access to neighbours of factor
        const _nb_t & nbF( size_t i2 ) const { return nb2(i2); }
        /// Provides full access to neighbours of factor
        _nb_t & nbF( size_t i2 ) { return nb2(i2); }

        size_t findVar(const Var & n) const {
            size_t i = find( vars().begin(), vars().end(), n ) - vars().begin();
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

        VarSet delta(const Var & n) const;
        VarSet Delta(const Var & n) const;
        virtual void makeFactorCavity(size_t I);
        virtual void makeCavity(const Var & n);

        long ReadFromFile(const char *filename);
        long WriteToFile(const char *filename) const;
        long WriteToDotFile(const char *filename) const;

        Factor ExactMarginal(const VarSet & x) const;
        Real ExactlogZ() const;

        virtual void clamp( const Var & n, size_t i );
        
        bool hasNegatives() const;
        Prob::NormType NormType() const { return _normtype; }
        
        std::vector<VarSet> Cliques() const;

        virtual void undoProbs( const VarSet &ns );
        void saveProbs( const VarSet &ns );
        virtual void undoProb( size_t I );
        void saveProb( size_t I );

        bool isConnected() const;

        virtual void updatedFactor( size_t /*I*/ ) {};

    private:
        /// Part of constructors (creates edges, neighbours and adjacency matrix)
        void createGraph( size_t nrEdges );
};


// assumes that the set of variables in [var_begin,var_end) is the union of the variables in the factors in [fact_begin, fact_end)
template<typename FactorInputIterator, typename VarInputIterator>
FactorGraph::FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint, size_t nr_var_hint ) : BipartiteGraph<Var,Factor>(), _undoProbs(), _normtype(Prob::NORMPROB) {
    // add factors
    size_t nrEdges = 0;
    V2s().reserve( nr_fact_hint );
    for( FactorInputIterator p2 = fact_begin; p2 != fact_end; ++p2 ) {
        V2s().push_back( *p2 );
	nrEdges += p2->vars().size();
    }
 
    // add variables
    V1s().reserve( nr_var_hint );
    for( VarInputIterator p1 = var_begin; p1 != var_end; ++p1 )
	V1s().push_back( *p1 );

    // create graph structure
    createGraph( nrEdges );
}


} // end of namespace dai


#endif
