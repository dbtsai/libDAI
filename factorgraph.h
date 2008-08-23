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


#ifndef __FACTORGRAPH_H__
#define __FACTORGRAPH_H__


#include <iostream>
#include <map>
#include "bipgraph.h"
#include "factor.h"


using namespace std;


class FactorGraph : public BipartiteGraph<Var,Factor> {
    protected:
        map<size_t,Prob>    _undoProbs;
        bool                _hasNegatives;
        Prob::NormType      _normtype;

    public:
        /// Default constructor
        FactorGraph() : BipartiteGraph<Var,Factor>(), _undoProbs(), _hasNegatives(false), _normtype(Prob::NORMPROB) {};
        /// Copy constructor
        FactorGraph(const FactorGraph & x) : BipartiteGraph<Var,Factor>(x), _undoProbs(), _hasNegatives(x._hasNegatives), _normtype(x._normtype) {};
        /// Construct FactorGraph from vector of Factors
        FactorGraph(const vector<Factor> &P);
        /// Assignment operator
        FactorGraph & operator=(const FactorGraph & x) {
            if(this!=&x) {
                BipartiteGraph<Var,Factor>::operator=(x);
                _undoProbs      = x._undoProbs;
                _hasNegatives   = x._hasNegatives;
                _normtype       = x._normtype;
            }
            return *this;
        }
        virtual ~FactorGraph() {}

        // aliases
        Var & var(size_t i) { return V1(i); }
        const Var & var(size_t i) const { return V1(i); }
        const vector<Var> & vars() const { return V1s(); }
        vector<Var> & vars() { return V1s(); }
        size_t nrVars() const { return V1s().size(); }
        Factor & factor(size_t I) { return V2(I); }
        const Factor & factor(size_t I) const { return V2(I); }
        const vector<Factor> & factors() const { return V2s(); }
        vector<Factor> & factors() { return V2s(); }
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

        friend ostream& operator << (ostream& os, const FactorGraph& fg);
        friend istream& operator >> (istream& is, FactorGraph& fg);

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
        
        bool hasNegatives() const { return _hasNegatives; }
        Prob::NormType NormType() const { return _normtype; }
        
        vector<VarSet> Cliques() const;

        virtual void undoProbs( const VarSet &ns );
        void saveProbs( const VarSet &ns );
        virtual void undoProb( size_t I );
        void saveProb( size_t I );

        bool isConnected() const;

        virtual void updatedFactor( size_t I ) {};
};


bool hasShortLoops(const vector<Factor> &P);
void RemoveShortLoops(vector<Factor> &P);


#endif
