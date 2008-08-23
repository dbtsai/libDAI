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


#ifndef __TREEEP_H__
#define __TREEEP_H__


#include <vector>
#include "daialg.h"
#include "varset.h"
#include "regiongraph.h"
#include "factorgraph.h"
#include "clustergraph.h"
#include "weightedgraph.h"
#include "jtree.h"
#include "enum.h"


using namespace std;


class TreeEPSubTree {
    protected:
        vector<Factor>  _Qa;
        vector<Factor>  _Qb;
        DEdgeVec        _RTree;
        vector<size_t>  _a;             // _Qa[alpha]  <->  superTree._Qa[_a[alpha]]
        vector<size_t>  _b;             // _Qb[beta]   <->  superTree._Qb[_b[beta]]
                                        // _Qb[beta]   <->  _RTree[beta]    
        const Factor *  _I;
        VarSet          _ns;
        VarSet          _nsrem;
        double          _logZ;
        
        
    public:
        TreeEPSubTree() : _Qa(), _Qb(), _RTree(), _a(), _b(), _I(NULL), _ns(), _nsrem(), _logZ(0.0) {}
        TreeEPSubTree( const TreeEPSubTree &x) : _Qa(x._Qa), _Qb(x._Qb), _RTree(x._RTree), _a(x._a), _b(x._b), _I(x._I), _ns(x._ns), _nsrem(x._nsrem), _logZ(x._logZ) {}
        TreeEPSubTree & operator=( const TreeEPSubTree& x ) {
            if( this != &x ) {
                _Qa         = x._Qa;
                _Qb         = x._Qb;
                _RTree      = x._RTree;
                _a          = x._a;
                _b          = x._b;
                _I          = x._I;
                _ns         = x._ns;
                _nsrem      = x._nsrem;
                _logZ       = x._logZ;
            }
            return *this;
        }

        TreeEPSubTree( const DEdgeVec &subRTree, const DEdgeVec &jt_RTree, const vector<Factor> &jt_Qa, const vector<Factor> &jt_Qb, const Factor *I );
        void init();
        void InvertAndMultiply( const vector<Factor> &Qa, const vector<Factor> &Qb );
        void HUGIN_with_I( vector<Factor> &Qa, vector<Factor> &Qb );
        double logZ( const vector<Factor> &Qa, const vector<Factor> &Qb ) const;
        const Factor *& I() { return _I; }
};


class TreeEP : public JTree {
    protected:
        map<size_t, TreeEPSubTree>  _Q;

    public:
        ENUM2(TypeType,ORG,ALT)
        TypeType Type() const { return GetPropertyAs<TypeType>("type"); }
        
        TreeEP() : JTree(), _Q() {};
        TreeEP( const TreeEP& x ) : JTree(x), _Q(x._Q) {
            for( size_t I = 0; I < nrFactors(); I++ )
                if( offtree( I ) )
                    _Q[I].I() = &factor(I);
        }
        TreeEP* clone() const { return new TreeEP(*this); }
        TreeEP & operator=( const TreeEP& x ) {
            if( this != &x ) {
                JTree::operator=(x);
                _Q   = x._Q;
                for( size_t I = 0; I < nrFactors(); I++ )
                    if( offtree( I ) )
                        _Q[I].I() = &factor(I);
            }
            return *this;
        }
        TreeEP( const FactorGraph &fg, const Properties &opts );
        void ConstructRG( const DEdgeVec &tree );

        static const char *Name;
        string identify() const;
        void init();
        double run();
        Complex logZ() const;

        bool offtree(size_t I) const { return !_fac2OR.count(I); }

        void init( const VarSet &ns ) { init(); }
        void undoProbs( const VarSet &ns ) { RegionGraph::undoProbs( ns ); init( ns ); }
        bool checkProperties();
};


#endif
