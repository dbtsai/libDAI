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


#ifndef __JTREE_H__
#define __JTREE_H__


#include <vector>
#include "daialg.h"
#include "varset.h"
#include "regiongraph.h"
#include "factorgraph.h"
#include "clustergraph.h"
#include "weightedgraph.h"
#include "enum.h"


namespace dai {


using namespace std;


class JTree : public DAIAlgRG {
    protected:
        DEdgeVec        _RTree;     // rooted tree
        vector<Factor>  _Qa;
        vector<Factor>  _Qb;
        vector<Factor>  _mes;
        double          _logZ;


    public:
        ENUM2(UpdateType,HUGIN,SHSH)
        UpdateType Updates() const { return GetPropertyAs<UpdateType>("updates"); }

        JTree() : DAIAlgRG(), _RTree(), _Qa(), _Qb(), _mes(), _logZ() {};
        JTree( const JTree& x ) : DAIAlgRG(x), _RTree(x._RTree), _Qa(x._Qa), _Qb(x._Qb), _mes(x._mes), _logZ(x._logZ) {};
        JTree* clone() const { return new JTree(*this); }
        JTree & operator=( const JTree& x ) {
            if( this != &x ) {
                DAIAlgRG::operator=(x);
                _RTree  = x._RTree;
                _Qa     = x._Qa;
                _Qb     = x._Qb;
                _mes    = x._mes;
                _logZ   = x._logZ;
            }
            return *this;
        }
        JTree( const FactorGraph &fg, const Properties &opts, bool automatic=true );
        void GenerateJT( const vector<VarSet> &Cliques );

        Factor & message(size_t i1, size_t i2) { return( _mes[ORIR2E(i1,i2)] ); }   
        const Factor & message(size_t i1, size_t i2) const { return( _mes[ORIR2E(i1,i2)] ); }   

        static const char *Name;
        string identify() const;
//      void Regenerate();
        void init() {
            assert( checkProperties() );
        }
        void runHUGIN();
        void runShaferShenoy();
        double run();
        Factor belief( const Var &n ) const;
        Factor belief( const VarSet &ns ) const;
        vector<Factor> beliefs() const;
        Complex logZ() const;

        void init( const VarSet &ns ) {}
        void undoProbs( const VarSet &ns ) { RegionGraph::undoProbs( ns ); init( ns ); }

        size_t findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot=(size_t)-1 ) const;
        Factor calcMarginal( const VarSet& ns );
        bool checkProperties();
};


}


#endif
