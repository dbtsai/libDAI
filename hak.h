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


#ifndef __HAK_H__
#define __HAK_H__


#include "daialg.h"
#include "regiongraph.h"
#include "enum.h"


namespace dai {


using namespace std;


/// HAK provides an implementation of the single and double-loop algorithms by Heskes, Albers and Kappen
class HAK : public DAIAlgRG {
    protected:
        vector<Factor>          _Qa;
        vector<Factor>          _Qb;
        vector<Factor>          _muab;
        vector<Factor>          _muba;
        
    public:
        /// Default constructor
        HAK() : DAIAlgRG() {};

        /// Copy constructor
        HAK(const HAK & x) : DAIAlgRG(x), _Qa(x._Qa), _Qb(x._Qb), _muab(x._muab), _muba(x._muba) {};

        /// Clone function
        HAK* clone() const { return new HAK(*this); }
        
        /// Construct from RegionGraph
        HAK(const RegionGraph & rg, const Properties &opts);

        /// Construct from RactorGraph using "clusters" option
        HAK(const FactorGraph & fg, const Properties &opts);

        /// Assignment operator
        HAK & operator=(const HAK & x) {
            if( this != &x ) {
                DAIAlgRG::operator=(x);
                _Qa         = x._Qa;
                _Qb         = x._Qb;
                _muab       = x._muab;
                _muba       = x._muba;
            }
            return *this;
        }
        
        static const char *Name;

        ENUM3(ClustersType,MIN,DELTA,LOOP)
        ClustersType Clusters() const { return GetPropertyAs<ClustersType>("clusters"); }
        bool DoubleLoop() { return GetPropertyAs<bool>("doubleloop"); }
        size_t LoopDepth() { return GetPropertyAs<size_t>("loopdepth"); }

        Factor & muab( size_t alpha, size_t beta ) { return _muab[ORIR2E(alpha,beta)]; }
        Factor & muba( size_t beta, size_t alpha ) { return _muba[ORIR2E(alpha,beta)]; }
        const Factor& Qa( size_t alpha ) const { return _Qa[alpha]; };
        const Factor& Qb( size_t beta ) const { return _Qb[beta]; };

//      void Regenerate();
        double doGBP();
        double doDoubleLoop();
        double run();
        void init();
        string identify() const;
        Factor belief( const Var &n ) const;
        Factor belief( const VarSet &ns ) const;
        vector<Factor> beliefs() const;
        Complex logZ () const;

        void init( const VarSet &ns );
        void undoProbs( const VarSet &ns ) { RegionGraph::undoProbs( ns ); init( ns ); }
        bool checkProperties();

    private:
        void constructMessages();
        void findLoopClusters( const FactorGraph &fg, set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars );
};


}


#endif
