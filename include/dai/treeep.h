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


#ifndef __defined_libdai_treeep_h
#define __defined_libdai_treeep_h


#include <vector>
#include <string>
#include <dai/daialg.h>
#include <dai/varset.h>
#include <dai/regiongraph.h>
#include <dai/factorgraph.h>
#include <dai/clustergraph.h>
#include <dai/weightedgraph.h>
#include <dai/jtree.h>
#include <dai/properties.h>
#include <dai/enum.h>


namespace dai {


class TreeEPSubTree {
    protected:
        std::vector<Factor>  _Qa;
        std::vector<Factor>  _Qb;
        DEdgeVec             _RTree;
        std::vector<size_t>  _a;        // _Qa[alpha]  <->  superTree._Qa[_a[alpha]]
        std::vector<size_t>  _b;        // _Qb[beta]   <->  superTree._Qb[_b[beta]]
                                        // _Qb[beta]   <->  _RTree[beta]    
        const Factor *       _I;
        VarSet               _ns;
        VarSet               _nsrem;
        double               _logZ;
        
        
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

        TreeEPSubTree( const DEdgeVec &subRTree, const DEdgeVec &jt_RTree, const std::vector<Factor> &jt_Qa, const std::vector<Factor> &jt_Qb, const Factor *I );
        void init();
        void InvertAndMultiply( const std::vector<Factor> &Qa, const std::vector<Factor> &Qb );
        void HUGIN_with_I( std::vector<Factor> &Qa, std::vector<Factor> &Qb );
        double logZ( const std::vector<Factor> &Qa, const std::vector<Factor> &Qb ) const;
        const Factor *& I() { return _I; }
};


class TreeEP : public JTree {
    protected:
        std::map<size_t, TreeEPSubTree>  _Q;

    public:
        struct Properties {
            size_t verbose;
            size_t maxiter;
            double tol;
            ENUM2(TypeType,ORG,ALT)
            TypeType type;
        } props;
        double maxdiff;

    public:
        /// Default constructor
        TreeEP() : JTree(), _Q(), props(), maxdiff() {};
        /// Copy constructor
        TreeEP( const TreeEP& x ) : JTree(x), _Q(x._Q), props(x.props), maxdiff(x.maxdiff) {
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
                props = x.props;
                maxdiff = x.maxdiff;
            }
            return *this;
        }
        TreeEP( const FactorGraph &fg, const PropertySet &opts );
        void ConstructRG( const DEdgeVec &tree );

        static const char *Name;
        std::string identify() const;
        void init();
        double run();
        Real logZ() const;

        bool offtree( size_t I ) const { return (fac2OR[I] == -1U); }

        void init( const VarSet &/*ns*/ ) { init(); }
        void undoProbs( const VarSet &ns ) { RegionGraph::undoProbs( ns ); init( ns ); }

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        double maxDiff() const { return maxdiff; }
};


} // end of namespace dai


#endif
