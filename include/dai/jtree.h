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


#ifndef __defined_libdai_jtree_h
#define __defined_libdai_jtree_h


#include <vector>
#include <string>
#include <dai/daialg.h>
#include <dai/varset.h>
#include <dai/regiongraph.h>
#include <dai/factorgraph.h>
#include <dai/clustergraph.h>
#include <dai/weightedgraph.h>
#include <dai/enum.h>
#include <dai/properties.h>


namespace dai {


class JTree : public DAIAlgRG {
    protected:
        DEdgeVec             _RTree;     // rooted tree
        std::vector<Factor>  _Qa;
        std::vector<Factor>  _Qb;
        std::vector<std::vector<Factor> >  _mes;
        double               _logZ;

    public:
        struct Properties {
            size_t verbose;
            DAI_ENUM(UpdateType,HUGIN,SHSH)
            UpdateType updates;
        } props;

    public:
        JTree() : DAIAlgRG(), _RTree(), _Qa(), _Qb(), _mes(), _logZ(), props() {}
        JTree( const JTree& x ) : DAIAlgRG(x), _RTree(x._RTree), _Qa(x._Qa), _Qb(x._Qb), _mes(x._mes), _logZ(x._logZ), props(x.props) {}
        JTree* clone() const { return new JTree(*this); }
        JTree & operator=( const JTree& x ) {
            if( this != &x ) {
                DAIAlgRG::operator=(x);
                _RTree  = x._RTree;
                _Qa     = x._Qa;
                _Qb     = x._Qb;
                _mes    = x._mes;
                _logZ   = x._logZ;
                props   = x.props;
            }
            return *this;
        }
        JTree( const FactorGraph &fg, const PropertySet &opts, bool automatic=true );
        void GenerateJT( const std::vector<VarSet> &Cliques );

        Factor & message( size_t alpha, size_t _beta ) { return _mes[alpha][_beta]; }   
        const Factor & message( size_t alpha, size_t _beta ) const { return _mes[alpha][_beta]; }   

        static const char *Name;
        std::string identify() const;
        void init() {}
        void runHUGIN();
        void runShaferShenoy();
        double run();
        Factor belief( const Var &n ) const;
        Factor belief( const VarSet &ns ) const;
        std::vector<Factor> beliefs() const;
        Real logZ() const;

        void init( const VarSet &/*ns*/ ) {}
        void restoreFactors( const VarSet &ns ) { RegionGraph::restoreFactors( ns ); init( ns ); }

        size_t findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot=(size_t)-1 ) const;
        Factor calcMarginal( const VarSet& ns );
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
        double maxDiff() const { return 0.0; }
};


} // end of namespace dai


#endif
