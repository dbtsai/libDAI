/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

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
    private:
        std::vector<std::vector<Factor> >  _mes;
        double               _logZ;

    public:
        DEdgeVec             RTree;     // rooted tree
        std::vector<Factor>  Qa;
        std::vector<Factor>  Qb;
        struct Properties {
            size_t verbose;
            DAI_ENUM(UpdateType,HUGIN,SHSH)
            UpdateType updates;
        } props;
        /// Name of this inference method
        static const char *Name;

    public:
        /// Default constructor
        JTree() : DAIAlgRG(), _mes(), _logZ(), RTree(), Qa(), Qb(), props() {}

        /// Construct from FactorGraph fg and PropertySet opts
        JTree( const FactorGraph &fg, const PropertySet &opts, bool automatic=true );

        /// Copy constructor
        JTree( const JTree &x ) : DAIAlgRG(x), _mes(x._mes), _logZ(x._logZ), RTree(x.RTree), Qa(x.Qa), Qb(x.Qb), props(x.props) {}

        /// Clone *this (virtual copy constructor)
        virtual JTree* clone() const { return new JTree(*this); }

        /// Create (virtual default constructor)
        virtual JTree* create() const { return new JTree(); }

        /// Assignment operator
        JTree& operator=( const JTree &x ) {
            if( this != &x ) {
                DAIAlgRG::operator=( x );
                _mes    = x._mes;
                _logZ   = x._logZ;
                RTree   = x.RTree;
                Qa      = x.Qa;
                Qb      = x.Qb;
                props   = x.props;
            }
            return *this;
        }

        /// Identifies itself for logging purposes
        virtual std::string identify() const;

        /// Get single node belief
        virtual Factor belief( const Var &n ) const;

        /// Get general belief
        virtual Factor belief( const VarSet &ns ) const;

        /// Get all beliefs
        virtual std::vector<Factor> beliefs() const;

        /// Get log partition sum
        virtual Real logZ() const;

        /// Clear messages and beliefs
        virtual void init() {}

        /// Clear messages and beliefs corresponding to the nodes in ns
        virtual void init( const VarSet &/*ns*/ ) {}

        /// The actual approximate inference algorithm
        virtual double run();

        /// Return maximum difference between single node beliefs in the last pass
        virtual double maxDiff() const { return 0.0; }

        /// Return number of passes over the factorgraph
        virtual size_t Iterations() const { return 1UL; }


        void GenerateJT( const std::vector<VarSet> &Cliques );

        Factor & message( size_t alpha, size_t _beta ) { return _mes[alpha][_beta]; }   
        const Factor & message( size_t alpha, size_t _beta ) const { return _mes[alpha][_beta]; }   

        void runHUGIN();
        void runShaferShenoy();
        size_t findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot=(size_t)-1 ) const;
        Factor calcMarginal( const VarSet& ns );

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


std::pair<size_t,size_t> treewidth( const FactorGraph & fg );


} // end of namespace dai


#endif
