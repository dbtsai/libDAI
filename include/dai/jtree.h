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


/// \file
/// \brief Defines class JTree


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


/// Exact inference algorithm using junction tree
class JTree : public DAIAlgRG {
    private:
        std::vector<std::vector<Factor> >  _mes;
        double               _logZ;

    public:
        /// Rooted tree
        DEdgeVec             RTree;
        
        /// Outer region beliefs
        std::vector<Factor>  Qa;
        
        /// Inner region beliefs
        std::vector<Factor>  Qb;

        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible JTree updates
            DAI_ENUM(UpdateType,HUGIN,SHSH)

            /// Verbosity
            size_t verbose;

            /// Type of updates
            UpdateType updates;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        JTree() : DAIAlgRG(), _mes(), _logZ(), RTree(), Qa(), Qb(), props() {}

        /// Copy constructor
        JTree( const JTree &x ) : DAIAlgRG(x), _mes(x._mes), _logZ(x._logZ), RTree(x.RTree), Qa(x.Qa), Qb(x.Qb), props(x.props) {}

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

        /// Construct from FactorGraph fg and PropertySet opts
        JTree( const FactorGraph &fg, const PropertySet &opts, bool automatic=true );


        /// @name General InfAlg interface
        //@{
        virtual JTree* clone() const { return new JTree(*this); }
        virtual JTree* create() const { return new JTree(); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init() {}
        virtual void init( const VarSet &/*ns*/ ) {}
        virtual double run();
        virtual double maxDiff() const { return 0.0; }
        virtual size_t Iterations() const { return 1UL; }
        //@}


        /// @name Additional interface specific for JTree
        //@{ 
        void GenerateJT( const std::vector<VarSet> &Cliques );

        /// Returns reference the message from outer region alpha to its _beta'th neighboring inner region
        Factor & message( size_t alpha, size_t _beta ) { return _mes[alpha][_beta]; }   
        /// Returns const reference to the message from outer region alpha to its _beta'th neighboring inner region
        const Factor & message( size_t alpha, size_t _beta ) const { return _mes[alpha][_beta]; }   

        /// Runs junction-tree with HUGIN updates
        void runHUGIN();

        /// Runs junction-tree with Shafer-Shenoy updates
        void runShaferShenoy();

        /// Finds an efficient tree for calculating the marginal of some variables
        size_t findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot=(size_t)-1 ) const;

        /// Calculates the marginal of a set of variables
        Factor calcMarginal( const VarSet& ns );
        //@}

    private:
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


/// Calculates upper bound to the treewidth of a FactorGraph
/** \relates JTree
 *  \return a pair (number of variables in largest clique, number of states in largest clique)
 */
std::pair<size_t,size_t> treewidth( const FactorGraph & fg );


} // end of namespace dai


#endif
