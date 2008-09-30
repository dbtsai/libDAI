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
/// \brief Defines class TreeEP


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


/// Approximate inference algorithm "TreeEP" by Minka and Qi
class TreeEP : public JTree {
    private:
        /// Maximum difference encountered so far
        double                  _maxdiff;
        /// Number of iterations needed
        size_t                  _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible choices for the tree
            DAI_ENUM(TypeType,ORG,ALT)

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            double tol;

            /// How to choose the tree
            TypeType type;
        } props; // FIXME: should be props2 because of conflict with JTree::props?

        /// Name of this inference method
        static const char *Name;

    private:
        class TreeEPSubTree {
            private:
                std::vector<Factor>  _Qa;
                std::vector<Factor>  _Qb;
                DEdgeVec             _RTree;
                std::vector<size_t>  _a;        // _Qa[alpha]  <->  superTree.Qa[_a[alpha]]
                std::vector<size_t>  _b;        // _Qb[beta]   <->  superTree.Qb[_b[beta]]
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

        std::map<size_t, TreeEPSubTree>  _Q;

    public:
        /// Default constructor
        TreeEP() : JTree(), _maxdiff(0.0), _iters(0), props(), _Q() {}

        /// Copy constructor
        TreeEP( const TreeEP &x ) : JTree(x), _maxdiff(x._maxdiff), _iters(x._iters), props(x.props), _Q(x._Q) {
            for( size_t I = 0; I < nrFactors(); I++ )
                if( offtree( I ) )
                    _Q[I].I() = &factor(I);
        }

        /// Assignment operator
        TreeEP& operator=( const TreeEP &x ) {
            if( this != &x ) {
                JTree::operator=( x );
                _maxdiff = x._maxdiff;
                _iters   = x._iters;
                props    = x.props;
                _Q       = x._Q;
                for( size_t I = 0; I < nrFactors(); I++ )
                    if( offtree( I ) )
                        _Q[I].I() = &factor(I);
            }
            return *this;
        }

        /// Construct from FactorGraph fg and PropertySet opts
        TreeEP( const FactorGraph &fg, const PropertySet &opts );


        /// @name General InfAlg interface
        //@{
        virtual TreeEP* clone() const { return new TreeEP(*this); }
        virtual TreeEP* create() const { return new TreeEP(); }
        virtual std::string identify() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual double run();
        virtual double maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        //@}


        /// @name Additional interface specific for TreeEP
        //@{ 
        //@}

    private:
        void ConstructRG( const DEdgeVec &tree );
        bool offtree( size_t I ) const { return (fac2OR[I] == -1U); }

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
