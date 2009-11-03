/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines class TreeEP
/// \todo Improve documentation


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
        Real                  _maxdiff;
        /// Number of iterations needed
        size_t                _iters;

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
            Real tol;

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
                Real                 _logZ;


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
                Real logZ( const std::vector<Factor> &Qa, const std::vector<Factor> &Qb ) const;
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


    /// \name General InfAlg interface
    //@{
        virtual TreeEP* clone() const { return new TreeEP(*this); }
        virtual std::string identify() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    //@}


    /// \name Managing parameters (which are stored in TreeEP::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding TreeEP::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        void ConstructRG( const DEdgeVec &tree );
        bool offtree( size_t I ) const { return (fac2OR[I] == -1U); }
};


} // end of namespace dai


#endif
