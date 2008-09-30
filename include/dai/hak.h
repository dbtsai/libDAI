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
/// \brief Defines class HAK.


#ifndef __defined_libdai_hak_h
#define __defined_libdai_hak_h


#include <string>
#include <dai/daialg.h>
#include <dai/regiongraph.h>
#include <dai/enum.h>
#include <dai/properties.h>


namespace dai {


/// Approximate inference algorithm: implementation of single-loop ("Generalized Belief Propagation") and double-loop algorithms by Heskes, Albers and Kappen
class HAK : public DAIAlgRG {
    private:
        std::vector<Factor>                _Qa;
        std::vector<Factor>                _Qb;
        std::vector<std::vector<Factor> >  _muab;
        std::vector<std::vector<Factor> >  _muba;
        /// Maximum difference encountered so far
        double _maxdiff;
        /// Number of iterations needed
        size_t _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible cluster choices
            DAI_ENUM(ClustersType,MIN,DELTA,LOOP)

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            double tol;

            /// Damping constant
            double damping;

            /// How to choose the clusters
            ClustersType clusters;

            /// Use single-loop (GBP) or double-loop (HAK)
            bool doubleloop;

            /// Depth of loops (only relevant for clusters == ClustersType::LOOP)
            size_t loopdepth;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        HAK() : DAIAlgRG(), _Qa(), _Qb(), _muab(), _muba(), _maxdiff(0.0), _iters(0U), props() {}

        /// Copy constructor
        HAK( const HAK &x ) : DAIAlgRG(x), _Qa(x._Qa), _Qb(x._Qb), _muab(x._muab), _muba(x._muba), _maxdiff(x._maxdiff), _iters(x._iters), props(x.props) {}

        /// Assignment operator
        HAK& operator=( const HAK &x ) {
            if( this != &x ) {
                DAIAlgRG::operator=( x );
                _Qa      = x._Qa;
                _Qb      = x._Qb;
                _muab    = x._muab;
                _muba    = x._muba;
                _maxdiff = x._maxdiff;
                _iters   = x._iters;
                props    = x.props;
            }
            return *this;
        }

        /// Construct from FactorGraph fg and PropertySet opts
        HAK( const FactorGraph &fg, const PropertySet &opts );

        /// Construct from RegionGraph rg and PropertySet opts
        HAK( const RegionGraph &rg, const PropertySet &opts );


        /// @name General InfAlg interface
        //@{
        virtual HAK* clone() const { return new HAK(*this); }
        virtual HAK* create() const { return new HAK(); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual double run();
        virtual double maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        //@}


        /// @name Additional interface specific for HAK
        //@{ 
        Factor & muab( size_t alpha, size_t _beta ) { return _muab[alpha][_beta]; }
        Factor & muba( size_t alpha, size_t _beta ) { return _muba[alpha][_beta]; }
        const Factor& Qa( size_t alpha ) const { return _Qa[alpha]; };
        const Factor& Qb( size_t beta ) const { return _Qb[beta]; };

        double doGBP();
        double doDoubleLoop();
        //@}

    private:
        void constructMessages();
        void findLoopClusters( const FactorGraph &fg, std::set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars );

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
