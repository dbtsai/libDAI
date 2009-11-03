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
/// \brief Defines class HAK.
/// \todo Improve documentation


#ifndef __defined_libdai_hak_h
#define __defined_libdai_hak_h


#include <string>
#include <dai/daialg.h>
#include <dai/regiongraph.h>
#include <dai/enum.h>
#include <dai/properties.h>


namespace dai {


/// Approximate inference algorithm: implementation of single-loop ("Generalized Belief Propagation") and double-loop algorithms by Heskes, Albers and Kappen
/** \todo Optimize HAK with precalculated indices, similarly to BP.
 */
class HAK : public DAIAlgRG {
    private:
        std::vector<Factor>                _Qa;
        std::vector<Factor>                _Qb;
        std::vector<std::vector<Factor> >  _muab;
        std::vector<std::vector<Factor> >  _muba;
        /// Maximum difference encountered so far
        Real _maxdiff;
        /// Number of iterations needed
        size_t _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible cluster choices
            DAI_ENUM(ClustersType,MIN,DELTA,LOOP);

            /// Enumeration of possible message initializations
            DAI_ENUM(InitType,UNIFORM,RANDOM);

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            Real tol;

            /// Damping constant
            Real damping;

            /// How to choose the clusters
            ClustersType clusters;

            /// How to initialize the messages
            InitType init;

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

        /// Construct from FactorGraph fg and PropertySet opts
        HAK( const FactorGraph &fg, const PropertySet &opts );

        /// Construct from RegionGraph rg and PropertySet opts
        HAK( const RegionGraph &rg, const PropertySet &opts );


    /// \name General InfAlg interface
    //@{
        virtual HAK* clone() const { return new HAK(*this); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    //@}


    /// \name Additional interface specific for HAK
    //@{
        Factor & muab( size_t alpha, size_t _beta ) { return _muab[alpha][_beta]; }
        Factor & muba( size_t alpha, size_t _beta ) { return _muba[alpha][_beta]; }
        const Factor& Qa( size_t alpha ) const { return _Qa[alpha]; };
        const Factor& Qb( size_t beta ) const { return _Qb[beta]; };

        Real doGBP();
        Real doDoubleLoop();
    //@}

    /// \name Managing parameters (which are stored in HAK::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding HAK::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        void constructMessages();
        void findLoopClusters( const FactorGraph &fg, std::set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars );
};


} // end of namespace dai


#endif
