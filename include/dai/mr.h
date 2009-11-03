/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2007       Bastian Wemmenhove
 *  Copyright (C) 2007-2009  Joris Mooij         [joris dot mooij at libdai dot org]
 *  Copyright (C) 2007       Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines class MR
/// \todo Improve documentation


#ifndef __defined_libdai_mr_h
#define __defined_libdai_mr_h


#include <vector>
#include <string>
#include <dai/factorgraph.h>
#include <dai/daialg.h>
#include <dai/enum.h>
#include <dai/properties.h>
#include <dai/exceptions.h>
#include <boost/dynamic_bitset.hpp>


namespace dai {


/// Approximate inference algorithm by Montanari and Rizzo
class MR : public DAIAlgFG {
    private:
        bool supported;                                            // is the underlying factor graph supported?

        std::vector<size_t>                             con;       // con[i] = connectivity of spin i
        std::vector<std::vector<size_t> >               nb;        // nb[i] are the neighbours of spin i
        std::vector<std::vector<Real> >                 tJ;        // tJ[i][_j] is the tanh of the interaction between spin i and its neighbour nb[i][_j]
        std::vector<Real>                               theta;     // theta[i] is the local field on spin i
        std::vector<std::vector<Real> >                 M;         // M[i][_j] is M^{(i)}_j
        std::vector<std::vector<size_t> >               kindex;    // the _j'th neighbour of spin i has spin i as its kindex[i][_j]'th neighbour
        std::vector<std::vector<std::vector<Real> > >   cors;

        static const size_t kmax = 31;
        typedef boost::dynamic_bitset<> sub_nb;

        size_t N;

        std::vector<Real> Mag;

        Real _maxdiff;
        size_t _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of different types of update equations
            DAI_ENUM(UpdateType,FULL,LINEAR);

            /// Enumeration of different ways of initializing the cavity correlations
            DAI_ENUM(InitType,RESPPROP,CLAMPING,EXACT);

            /// Verbosity
            size_t verbose;

            /// Tolerance
            Real tol;

            /// Update equations
            UpdateType updates;

            /// How to initialize the cavity correlations
            InitType inits;
        } props;

        /// Name of this inference method
        static const char *Name;

    public:
        /// Default constructor
        MR() : DAIAlgFG(), supported(), con(), nb(), tJ(), theta(), M(), kindex(), cors(), N(), Mag(), _maxdiff(), _iters(), props() {}

        /// Construct from FactorGraph fg and PropertySet opts
        MR( const FactorGraph &fg, const PropertySet &opts );


    /// \name General InfAlg interface
    //@{
        virtual MR* clone() const { return new MR(*this); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &/*ns*/ ) const { DAI_THROW(NOT_IMPLEMENTED); return Factor(); }
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init() {}
        virtual void init( const VarSet &/*ns*/ ) { DAI_THROW(NOT_IMPLEMENTED); }
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    //@}

    /// \name Managing parameters (which are stored in MR::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding MR::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        void init(size_t Nin, Real *_w, Real *_th);
        void makekindex();
        void init_cor();
        Real init_cor_resp();
        void solvemcav();
        void solveM();

        Real _tJ(size_t i, sub_nb A);

        Real Omega(size_t i, size_t _j, size_t _l);
        Real T(size_t i, sub_nb A);
        Real T(size_t i, size_t _j);
        Real Gamma(size_t i, size_t _j, size_t _l1, size_t _l2);
        Real Gamma(size_t i, size_t _l1, size_t _l2);

        Real appM(size_t i, sub_nb A);
        void sum_subs(size_t j, sub_nb A, Real *sum_even, Real *sum_odd);

        Real sign(Real a) { return (a >= 0) ? 1.0 : -1.0; }
};


} // end of namespace dai


#endif
