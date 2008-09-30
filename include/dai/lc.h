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
/// \brief Defines class LC


#ifndef __defined_libdai_lc_h
#define __defined_libdai_lc_h


#include <string>
#include <dai/daialg.h>
#include <dai/enum.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/exceptions.h>


namespace dai {


/// Approximate inference algorithm "Loop Corrected Belief Propagation" by Mooij and Kappen
class LC : public DAIAlgFG {
    private:
        std::vector<Factor>      _pancakes;      // used by all LC types (psi_I is stored in the pancake)
        std::vector<Factor>      _cavitydists;   // used by all LC types to store the approximate cavity distribution
        /// _phis[i][_I] corresponds to \f$ \phi^{\setminus i}_I(x_{I \setminus i}) \f$
        std::vector<std::vector<Factor> >      _phis;

        /// Single variable beliefs
        std::vector<Factor>      _beliefs;

        /// Maximum difference encountered so far
        double                  _maxdiff;
        /// Number of iterations needed
        size_t                  _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible ways to initialize the cavities
            DAI_ENUM(CavityType,FULL,PAIR,PAIR2,UNIFORM)

            /// Enumeration of different update schedules
            DAI_ENUM(UpdateType,SEQFIX,SEQRND,NONE)

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            double tol;

            /// Complete or partial reinit of cavity graphs?
            bool reinit;

            /// Damping constant
            double damping;

            /// How to initialize the cavities
            CavityType cavity;

            /// What update schedule to use
            UpdateType updates;

            /// Name of the algorithm used to initialize the cavity distributions
            std::string cavainame;      // FIXME: needs assignment operator?

            /// Parameters for the algorithm used to initialize the cavity distributions
            PropertySet cavaiopts;      // FIXME: needs assignment operator?
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        LC() : DAIAlgFG(), _pancakes(), _cavitydists(), _phis(), _beliefs(), _maxdiff(), _iters(), props() {}

        /// Copy constructor
        LC( const LC &x ) : DAIAlgFG(x), _pancakes(x._pancakes), _cavitydists(x._cavitydists), _phis(x._phis), _beliefs(x._beliefs), _maxdiff(x._maxdiff), _iters(x._iters), props(x.props) {}

        /// Assignment operator
        LC& operator=( const LC &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _pancakes     = x._pancakes;
                _cavitydists  = x._cavitydists;
                _phis         = x._phis;
                _beliefs      = x._beliefs;
                _maxdiff      = x._maxdiff;
                _iters        = x._iters;
                props         = x.props;
            }
            return *this;
        }

        /// Construct from FactorGraph fg and PropertySet opts
        LC( const FactorGraph &fg, const PropertySet &opts );


        /// @name General InfAlg interface
        //@{
        virtual LC* clone() const { return new LC(*this); }
        virtual LC* create() const { return new LC(); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const { return( _beliefs[findVar(n)] ); }
        virtual Factor belief( const VarSet &/*ns*/ ) const { DAI_THROW(NOT_IMPLEMENTED); return Factor(); }
        virtual std::vector<Factor> beliefs() const { return _beliefs; }
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual double run();
        virtual double maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        //@}


        /// @name Additional interface specific for LC
        //@{ 
        double CalcCavityDist( size_t i, const std::string &name, const PropertySet &opts );
        double InitCavityDists( const std::string &name, const PropertySet &opts );
        long SetCavityDists( std::vector<Factor> &Q );

        Factor NewPancake (size_t i, size_t _I, bool & hasNaNs);

        void CalcBelief (size_t i);
        const Factor &belief (size_t i) const { return _beliefs[i]; };
        const Factor &pancake (size_t i) const { return _pancakes[i]; };
        const Factor &cavitydist (size_t i) const { return _cavitydists[i]; };
        //@}

    private:
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
