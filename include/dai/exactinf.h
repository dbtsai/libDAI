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
/// \brief Defines ExactInf class
/// \todo Improve documentation


#ifndef __defined_libdai_exactinf_h
#define __defined_libdai_exactinf_h


#include <dai/daialg.h>
#include <dai/properties.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>


namespace dai {


/// Exact inference algorithm using brute force enumeration (mainly useful for testing purposes)
class ExactInf : public DAIAlgFG {
    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Verbosity
            size_t verbose;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    private:
        std::vector<Factor> _beliefsV;
        std::vector<Factor> _beliefsF;
        Real                _logZ;

    public:
        /// Default constructor
        ExactInf() : DAIAlgFG(), props(), _beliefsV(), _beliefsF(), _logZ(0) {}

        /// Construct from FactorGraph fg and PropertySet opts
        ExactInf( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), props(), _beliefsV(), _beliefsF(), _logZ() {
            setProperties( opts );
            construct();
        }


        /// @name General InfAlg interface
        //@{
        virtual ExactInf* clone() const { return new ExactInf(*this); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const { return beliefV( findVar( n ) ); }
        virtual Factor belief( const VarSet &ns ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { return _logZ; }
        virtual void init();
        virtual void init( const VarSet &/*ns*/ ) { DAI_THROW(NOT_IMPLEMENTED); }
        virtual double run();
        virtual double maxDiff() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual size_t Iterations() const { DAI_THROW(NOT_IMPLEMENTED); return 0; }
        //@}


        /// @name Additional interface specific for ExactInf
        //@{
        Factor beliefV( size_t i ) const { return _beliefsV[i]; }
        Factor beliefF( size_t I ) const { return _beliefsF[I]; }
        //@}

    private:
        void construct();
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
