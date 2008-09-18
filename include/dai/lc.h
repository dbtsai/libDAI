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


#ifndef __defined_libdai_lc_h
#define __defined_libdai_lc_h


#include <string>
#include <dai/daialg.h>
#include <dai/enum.h>
#include <dai/factorgraph.h>


namespace dai {


class LC : public DAIAlgFG {
    protected:
        std::vector<Factor>      _pancakes;      // used by all LC types (psi_I is stored in the pancake)
        std::vector<Factor>      _cavitydists;   // used by all LC types to store the approximate cavity distribution
        /// _phis[i][_I] corresponds to \f$ \phi^{\setminus i}_I(x_{I \setminus i}) \f$
        std::vector<std::vector<Factor> >      _phis;

        /// Single variable beliefs
        std::vector<Factor>      _beliefs;

    public:
        ENUM6(CavityType,FULL,PAIR,PAIR2,PAIRINT,PAIRCUM,UNIFORM)
        ENUM3(UpdateType,SEQFIX,SEQRND,NONE)

        CavityType Cavity() const { return GetPropertyAs<CavityType>("cavity"); }
        UpdateType Updates() const { return GetPropertyAs<UpdateType>("updates"); }
        bool reInit() const { return GetPropertyAs<bool>("reinit"); }
        
        /// Default constructor
        LC() : DAIAlgFG() {};
        /// Copy constructor
        LC(const LC & x) : DAIAlgFG(x), _pancakes(x._pancakes), _cavitydists(x._cavitydists), _phis(x._phis), _beliefs(x._beliefs) {};
        /// Clone function
        LC* clone() const { return new LC(*this); }
        /// Construct LC object from a FactorGraph and parameters
        LC(const FactorGraph & fg, const Properties &opts);
        /// Assignment operator
        LC& operator=(const LC & x) {
            if( this != &x ) {
                DAIAlgFG::operator=(x);
                _pancakes       = x._pancakes;
                _cavitydists    = x._cavitydists;
                _phis           = x._phis;
                _beliefs        = x._beliefs;
            }
            return *this;
        }

        static const char *Name;
        double CalcCavityDist( size_t i, const std::string &name, const Properties &opts );
        double InitCavityDists( const std::string &name, const Properties &opts );
        long SetCavityDists( std::vector<Factor> &Q );

        void init();
        Factor NewPancake (size_t i, size_t _I, bool & hasNaNs);
        double run();

        std::string identify() const;
        Factor belief (const Var &n) const { return( _beliefs[findVar(n)] ); }
        Factor belief (const VarSet &/*ns*/) const { assert( 0 == 1 ); return Factor(); }
        std::vector<Factor> beliefs() const { return _beliefs; }
        Complex logZ() const { /*assert( 0 == 1 );*/ return 0.0; }
        void CalcBelief (size_t i);
        const Factor & belief (size_t i) const { return _beliefs[i]; };
        const Factor & pancake (size_t i) const { return _pancakes[i]; };
        const Factor & cavitydist (size_t i) const { return _cavitydists[i]; };

        void clamp( const Var &/*n*/, size_t /*i*/ ) { assert( 0 == 1 ); }
        void undoProbs( const VarSet &/*ns*/ ) { assert( 0 == 1 ); }
        void saveProbs( const VarSet &/*ns*/ ) { assert( 0 == 1 ); }
        virtual void makeCavity(const Var & /*n*/) { assert( 0 == 1 ); }
        bool checkProperties();
};


} // end of namespace dai


#endif
