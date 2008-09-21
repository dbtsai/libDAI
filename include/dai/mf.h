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


#ifndef __defined_libdai_mf_h
#define __defined_libdai_mf_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


class MF : public DAIAlgFG {
    protected:
        std::vector<Factor>  _beliefs;

    public:
        struct Properties {
            size_t verbose;
            size_t maxiter;
            double tol;
        } props;
        double maxdiff;
        
    public:
        // default constructor
        MF() : DAIAlgFG(), _beliefs(), props(), maxdiff(0.0) {}
        // copy constructor
        MF( const MF& x ) : DAIAlgFG(x), _beliefs(x._beliefs), props(x.props), maxdiff(x.maxdiff) {}
        MF* clone() const { return new MF(*this); }
        // construct MF object from FactorGraph
        MF( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg), _beliefs(), props(), maxdiff(0.0) {
            setProperties( opts );
            create();
        }
        // assignment operator
        MF& operator=( const MF &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _beliefs = x._beliefs;
                props = x.props;
                maxdiff = x.maxdiff;
            }
            return *this;
        }

        static const char *Name;
        std::string identify() const;
        void create();
        void init();
        double run();
        Factor beliefV (size_t i) const;
        Factor belief (const Var &n) const;
        Factor belief (const VarSet &ns) const;
        std::vector<Factor> beliefs() const;
        Complex logZ() const;

        void init( const VarSet &ns );
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs(ns); init(ns); }
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        double maxDiff() const { return maxdiff; }
};


} // end of namespace dai


#endif
