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


#ifndef __defined_libdai_bp_h
#define __defined_libdai_bp_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>


namespace dai {


class BP : public DAIAlgFG {
    protected:
        typedef std::vector<size_t> ind_t;
        struct EdgeProp {
            ind_t  index;
            Prob   message;
            Prob   newMessage;
            double residual;
        };
        std::vector<std::vector<EdgeProp> > edges;
    
    public:
        struct Properties {
            size_t verbose;
            size_t maxiter;
            double tol;
            bool logdomain;
            ENUM4(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL)
            UpdateType updates;
        } props;
        double maxdiff;

    public:
        /// Default constructor
        BP() : DAIAlgFG(), edges(), props(), maxdiff(0.0) {};
        /// Copy constructor
        BP( const BP & x ) : DAIAlgFG(x), edges(x.edges), props(x.props), maxdiff(x.maxdiff) {};
        /// Clone *this
        BP* clone() const { return new BP(*this); }
        /// Construct from FactorGraph fg and PropertySet opts
        BP( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg), edges(), props(), maxdiff(0.0) {
            setProperties( opts );
            create();
        }
        /// Assignment operator
        BP& operator=( const BP & x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                edges = x.edges;
                props = x.props;
                maxdiff = x.maxdiff;
            }
            return *this;
        }

        static const char *Name;

        Prob & message(size_t i, size_t _I) { return edges[i][_I].message; }
        const Prob & message(size_t i, size_t _I) const { return edges[i][_I].message; }
        Prob & newMessage(size_t i, size_t _I) { return edges[i][_I].newMessage; }
        const Prob & newMessage(size_t i, size_t _I) const { return edges[i][_I].newMessage; }
        ind_t & index(size_t i, size_t _I) { return edges[i][_I].index; }
        const ind_t & index(size_t i, size_t _I) const { return edges[i][_I].index; }
        double & residual(size_t i, size_t _I) { return edges[i][_I].residual; }
        const double & residual(size_t i, size_t _I) const { return edges[i][_I].residual; }

        std::string identify() const;
        void create();
        void init();
        double run();

        void findMaxResidual( size_t &i, size_t &_I );
        void calcNewMessage( size_t i, size_t _I );
        Factor beliefV (size_t i) const;
        Factor beliefF (size_t I) const;
        Factor belief (const Var &n) const;
        Factor belief (const VarSet &n) const;
        std::vector<Factor> beliefs() const;
        Real logZ() const;

        void init( const VarSet &ns );
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs(ns); init(ns); }

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        double maxDiff() const { return maxdiff; }
};


} // end of namespace dai


#endif
