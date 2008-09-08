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
        ENUM4(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL)
        UpdateType Updates() const { return GetPropertyAs<UpdateType>("updates"); }

        // default constructor
        BP() : DAIAlgFG() {};
        // copy constructor
        BP(const BP & x) : DAIAlgFG(x), edges(x.edges) {};
        BP* clone() const { return new BP(*this); }
        // construct BP object from FactorGraph
        BP(const FactorGraph & fg, const Properties &opts) : DAIAlgFG(fg, opts) {
            assert( checkProperties() );
            create();
        }
        // assignment operator
        BP & operator=(const BP & x) {
            if(this!=&x) {
                DAIAlgFG::operator=(x);
                edges = x.edges;
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
        void findMaxResidual( size_t &i, size_t &_I );

        std::string identify() const;
        void create();
        void init();
        void calcNewMessage( size_t i, size_t _I );
        double run();
        Factor beliefV (size_t i) const;
        Factor beliefF (size_t I) const;
        Factor belief (const Var &n) const;
        Factor belief (const VarSet &n) const;
        std::vector<Factor> beliefs() const;
        Complex logZ() const;

        void init( const VarSet &ns );
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs(ns); init(ns); }
        bool checkProperties();
};


} // end of namespace dai


#endif
