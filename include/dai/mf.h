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
        /// Maximum difference encountered so far
        double _maxdiff;
        /// Number of iterations needed
        size_t _iters;

    public:
        struct Properties {
            size_t verbose;
            size_t maxiter;
            double tol;
            double damping;
        } props;
        static const char *Name;

    public:
        /// Default constructor
        MF() : DAIAlgFG(), _beliefs(), _maxdiff(0.0), _iters(0U), props() {}

        // construct MF object from FactorGraph
        MF( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), _beliefs(), _maxdiff(0.0), _iters(0U), props() {
            setProperties( opts );
            construct();
        }

        /// Copy constructor
        MF( const MF &x ) : DAIAlgFG(x), _beliefs(x._beliefs), _maxdiff(x._maxdiff), _iters(x._iters), props(x.props) {}

        /// Assignment operator
        MF & operator=( const MF &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _beliefs = x._beliefs;
                _maxdiff = x._maxdiff;
                _iters   = x._iters;
                props    = x.props;
            }
            return *this;
        }

        /// Clone *this (virtual copy constructor)
        virtual MF* clone() const { return new MF(*this); }

        /// Create (virtual constructor)
        virtual MF* create() const { return new MF(); }

        /// Return number of passes over the factorgraph needed
        virtual size_t Iterations() const { return _iters; }

        /// Return maximum difference between single node beliefs for two consecutive iterations
        double maxDiff() const { return _maxdiff; }

        /// Identify *this for logging purposes
        std::string identify() const;

        /// Get single node belief
        Factor belief( const Var &n ) const;

        /// Get general belief
        Factor belief( const VarSet &ns ) const;

        /// Get all beliefs
        std::vector<Factor> beliefs() const;

        /// Get log partition sum
        Real logZ() const;

        void construct();

        void init();

        /// Clear messages and beliefs corresponding to the nodes in ns
        virtual void init( const VarSet &ns );

        /// The actual approximate inference algorithm
        double run();

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;

        Factor beliefV( size_t i ) const;
};


} // end of namespace dai


#endif
