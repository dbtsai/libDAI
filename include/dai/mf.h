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


#ifndef __defined_libdai_mf_h
#define __defined_libdai_mf_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


class MF : public DAIAlgFG {
    private:
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

        /// Construct from FactorGraph fg and PropertySet opts
        MF( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), _beliefs(), _maxdiff(0.0), _iters(0U), props() {
            setProperties( opts );
            construct();
        }

        /// Copy constructor
        MF( const MF &x ) : DAIAlgFG(x), _beliefs(x._beliefs), _maxdiff(x._maxdiff), _iters(x._iters), props(x.props) {}

        /// Clone *this (virtual copy constructor)
        virtual MF* clone() const { return new MF(*this); }

        /// Create (virtual default constructor)
        virtual MF* create() const { return new MF(); }

        /// Assignment operator
        MF& operator=( const MF &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _beliefs = x._beliefs;
                _maxdiff = x._maxdiff;
                _iters   = x._iters;
                props    = x.props;
            }
            return *this;
        }

        /// Identifies itself for logging purposes
        virtual std::string identify() const;

        /// Get single node belief
        virtual Factor belief( const Var &n ) const;

        /// Get general belief
        virtual Factor belief( const VarSet &ns ) const;

        /// Get all beliefs
        virtual std::vector<Factor> beliefs() const;

        /// Get log partition sum
        virtual Real logZ() const;

        /// Clear messages and beliefs
        virtual void init();

        /// Clear messages and beliefs corresponding to the nodes in ns
        virtual void init( const VarSet &ns );

        /// The actual approximate inference algorithm
        virtual double run();

        /// Return maximum difference between single node beliefs in the last pass
        virtual double maxDiff() const { return _maxdiff; }

        /// Return number of passes over the factorgraph
        virtual size_t Iterations() const { return _iters; }


        void construct();

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;

        Factor beliefV( size_t i ) const;
};


} // end of namespace dai


#endif
