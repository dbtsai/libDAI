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


#ifndef __defined_libdai_exactinf_h
#define __defined_libdai_exactinf_h


#include <dai/daialg.h>
#include <dai/properties.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>


namespace dai {


class ExactInf : public DAIAlgFG {
    public:
        struct Properties {
            size_t verbose;
        } props;

    private:
        std::vector<Factor> _beliefsV;
        std::vector<Factor> _beliefsF;
        Real                _logZ;

    public:
        /// Default constructor
        ExactInf() : DAIAlgFG(), props(), _beliefsV(), _beliefsF(), _logZ(0) {}
        
        /// Copy constructor
        ExactInf( const ExactInf &x ) : DAIAlgFG(x), props(x.props), _beliefsV(x._beliefsV), _beliefsF(x._beliefsF), _logZ(x._logZ) {}

        /// Clone (virtual copy constructor)
        virtual ExactInf* clone() const {
            return new ExactInf(*this);
        }

        /// Construct from FactorGraph fg and PropertySet opts
        ExactInf( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), props(), _beliefsV(), _beliefsF(), _logZ() {
            setProperties( opts );
            create();
        }
        
        /// Assignment operator
        ExactInf& operator=( const ExactInf &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                props     = x.props;
                _beliefsV = x._beliefsV;
                _beliefsF = x._beliefsF;
                _logZ     = x._logZ;
            }
            return *this;
        }

/*        /// Create (virtual constructor)
        virtual ExactInf* create() const {
            return new ExactInf();
        }
*/
        /// Return maximum difference between single node 
        /// beliefs for two consecutive iterations
        virtual double maxDiff() const {
            DAI_THROW(NOT_IMPLEMENTED);
            return 0.0;
        }
        
        /// Identifies itself for logging purposes
        virtual std::string identify() const;

        /// Get single node belief
        virtual Factor belief( const Var &n ) const {
            return beliefV( findVar( n ) ); 
        }

        /// Get general belief
        virtual Factor belief( const VarSet &n ) const;

        /// Get all beliefs
        virtual std::vector<Factor> beliefs() const;

        /// Get log partition sum
        virtual Real logZ() const {
            return _logZ;
        }

        /// Clear messages and beliefs
        virtual void init();

        /// Clear messages and beliefs corresponding to the nodes in ns
        virtual void init( const VarSet &/*ns*/ ) {}

        /// The actual approximate inference algorithm
        virtual double run();

        /// Name of this inference method
        static const char *Name;

        void create();
        void restoreFactors( const VarSet &ns ) { FactorGraph::restoreFactors(ns); init(ns); }
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;

        Factor beliefV( size_t i ) const { 
            return _beliefsV[i]; 
        }

        Factor beliefF( size_t I ) const { 
            return _beliefsF[I]; 
        }
};


} // end of namespace dai


#endif
