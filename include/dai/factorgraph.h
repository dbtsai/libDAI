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


#ifndef __defined_libdai_factorgraph_h
#define __defined_libdai_factorgraph_h


#include <iostream>
#include <map>
#include <dai/bipgraph.h>
#include <dai/factor.h>


namespace dai {


class FactorGraph {
    public:
        BipartiteGraph         G;
        std::vector<Var>       vars;
        typedef BipartiteGraph::Neighbor  Neighbor;
        typedef BipartiteGraph::Neighbors Neighbors;
        typedef BipartiteGraph::Edge      Edge;

    private:
        std::vector<Factor>      _factors;
        std::map<size_t,Factor>  _backup;

    public:
        /// Default constructor
        FactorGraph() : G(), vars(), _factors(), _backup() {}
        /// Copy constructor
        FactorGraph(const FactorGraph & x) : G(x.G), vars(x.vars), _factors(x._factors), _backup(x._backup) {}
        /// Construct FactorGraph from vector of Factors
        FactorGraph(const std::vector<Factor> &P);
        // Construct a FactorGraph from given factor and variable iterators
        template<typename FactorInputIterator, typename VarInputIterator>
        FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint = 0, size_t nr_var_hint = 0 );
        
        /// Assignment operator
        FactorGraph & operator=(const FactorGraph & x) {
            if( this != &x ) {
                G          = x.G;
                vars       = x.vars;
                _factors   = x._factors;
                _backup    = x._backup;
            }
            return *this;
        }
        virtual ~FactorGraph() {}

        /// Create (virtual default constructor)
        virtual FactorGraph* create() const { return new FactorGraph(*this); }

        /// Clone (virtual copy constructor)
        virtual FactorGraph* clone() const { return new FactorGraph(); }

        // aliases
        Var & var(size_t i) { return vars[i]; }
        /// Get const reference to i'th variable
        const Var & var(size_t i) const { return vars[i]; }
        /// Get const reference to I'th factor
        Factor & factor(size_t I) { return _factors[I]; }
        /// Get const reference to I'th factor
        const Factor & factor(size_t I) const { return _factors[I]; }
        /// Get const reference to all factors
        const std::vector<Factor> & factors() const { return _factors; }

        /// Get number of variables
        size_t nrVars() const { return vars.size(); }
        /// Get number of factors
        size_t nrFactors() const { return _factors.size(); }
        size_t nrEdges() const { return G.nrEdges(); }

        /// Provides read access to neighbors of variable
        const Neighbors & nbV( size_t i ) const { return G.nb1(i); }
        /// Provides full access to neighbors of variable
        Neighbors & nbV( size_t i ) { return G.nb1(i); }
        /// Provides read access to neighbors of factor
        const Neighbors & nbF( size_t I ) const { return G.nb2(I); }
        /// Provides full access to neighbors of factor
        Neighbors & nbF( size_t I ) { return G.nb2(I); }
        /// Provides read access to neighbor of variable
        const Neighbor & nbV( size_t i, size_t _I ) const { return G.nb1(i)[_I]; }
        /// Provides full access to neighbor of variable
        Neighbor & nbV( size_t i, size_t _I ) { return G.nb1(i)[_I]; }
        /// Provides read access to neighbor of factor
        const Neighbor & nbF( size_t I, size_t _i ) const { return G.nb2(I)[_i]; }
        /// Provides full access to neighbor of factor
        Neighbor & nbF( size_t I, size_t _i ) { return G.nb2(I)[_i]; }

        /// Get index of variable n
        size_t findVar( const Var & n ) const {
            size_t i = find( vars.begin(), vars.end(), n ) - vars.begin();
            assert( i != nrVars() );
            return i;
        }

        /// Get set of indexes for set of variables
        std::set<size_t> findVars( VarSet &ns ) const {
            std::set<size_t> indexes;
            for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
                indexes.insert( findVar( *n ) );
            return indexes;
        }

        /// Get index of first factor involving ns
        size_t findFactor(const VarSet &ns) const {
            size_t I;
            for( I = 0; I < nrFactors(); I++ )
                if( factor(I).vars() == ns )
                    break;
            assert( I != nrFactors() );
            return I;
        }

        /// Return all variables that occur in a factor involving the i'th variable, itself included
        VarSet Delta( unsigned i ) const;

        /// Return all variables that occur in a factor involving some variable in ns, ns itself included
        VarSet Delta( const VarSet &ns ) const;

        /// Return all variables that occur in a factor involving the i'th variable, n itself excluded
        VarSet delta( unsigned i ) const;

        /// Return all variables that occur in a factor involving some variable in ns, ns itself excluded
        VarSet delta( const VarSet & ns ) const {
            return Delta( ns ) / ns;
        }

        /// Set the content of the I'th factor and make a backup of its old content if backup == true
        virtual void setFactor( size_t I, const Factor &newFactor, bool backup = false ) {
            assert( newFactor.vars() == factor(I).vars() ); 
            if( backup )
                backupFactor( I );
            _factors[I] = newFactor; 
        }

        /// Set the contents of all factors as specified by facs and make a backup of the old contents if backup == true
        virtual void setFactors( const std::map<size_t, Factor> & facs, bool backup = false ) {
            for( std::map<size_t, Factor>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ ) {
                if( backup )
                    backupFactor( fac->first );
                setFactor( fac->first, fac->second );
            }
        }

        /// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$);
        /// If backup == true, make a backup of all factors that are changed
        virtual void clamp( const Var & n, size_t i, bool backup = false );

        /// Set all factors interacting with the i'th variable 1
        virtual void makeCavity( unsigned i, bool backup = false );

        /// Backup the factors specified by indices in facs
        virtual void backupFactors( const std::set<size_t> & facs );

        /// Restore all factors to the backup copies
        virtual void restoreFactors();

        bool isConnected() const { return G.isConnected(); }
        bool isTree() const { return G.isTree(); }

        friend std::ostream& operator << (std::ostream& os, const FactorGraph& fg);
        friend std::istream& operator >> (std::istream& is, FactorGraph& fg);

        void ReadFromFile(const char *filename);
        void WriteToFile(const char *filename) const;
        void printDot( std::ostream& os ) const;
        
        std::vector<VarSet> Cliques() const;

        // Clamp variable v_i to value state (i.e. multiply with a Kronecker delta \f$\delta_{x_{v_i},x}\f$);
        // This version changes the factor graph structure and thus returns a newly constructed FactorGraph
        // and keeps the current one constant, contrary to clamp()
        FactorGraph clamped( const Var & v_i, size_t x ) const;

        FactorGraph maximalFactors() const;

        bool isPairwise() const;
        bool isBinary() const;

        void restoreFactor( size_t I );
        void backupFactor( size_t I );
        void restoreFactors( const VarSet &ns );
        void backupFactors( const VarSet &ns );
        /// Part of constructors (creates edges, neighbors and adjacency matrix)
        void constructGraph( size_t nrEdges );
};


// assumes that the set of variables in [var_begin,var_end) is the union of the variables in the factors in [fact_begin, fact_end)
template<typename FactorInputIterator, typename VarInputIterator>
FactorGraph::FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint, size_t nr_var_hint ) : G(), _backup() {
    // add factors
    size_t nrEdges = 0;
    _factors.reserve( nr_fact_hint );
    for( FactorInputIterator p2 = fact_begin; p2 != fact_end; ++p2 ) {
        _factors.push_back( *p2 );
        nrEdges += p2->vars().size();
    }

    // add variables
    vars.reserve( nr_var_hint );
    for( VarInputIterator p1 = var_begin; p1 != var_end; ++p1 )
        vars.push_back( *p1 );

    // create graph structure
    constructGraph( nrEdges );
}


} // end of namespace dai


#endif
