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
/// \brief Defines the FactorGraph class


#ifndef __defined_libdai_factorgraph_h
#define __defined_libdai_factorgraph_h


#include <iostream>
#include <map>
#include <dai/bipgraph.h>
#include <dai/factor.h>


namespace dai {


/// Represents a factor graph.
/** Both Bayesian Networks and Markov random fields can be represented in a 
 *  unifying representation, called <em>factor graph</em> [\ref KFL01], 
 *  implemented in libDAI by the FactorGraph class.
 *  
 *  Consider a probability distribution over \f$N\f$ discrete random variables 
 *  \f$x_0,x_1,\dots,x_N\f$ that factorizes as a product of factors, each of
 *  which depends on some subset of the variables:
 *  \f[
 *    P(x_0,x_1,\dots,x_N) = \frac{1}{Z} \prod_{I=0}^M f_I(x_I), \qquad
 *    Z = \sum_{x_0}\dots\sum_{x_N} \prod_{I=0}^M f_I(X_I).
 *  \f]
 *  Each factor \f$f_I\f$ is a function from an associated subset
 *  of variables \f$X_I \subset \{x_0,x_1,\dots,x_N\}\f$ to the nonnegative
 *  real numbers.
 * 
 *  For a Bayesian network, each factor corresponds to a (conditional) 
 *  probability table, whereas for a Markov random field, each factor 
 *  corresponds to a maximal clique of the undirected graph.
 *
 *  Factor graphs explicitly express the factorization structure of the
 *  corresponding probability distribution.
 */ 
class FactorGraph {
    public:
        /// Stores the neighborhood structure
        BipartiteGraph                    G;

        /// Shorthand for BipartiteGraph::Neighbor
        typedef BipartiteGraph::Neighbor  Neighbor;

        /// Shorthand for BipartiteGraph::Neighbors
        typedef BipartiteGraph::Neighbors Neighbors;

        /// Shorthand for BipartiteGraph::Edge
        typedef BipartiteGraph::Edge      Edge;

    private:
        std::vector<Var>         _vars;
        std::vector<Factor>      _factors;
        std::map<size_t,Factor>  _backup;

    public:
        /// Default constructor
        FactorGraph() : G(), _vars(), _factors(), _backup() {}

        /// Copy constructor
        FactorGraph(const FactorGraph & x) : G(x.G), _vars(x._vars), _factors(x._factors), _backup(x._backup) {}

        /// Assignment operator
        FactorGraph & operator=(const FactorGraph & x) {
            if( this != &x ) {
                G          = x.G;
                _vars      = x._vars;
                _factors   = x._factors;
                _backup    = x._backup;
            }
            return *this;
        }

        /// Constructs a FactorGraph from a vector of factors
        FactorGraph(const std::vector<Factor> &P);

        /// Constructs a FactorGraph from given factor and variable iterators
        /** \tparam FactorInputIterator Iterator with value_type Factor
         *  \tparam VarInputIterator Iterator with value_type Var
         *  \pre Assumes that the set of variables in [var_begin,var_end) is the union of the variables in the factors in [fact_begin, fact_end)
         */
        template<typename FactorInputIterator, typename VarInputIterator>
        FactorGraph(FactorInputIterator fact_begin, FactorInputIterator fact_end, VarInputIterator var_begin, VarInputIterator var_end, size_t nr_fact_hint = 0, size_t nr_var_hint = 0 );

        /// Destructor        
        virtual ~FactorGraph() {}

        /// Clone *this (virtual copy constructor)
        virtual FactorGraph* clone() const { return new FactorGraph(); }

        /// Create (virtual default constructor)
        virtual FactorGraph* create() const { return new FactorGraph(*this); }

        /// Returns const reference to i'th variable
        const Var & var(size_t i) const { return _vars[i]; }
        /// Returns const reference to all factors
        const std::vector<Var> & vars() const { return _vars; }
        /// Returns reference to I'th factor
        Factor & factor(size_t I) { return _factors[I]; }
        /// Returns const reference to I'th factor
        const Factor & factor(size_t I) const { return _factors[I]; }
        /// Returns const reference to all factors
        const std::vector<Factor> & factors() const { return _factors; }

        /// Returns number of variables
        size_t nrVars() const { return vars().size(); }
        /// Returns number of factors
        size_t nrFactors() const { return factors().size(); }
        /// Calculates number of edges
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

        /// Returns the index of a particular variable
        size_t findVar( const Var & n ) const {
            size_t i = find( vars().begin(), vars().end(), n ) - vars().begin();
            assert( i != nrVars() );
            return i;
        }

        /// Returns a set of indexes corresponding to a set of variables
        std::set<size_t> findVars( VarSet &ns ) const {
            std::set<size_t> indexes;
            for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
                indexes.insert( findVar( *n ) );
            return indexes;
        }

        /// Returns index of the first factor that depends on the variables
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

        /// Returns true if the FactorGraph is connected
        bool isConnected() const { return G.isConnected(); }

        /// Returns true if the FactorGraph is a tree
        bool isTree() const { return G.isTree(); }

        /// Returns true if each factor depends on at most two variables
        bool isPairwise() const;

        /// Returns true if each variable has only two possible values
        bool isBinary() const;

        /// Reads a FactorGraph from a file
        void ReadFromFile(const char *filename);

        /// Writes a FactorGraph to a file
        void WriteToFile(const char *filename) const;

        /// Writes a FactorGraph to a GraphViz .dot file
        void printDot( std::ostream& os ) const;
        
        /// Returns the cliques in this FactorGraph
        std::vector<VarSet> Cliques() const;

        /// Clamp variable v_i to value state (i.e. multiply with a Kronecker delta \f$\delta_{x_{v_i},x}\f$);
        /** This version changes the factor graph structure and thus returns a newly constructed FactorGraph
         *  and keeps the current one constant, contrary to clamp()
         */
        FactorGraph clamped( const Var & v_i, size_t x ) const;

        /// Returns a copy of *this, where all factors that are subsumed by some larger factor are merged with the larger factors.
        FactorGraph maximalFactors() const;

        /// Makes a backup of the I'th Factor
        void restoreFactor( size_t I );

        /// Restores the I'th Factor from the backup (it should be backed up first)
        void backupFactor( size_t I );

        /// Makes a backup of all factors connected to a set of variables
        void backupFactors( const VarSet &ns );
        /// Restores all factors connected to a set of variables from their backups
        void restoreFactors( const VarSet &ns );

        // Friends
        friend std::ostream& operator << (std::ostream& os, const FactorGraph& fg);
        friend std::istream& operator >> (std::istream& is, FactorGraph& fg);

    private:
        /// Part of constructors (creates edges, neighbors and adjacency matrix)
        void constructGraph( size_t nrEdges );
};


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
    _vars.reserve( nr_var_hint );
    for( VarInputIterator p1 = var_begin; p1 != var_end; ++p1 )
        _vars.push_back( *p1 );

    // create graph structure
    constructGraph( nrEdges );
}


} // end of namespace dai


#endif
