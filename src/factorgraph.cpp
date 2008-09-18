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


#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <dai/factorgraph.h>
#include <dai/util.h>


namespace dai {


using namespace std;


FactorGraph::FactorGraph( const std::vector<Factor> &P ) : G(), _undoProbs(), _normtype(Prob::NORMPROB) {
    // add factors, obtain variables
    set<Var> _vars;
    factors.reserve( P.size() );
    size_t nrEdges = 0;
    for( vector<Factor>::const_iterator p2 = P.begin(); p2 != P.end(); p2++ ) {
        factors.push_back( *p2 );
        copy( p2->vars().begin(), p2->vars().end(), inserter( _vars, _vars.begin() ) );
        nrEdges += p2->vars().size();
    }

    // add _vars
    vars.reserve( _vars.size() );
    for( set<Var>::const_iterator p1 = _vars.begin(); p1 != _vars.end(); p1++ )
        vars.push_back( *p1 );
    
    // create graph structure
    createGraph( nrEdges );
}


/// Part of constructors (creates edges, neighbours and adjacency matrix)
void FactorGraph::createGraph( size_t nrEdges ) {
    // create a mapping for indices
    hash_map<size_t, size_t> hashmap;
    
    for( size_t i = 0; i < vars.size(); i++ )
        hashmap[var(i).label()] = i;
    
    // create edge list
    vector<Edge> edges;
    edges.reserve( nrEdges );
    for( size_t i2 = 0; i2 < nrFactors(); i2++ ) {
        const VarSet& ns = factor(i2).vars();
        for( VarSet::const_iterator q = ns.begin(); q != ns.end(); q++ )
            edges.push_back( Edge(hashmap[q->label()], i2) );
    }

    // create bipartite graph
    G.create( nrVars(), nrFactors(), edges.begin(), edges.end() );
}


ostream& operator << (ostream& os, const FactorGraph& fg) {
    os << fg.nrFactors() << endl;

    for( size_t I = 0; I < fg.nrFactors(); I++ ) {
        os << endl;
        os << fg.factor(I).vars().size() << endl;
        for( VarSet::const_iterator i = fg.factor(I).vars().begin(); i != fg.factor(I).vars().end(); i++ )
            os << i->label() << " ";
        os << endl;
        for( VarSet::const_iterator i = fg.factor(I).vars().begin(); i != fg.factor(I).vars().end(); i++ )
            os << i->states() << " ";
        os << endl;
        size_t nr_nonzeros = 0;
        for( size_t k = 0; k < fg.factor(I).states(); k++ )
            if( fg.factor(I)[k] != 0.0 )
                nr_nonzeros++;
        os << nr_nonzeros << endl;
        for( size_t k = 0; k < fg.factor(I).states(); k++ )
            if( fg.factor(I)[k] != 0.0 ) {
                char buf[20];
                sprintf(buf,"%18.14g", fg.factor(I)[k]);
                os << k << " " << buf << endl;
            }
    }

    return(os);
}


istream& operator >> (istream& is, FactorGraph& fg) {
    long verbose = 0;

    try {
        vector<Factor> factors;
        size_t nr_f;
        string line;
        
        while( (is.peek()) == '#' )
            getline(is,line);
        is >> nr_f;
        if( is.fail() )
            throw "ReadFromFile: unable to read number of Factors";
        if( verbose >= 2 )
            cout << "Reading " << nr_f << " factors..." << endl;

        getline (is,line);
        if( is.fail() )
            throw "ReadFromFile: empty line expected";

        for( size_t I = 0; I < nr_f; I++ ) {
            if( verbose >= 3 )
                cout << "Reading factor " << I << "..." << endl;
            size_t nr_members;
            while( (is.peek()) == '#' )
                getline(is,line);
            is >> nr_members;
            if( verbose >= 3 )
                cout << "  nr_members: " << nr_members << endl;

            vector<long> labels;
            for( size_t mi = 0; mi < nr_members; mi++ ) {
                long mi_label;
                while( (is.peek()) == '#' )
                    getline(is,line);
                is >> mi_label;
                labels.push_back(mi_label);
            }
            if( verbose >= 3 ) {
                cout << "  labels: ";
                copy (labels.begin(), labels.end(), ostream_iterator<int>(cout, " "));
                cout << endl;
            }

            vector<size_t> dims;
            for( size_t mi = 0; mi < nr_members; mi++ ) {
                size_t mi_dim;
                while( (is.peek()) == '#' )
                    getline(is,line);
                is >> mi_dim;
                dims.push_back(mi_dim);
            }
            if( verbose >= 3 ) {
                cout << "  dimensions: ";
                copy (dims.begin(), dims.end(), ostream_iterator<int>(cout, " "));
                cout << endl;
            }

            // add the Factor
            VarSet I_vars;
            for( size_t mi = 0; mi < nr_members; mi++ )
                I_vars |= Var(labels[mi], dims[mi]);
            factors.push_back(Factor(I_vars,0.0));
            
            // calculate permutation sigma (internally, members are sorted)
            vector<size_t> sigma(nr_members,0);
            VarSet::iterator j = I_vars.begin();
            for( size_t mi = 0; mi < nr_members; mi++,j++ ) {
                long search_for = j->label();
                vector<long>::iterator j_loc = find(labels.begin(),labels.end(),search_for);
                sigma[mi] = j_loc - labels.begin();
            }
            if( verbose >= 3 ) {
                cout << "  sigma: ";
                copy( sigma.begin(), sigma.end(), ostream_iterator<int>(cout," "));
                cout << endl;
            }

            // calculate multindices
            Permute permindex( dims, sigma );
            
            // read values
            size_t nr_nonzeros;
            while( (is.peek()) == '#' )
                getline(is,line);
            is >> nr_nonzeros;
            if( verbose >= 3 ) 
                cout << "  nonzeroes: " << nr_nonzeros << endl;
            for( size_t k = 0; k < nr_nonzeros; k++ ) {
                size_t li;
                double val;
                while( (is.peek()) == '#' )
                    getline(is,line);
                is >> li;
                while( (is.peek()) == '#' )
                    getline(is,line);
                is >> val;

                // store value, but permute indices first according
                // to internal representation
                factors.back()[permindex.convert_linear_index( li  )] = val;
            }
        }

        if( verbose >= 3 ) {
            cout << "factors:" << endl;
            copy(factors.begin(), factors.end(), ostream_iterator<Factor>(cout,"\n"));
        }

        fg = FactorGraph(factors);
    } catch (char *e) {
        cout << e << endl;
    }

    return is;
}


VarSet FactorGraph::delta( unsigned i ) const {
    // calculate Markov Blanket
    VarSet del;
    foreach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
        foreach( const Neighbor &j, nbF(I) ) // for all neighboring variables j of I
            if( j != i )
                del |= var(j);

    return del;
}


VarSet FactorGraph::Delta( unsigned i ) const {
    return( delta(i) | var(i) );
}


void FactorGraph::makeCavity( unsigned i ) {
    // fills all Factors that include var(i) with ones
    foreach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
        factor(I).fill( 1.0 );
}


bool FactorGraph::hasNegatives() const {
    bool result = false;
    for( size_t I = 0; I < nrFactors() && !result; I++ )
        if( factor(I).hasNegatives() )
            result = true;
    return result;
}
 

long FactorGraph::ReadFromFile(const char *filename) {
    ifstream infile;
    infile.open (filename);
    if (infile.is_open()) {
        infile >> *this;
        infile.close();
        return 0;
    } else {
        cout << "ERROR OPENING FILE" << endl;
        return 1;
    }
}


long FactorGraph::WriteToFile(const char *filename) const {
    ofstream outfile;
    outfile.open (filename);
    if (outfile.is_open()) {
        try {
            outfile << *this;
        } catch (char *e) {
            cout << e << endl;
            return 1;
        }
        outfile.close();
        return 0;
    } else {
        cout << "ERROR OPENING FILE" << endl;
        return 1;
    }
}


long FactorGraph::WriteToDotFile(const char *filename) const {
    ofstream outfile;
    outfile.open (filename);
    if (outfile.is_open()) {
        try {
            outfile << "graph G {" << endl;
            outfile << "graph[size=\"9,9\"];" << endl;
            outfile << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
            for( size_t i = 0; i < nrVars(); i++ )
                outfile << "\tx" << var(i).label() << ";" << endl;
            outfile << "node[shape=box,style=filled,color=lightgrey,width=0.3,height=0.3,fixedsize=true];" << endl;
            for( size_t I = 0; I < nrFactors(); I++ )
                outfile << "\tp" << I << ";" << endl;
            for( size_t i = 0; i < nrVars(); i++ )
                foreach( const Neighbor &I, nbV(i) )  // for all neighboring factors I of i
                    outfile << "\tx" << var(i).label() << " -- p" << I << ";" << endl;
            outfile << "}" << endl;
        } catch (char *e) {
            cout << e << endl;
            return 1;
        }
        outfile.close();
        return 0;
    } else {
        cout << "ERROR OPENING FILE" << endl;
        return 1;
    }
}


bool hasShortLoops( const vector<Factor> &P ) {
    bool found = false;
    vector<Factor>::const_iterator I, J;
    for( I = P.begin(); I != P.end(); I++ ) {
        J = I;
        J++;
        for( ; J != P.end(); J++ )
            if( (I->vars() & J->vars()).size() >= 2 ) {
                found = true;
                break;
            }
        if( found )
            break;
    }
    return found;
}


void RemoveShortLoops(vector<Factor> &P) {
    bool found = true;
    while( found ) {
        found = false;
        vector<Factor>::iterator I, J;
        for( I = P.begin(); I != P.end(); I++ ) {
            J = I;
            J++;
            for( ; J != P.end(); J++ )
                if( (I->vars() & J->vars()).size() >= 2 ) {
                    found = true;
                    break;
                }
            if( found )
                break;
        }
        if( found ) {
            cout << "Merging factors " << I->vars() << " and " << J->vars() << endl;
            *I *= *J;
            P.erase(J);
        }
    }
}


vector<VarSet> FactorGraph::Cliques() const {
    vector<VarSet> result;
    
    for( size_t I = 0; I < nrFactors(); I++ ) {
        bool maximal = true;
        for( size_t J = 0; (J < nrFactors()) && maximal; J++ )
            if( (factor(J).vars() >> factor(I).vars()) && !(factor(J).vars() == factor(I).vars()) )
                maximal = false;
        
        if( maximal )
            result.push_back( factor(I).vars() );
    }

    return result;
}


void FactorGraph::clamp( const Var & n, size_t i ) {
    assert( i <= n.states() );

    // Multiply each factor that contains the variable with a delta function

    Factor delta_n_i(n,0.0);
    delta_n_i[i] = 1.0;

    // For all factors that contain n
    for( size_t I = 0; I < nrFactors(); I++ ) 
        if( factor(I).vars().contains( n ) )
            // Multiply it with a delta function
            factor(I) *= delta_n_i;

    return;
}


void FactorGraph::saveProb( size_t I ) {
    map<size_t,Prob>::iterator it = _undoProbs.find( I );
    if( it != _undoProbs.end() )
        cout << "FactorGraph::saveProb:  WARNING: _undoProbs[I] already defined!" << endl;
    _undoProbs[I] = factor(I).p();
}


void FactorGraph::undoProb( size_t I ) {
    map<size_t,Prob>::iterator it = _undoProbs.find( I );
    if( it != _undoProbs.end() ) {
        factor(I).p() = (*it).second;
        _undoProbs.erase(it);
    }
}


void FactorGraph::saveProbs( const VarSet &ns ) {
    if( !_undoProbs.empty() )
        cout << "FactorGraph::saveProbs:  WARNING: _undoProbs not empy!" << endl;
    for( size_t I = 0; I < nrFactors(); I++ )
        if( factor(I).vars().intersects( ns ) )
            _undoProbs[I] = factor(I).p();
}


void FactorGraph::undoProbs( const VarSet &ns ) {
    for( map<size_t,Prob>::iterator uI = _undoProbs.begin(); uI != _undoProbs.end(); ) {
        if( factor((*uI).first).vars().intersects( ns ) ) {
//          cout << "undoing " << factor((*uI).first).vars() << endl;
//          cout << "from " << factor((*uI).first).p() << " to " << (*uI).second << endl;
            factor((*uI).first).p() = (*uI).second;
            _undoProbs.erase(uI++);
        } else
            uI++;
    }
}


} // end of namespace dai
