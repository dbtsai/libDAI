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
#include "factorgraph.h"


using namespace std;


FactorGraph::FactorGraph( const vector<Factor> &P ) : BipartiteGraph<Var,Factor>(), _undoProbs(), _hasNegatives(false), _normtype(Prob::NORMPROB) {
    // add Factors
    set<Var> _vars;
    for(vector<Factor>::const_iterator p2 = P.begin(); p2 != P.end(); p2++ ) {
        V2s().push_back(*p2);
        if( !_hasNegatives && p2->hasNegatives() )
            _hasNegatives = true;
        for( set<Var>::const_iterator i = p2->vars().begin(); i != p2->vars().end(); i++ )
            _vars.insert(*i);
    }

    // if negative factors are present, use LINF norm
    if( _hasNegatives )
        _normtype = Prob::NORMLINF;
    
    // add _vars
    for(VarSet::const_iterator p1 = _vars.begin(); p1 != _vars.end(); p1++ )
        vars().push_back(*p1);

    // create edges
    for(size_t i2 = 0; i2 < nrFactors(); i2++ ) {
        VarSet ns = factor(i2).vars();
        for(VarSet::const_iterator q = ns.begin(); q != ns.end(); q++ ) {
            for(size_t i1=0; i1 < nrVars(); i1++ ) {
                if (var(i1) == *q) {
                    edges().push_back(_edge_t(i1,i2));
                    break;
                }
            }
        }
    }

    // calc neighbours and adjacency matrix
    Regenerate();

    // Check for short loops
    if( hasShortLoops(P) )
        cerr << "FactorGraph::FactorGraph():  WARNING: short loops are present" << endl;
}


/*FactorGraph& FactorGraph::addFactor( const Factor &I ) {
    // add Factor
    _V2.push_back( I );
    if( !_hasNegatives && I.hasNegatives() )
        _hasNegatives = true;

    // if negative factors are present, use LINF norm
    if( _hasNegatives )
        _normtype = Prob::NORMLINF;

    // add new vars in Factor
    for( VarSet::const_iterator i = I.vars().begin(); i != I.vars().end(); i++ ) {
        size_t i_ind = find(vars().begin(), vars().end(), *i) - vars().begin();
        if( i_ind == vars().size() )
            _V1.push_back( *i );
        _E12.push_back( _edge_t( i_ind, nrFactors() - 1 ) );
    }

    Regenerate();
    return(*this);
}*/


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
        for( size_t k = 0; k < fg.factor(I).stateSpace(); k++ )
            if( fg.factor(I)[k] != 0.0 )
                nr_nonzeros++;
        os << nr_nonzeros << endl;
        for( size_t k = 0; k < fg.factor(I).stateSpace(); k++ )
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
                I_vars.insert( Var(labels[mi], dims[mi]) );
            factors.push_back(Factor(I_vars,0.0));
            
            // calculate permutation sigma (internally, members are sorted)
            vector<long> sigma(nr_members,0);
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
            vector<size_t> sdims(nr_members,0);
            for( size_t k = 0; k < nr_members; k++ ) {
                sdims[k] = dims[sigma[k]];
            }
            multind mi(dims);
            multind smi(sdims);
            if( verbose >= 3 ) {
                cout << "  mi.max(): " << mi.max() << endl;
                cout << "       ";
                for( size_t k=0; k < nr_members; k++ ) 
                    cout << labels[k] << " ";
                cout << "   ";
                for( size_t k=0; k < nr_members; k++ ) 
                    cout << labels[sigma[k]] << " ";
                cout << endl;
            }
            
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

                vector<size_t> vi = mi.vi(li);
                vector<size_t> svi(vi.size(),0);
                for( size_t k = 0; k < vi.size(); k++ )
                    svi[k] = vi[sigma[k]];
                size_t sli = smi.li(svi);
                if( verbose >= 3 ) {
                    cout << "    " << li << ": ";
                    copy( vi.begin(), vi.end(), ostream_iterator<size_t>(cout," "));
                    cout << "-> ";
                    copy( svi.begin(), svi.end(), ostream_iterator<size_t>(cout," "));
                    cout << ": " << sli << endl;
                }
                factors.back()[sli] = val;
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


VarSet FactorGraph::delta(const Var & n) const {
    // calculate Markov Blanket
    size_t i = findVar( n );

    VarSet del;
    for( _nb_cit I = nb1(i).begin(); I != nb1(i).end(); I++ )
        for( _nb_cit j = nb2(*I).begin(); j != nb2(*I).end(); j++ )
            if( *j != i )
                del |= var(*j);

    return del;
}


VarSet FactorGraph::Delta(const Var & n) const {
    return( delta(n) | n );
}


void FactorGraph::makeFactorCavity(size_t I) {
    // fill Factor I with ones
    factor(I).fill(1.0);
}


void FactorGraph::makeCavity(const Var & n) {
    // fills all Factors that include Var n with ones
    size_t i = findVar( n );

    for( _nb_cit I = nb1(i).begin(); I != nb1(i).end(); I++ )
        factor(*I).fill(1.0);
}


/*FactorGraph & FactorGraph::DeleteFactor(size_t I) {
    // Go through all edges
    for( vector<_edge_t>::iterator edge = _E12.begin(); edge != _E12.end(); edge++ )
        if( edge->second >= I ) {
            if( edge->second == I )
                edge->second = -1UL;
            else 
                (edge->second)--;
        }
    // Remove all edges containing I
    for( vector<_edge_t>::iterator edge = _E12.begin(); edge != _E12.end(); edge++ )
        if( edge->second == -1UL )
            edge = _E12.erase( edge );
//  vector<_edge_t>::iterator new_end = _E12.remove_if( _E12.begin(), _E12.end(), compose1( bind2nd(equal_to<size_t>(), -1), select2nd<_edge_t>() ) );
//  _E12.erase( new_end, _E12.end() );

    // Erase the factor
    _V2.erase( _V2.begin() + I );
    
    Regenerate();

    return *this;
}


FactorGraph & FactorGraph::DeleteVar(size_t i) {
    // Go through all edges
    for( vector<_edge_t>::iterator edge = _E12.begin(); edge != _E12.end(); edge++ )
        if( edge->first >= i ) {
            if( edge->first == i )
                edge->first = -1UL;
            else 
                (edge->first)--;
        }
    // Remove all edges containing i
    for( vector<_edge_t>::iterator edge = _E12.begin(); edge != _E12.end(); edge++ )
        if( edge->first == -1UL )
            edge = _E12.erase( edge );
                
//  vector<_edge_t>::iterator new_end = _E12.remove_if( _E12.begin(), _E12.end(), compose1( bind2nd(equal_to<size_t>(), -1), select1st<_edge_t>() ) );
//  _E12.erase( new_end, _E12.end() );

    // Erase the variable
    _V1.erase( _V1.begin() + i );
    
    Regenerate();

    return *this;
}*/


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
            for( size_t iI = 0; iI < nr_edges(); iI++ )
                outfile << "\tx" << var(edge(iI).first).label() << " -- p" << edge(iI).second << ";" << endl;
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


Factor FactorGraph::ExactMarginal(const VarSet & x) const {
    Factor P;
    for( size_t I = 0; I < nrFactors(); I++ )
        P *= factor(I);
    return P.marginal(x);
}


Real FactorGraph::ExactlogZ() const {
    Factor P;
    for( size_t I = 0; I < nrFactors(); I++ )
        P *= factor(I);
    return std::log(P.totalSum());
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

/*  if( do_surgery ) {
        size_t ni = find( vars().begin(), vars().end(), n) - vars().begin();

        if( ni != nrVars() ) {
            for( _nb_cit I = nb1(ni).begin(); I != nb1(ni).end(); I++ ) {
                if( factor(*I).size() == 1 )
                    // Remove this single-variable factor
    //              I = (_V2.erase(I))--;
                    _E12.erase( _E12.begin() + VV2E(ni, *I) );
                else {
                    // Replace it by the slice
                    Index ind_I_min_n( factor(*I), factor(*I) / n );
                    Index ind_n( factor(*I), n );
                    Factor slice_I( factor(*I) / n );
                    for( size_t ind_I = 0; ind_I < factor(*I).stateSpace(); ++ind_I, ++ind_I_min_n, ++ind_n )
                        if( ind_n == i )
                            slice_I[ind_I_min_n] = factor(*I)[ind_I];
                    factor(*I) = slice_I;

                    // Remove the edge between n and I
                    _E12.erase( _E12.begin() + VV2E(ni, *I) );
                }
            }

            Regenerate();
            
            // remove all unconnected factors
            for( size_t I = 0; I < nrFactors(); I++ )
                if( nb2(I).size() == 0 )
                    DeleteFactor(I--);

            DeleteVar( ni );

            // FIXME
        }
    } */

    // The cheap solution (at least in terms of coding time) is to multiply every factor
    // that contains the variable with a delta function

    Factor delta_n_i(n,0.0);
    delta_n_i[i] = 1.0;

    // For all factors that contain n
    for( size_t I = 0; I < nrFactors(); I++ ) 
        if( factor(I).vars() && n )
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
        if( factor(I).vars() && ns )
            _undoProbs[I] = factor(I).p();
}


void FactorGraph::undoProbs( const VarSet &ns ) {
    for( map<size_t,Prob>::iterator uI = _undoProbs.begin(); uI != _undoProbs.end(); ) {
        if( factor((*uI).first).vars() && ns ) {
//          cout << "undoing " << factor((*uI).first).vars() << endl;
//          cout << "from " << factor((*uI).first).p() << " to " << (*uI).second << endl;
            factor((*uI).first).p() = (*uI).second;
            _undoProbs.erase(uI++);
        } else
            uI++;
    }
}


bool FactorGraph::isConnected() const {
    if( nrVars() == 0 )
        return false;
    else {
        Var n = var( 0 );

        VarSet component = n;

        VarSet remaining;
        for( size_t i = 1; i < nrVars(); i++ )
            remaining |= var(i);

        bool found_new_vars = true;
        while( found_new_vars ) {
            VarSet new_vars;
            for( VarSet::const_iterator m = remaining.begin(); m != remaining.end(); m++ )
                if( delta(*m) && component )
                    new_vars |= *m;

            if( new_vars.empty() )
                found_new_vars = false;
            else 
                found_new_vars = true;

            component |= new_vars;
            remaining /= new_vars;
        };
        return remaining.empty();
    }
}
