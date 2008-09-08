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
#include <dai/jtree.h>


namespace dai {


using namespace std;


const char *JTree::Name = "JTREE";


bool JTree::checkProperties() {
    if (!HasProperty("verbose") )
        return false;
    if( !HasProperty("updates") )
        return false;
    
    ConvertPropertyTo<size_t>("verbose");
    ConvertPropertyTo<UpdateType>("updates");

    return true;
}


JTree::JTree( const FactorGraph &fg, const Properties &opts, bool automatic ) : DAIAlgRG(fg, opts), _RTree(), _Qa(), _Qb(), _mes(), _logZ() {
    assert( checkProperties() );

    if( automatic ) {
        ClusterGraph _cg;

        // Copy factors
        for( size_t I = 0; I < nrFactors(); I++ )
            _cg.insert( factor(I).vars() );
        if( Verbose() >= 3 )
            cout << "Initial clusters: " << _cg << endl;

        // Retain only maximal clusters
        _cg.eraseNonMaximal();
        if( Verbose() >= 3 )
            cout << "Maximal clusters: " << _cg << endl;

        vector<VarSet> ElimVec = _cg.VarElim_MinFill().eraseNonMaximal().toVector();
        if( Verbose() >= 3 ) {
            cout << "VarElim_MinFill result: {" << endl;
            for( size_t i = 0; i < ElimVec.size(); i++ ) {
                if( i != 0 )
                    cout << ", ";
                cout << ElimVec[i];
            }
            cout << "}" << endl;
        }

        GenerateJT( ElimVec );
    }
}


void JTree::GenerateJT( const vector<VarSet> &Cliques ) {
    // Construct a weighted graph (each edge is weighted with the cardinality 
    // of the intersection of the nodes, where the nodes are the elements of
    // Cliques).
    WeightedGraph<int> JuncGraph;
    for( size_t i = 0; i < Cliques.size(); i++ )
        for( size_t j = i+1; j < Cliques.size(); j++ ) {
            size_t w = (Cliques[i] & Cliques[j]).size();
            JuncGraph[UEdge(i,j)] = w;
        }
    
    // Construct maximal spanning tree using Prim's algorithm
    _RTree = MaxSpanningTreePrim( JuncGraph );

    // Construct corresponding region graph

    // Create outer regions
    ORs.reserve( Cliques.size() );
    for( size_t i = 0; i < Cliques.size(); i++ )
        ORs.push_back( FRegion( Factor(Cliques[i], 1.0), 1.0 ) );

    // For each factor, find an outer region that subsumes that factor.
    // Then, multiply the outer region with that factor.
    for( size_t I = 0; I < nrFactors(); I++ ) {
        size_t alpha;
        for( alpha = 0; alpha < nrORs(); alpha++ )
            if( OR(alpha).vars() >> factor(I).vars() ) {
//              OR(alpha) *= factor(I);
                fac2OR.push_back( alpha );
                break;
            }
        assert( alpha != nrORs() );
    }
    RecomputeORs();

    // Create inner regions and edges
    IRs.reserve( _RTree.size() );
    typedef pair<size_t,size_t> Edge;
    vector<Edge> edges;
    edges.reserve( 2 * _RTree.size() );
    for( size_t i = 0; i < _RTree.size(); i++ ) {
        edges.push_back( Edge( _RTree[i].n1, nrIRs() ) );
        edges.push_back( Edge( _RTree[i].n2, nrIRs() ) );
        // inner clusters have counting number -1
        IRs.push_back( Region( Cliques[_RTree[i].n1] & Cliques[_RTree[i].n2], -1.0 ) );
    }

    // create bipartite graph
    G.create( nrORs(), nrIRs(), edges.begin(), edges.end() );

    // Create messages and beliefs
    _Qa.clear();
    _Qa.reserve( nrORs() );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        _Qa.push_back( OR(alpha) );

    _Qb.clear();
    _Qb.reserve( nrIRs() );
    for( size_t beta = 0; beta < nrIRs(); beta++ ) 
        _Qb.push_back( Factor( IR(beta), 1.0 ) );

    _mes.clear();
    _mes.reserve( nrORs() );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        _mes.push_back( vector<Factor>() );
        _mes[alpha].reserve( nbOR(alpha).size() );
        foreach( const Neighbor &beta, nbOR(alpha) )
            _mes[alpha].push_back( Factor( IR(beta), 1.0 ) );
    }

    // Check counting numbers
    Check_Counting_Numbers();

    if( Verbose() >= 3 ) {
        cout << "Resulting regiongraph: " << *this << endl;
    }
}


string JTree::identify() const {
    stringstream result (stringstream::out);
    result << Name << GetProperties();
    return result.str();
}


Factor JTree::belief( const VarSet &ns ) const {
    vector<Factor>::const_iterator beta;
    for( beta = _Qb.begin(); beta != _Qb.end(); beta++ )
        if( beta->vars() >> ns )
            break;
    if( beta != _Qb.end() )
        return( beta->marginal(ns) );
    else {
        vector<Factor>::const_iterator alpha;
        for( alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
            if( alpha->vars() >> ns )
                break;
        assert( alpha != _Qa.end() );
        return( alpha->marginal(ns) );
    }
}


vector<Factor> JTree::beliefs() const {
    vector<Factor> result;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        result.push_back( _Qb[beta] );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        result.push_back( _Qa[alpha] );
    return result;
}


Factor JTree::belief( const Var &n ) const {
    return belief( (VarSet)n );
}


// Needs no init
void JTree::runHUGIN() {
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        _Qa[alpha] = OR(alpha);

    for( size_t beta = 0; beta < nrIRs(); beta++ )
        _Qb[beta].fill( 1.0 );

    // CollectEvidence
    _logZ = 0.0;
    for( size_t i = _RTree.size(); (i--) != 0; ) {
//      Make outer region _RTree[i].n1 consistent with outer region _RTree[i].n2
//      IR(i) = seperator OR(_RTree[i].n1) && OR(_RTree[i].n2)
        Factor new_Qb = _Qa[_RTree[i].n2].part_sum( IR( i ) );
        _logZ += log(new_Qb.normalize( Prob::NORMPROB ));
        _Qa[_RTree[i].n1] *= new_Qb.divided_by( _Qb[i] ); 
        _Qb[i] = new_Qb;
    }
    if( _RTree.empty() )
        _logZ += log(_Qa[0].normalize( Prob::NORMPROB ) );
    else
        _logZ += log(_Qa[_RTree[0].n1].normalize( Prob::NORMPROB ));

    // DistributeEvidence
    for( size_t i = 0; i < _RTree.size(); i++ ) {
//      Make outer region _RTree[i].n2 consistent with outer region _RTree[i].n1
//      IR(i) = seperator OR(_RTree[i].n1) && OR(_RTree[i].n2)
        Factor new_Qb = _Qa[_RTree[i].n1].marginal( IR( i ) );
        _Qa[_RTree[i].n2] *= new_Qb.divided_by( _Qb[i] ); 
        _Qb[i] = new_Qb;
    }

    // Normalize
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        _Qa[alpha].normalize( Prob::NORMPROB );
}


// Really needs no init! Initial messages can be anything.
void JTree::runShaferShenoy() {
    // First pass
    _logZ = 0.0;
    for( size_t e = nrIRs(); (e--) != 0; ) {
        // send a message from _RTree[e].n2 to _RTree[e].n1
        // or, actually, from the seperator IR(e) to _RTree[e].n1

        size_t i = nbIR(e)[1].node; // = _RTree[e].n2
        size_t j = nbIR(e)[0].node; // = _RTree[e].n1
        size_t _e = nbIR(e)[0].dual;
        
        Factor piet = OR(i);
        foreach( const Neighbor &k, nbOR(i) )
            if( k != e ) 
                piet *= message( i, k.iter );
        message( j, _e ) = piet.part_sum( IR(e) );
        _logZ += log( message(j,_e).normalize( Prob::NORMPROB ) );
    }

    // Second pass
    for( size_t e = 0; e < nrIRs(); e++ ) {
        size_t i = nbIR(e)[0].node; // = _RTree[e].n1
        size_t j = nbIR(e)[1].node; // = _RTree[e].n2
        size_t _e = nbIR(e)[1].dual;
        
        Factor piet = OR(i);
        foreach( const Neighbor &k, nbOR(i) )
            if(  k != e )
                piet *= message( i, k.iter );
        message( j, _e ) = piet.marginal( IR(e) );
    }

    // Calculate beliefs
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        Factor piet = OR(alpha);
        foreach( const Neighbor &k, nbOR(alpha) )
            piet *= message( alpha, k.iter );
        if( nrIRs() == 0 ) {
            _logZ += log( piet.normalize( Prob::NORMPROB ) );
            _Qa[alpha] = piet;
        } else if( alpha == nbIR(0)[0].node /*_RTree[0].n1*/ ) {
            _logZ += log( piet.normalize( Prob::NORMPROB ) );
            _Qa[alpha] = piet;
        } else
            _Qa[alpha] = piet.normalized( Prob::NORMPROB );
    }

    // Only for logZ (and for belief)...
    for( size_t beta = 0; beta < nrIRs(); beta++ ) 
        _Qb[beta] = _Qa[nbIR(beta)[0].node].marginal( IR(beta) );
}


double JTree::run() {
    if( Updates() == UpdateType::HUGIN )
        runHUGIN();
    else if( Updates() == UpdateType::SHSH )
        runShaferShenoy();
    return 0.0;
}


Complex JTree::logZ() const {
    Complex sum = 0.0;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        sum += Complex(IR(beta).c()) * _Qb[beta].entropy();
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        sum += Complex(OR(alpha).c()) * _Qa[alpha].entropy();
        sum += (OR(alpha).log0() * _Qa[alpha]).totalSum();
    }
    return sum;
}



size_t JTree::findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot ) const {
    // find new root clique (the one with maximal statespace overlap with ns)
    size_t maxval = 0, maxalpha = 0;
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        size_t val = (ns & OR(alpha).vars()).states();
        if( val > maxval ) {
            maxval = val;
            maxalpha = alpha;
        }
    }

//  for( size_t e = 0; e < _RTree.size(); e++ )
//      cout << OR(_RTree[e].n1).vars() << "->" << OR(_RTree[e].n2).vars() << ",  ";
//  cout << endl;
    // grow new tree
    Graph oldTree;
    for( DEdgeVec::const_iterator e = _RTree.begin(); e != _RTree.end(); e++ )
        oldTree.insert( UEdge(e->n1, e->n2) );
    DEdgeVec newTree = GrowRootedTree( oldTree, maxalpha );
//  cout << ns << ": ";
//  for( size_t e = 0; e < newTree.size(); e++ )
//      cout << OR(newTree[e].n1).vars() << "->" << OR(newTree[e].n2).vars() << ",  ";
//  cout << endl;
    
    // identify subtree that contains variables of ns which are not in the new root
    VarSet nsrem = ns / OR(maxalpha).vars();
//  cout << "nsrem:" << nsrem << endl;
    set<DEdge> subTree;
    // for each variable in ns that is not in the root clique
    for( VarSet::const_iterator n = nsrem.begin(); n != nsrem.end(); n++ ) {
        // find first occurence of *n in the tree, which is closest to the root
        size_t e = 0;
        for( ; e != newTree.size(); e++ ) {
            if( OR(newTree[e].n2).vars() && *n )
                break;
        }
        assert( e != newTree.size() );

        // track-back path to root and add edges to subTree
        subTree.insert( newTree[e] );
        size_t pos = newTree[e].n1;
        for( ; e > 0; e-- )
            if( newTree[e-1].n2 == pos ) {
                subTree.insert( newTree[e-1] );
                pos = newTree[e-1].n1;
            }
    }
    if( PreviousRoot != (size_t)-1 && PreviousRoot != maxalpha) {
        // find first occurence of PreviousRoot in the tree, which is closest to the new root
        size_t e = 0;
        for( ; e != newTree.size(); e++ ) {
            if( newTree[e].n2 == PreviousRoot )
                break;
        }
        assert( e != newTree.size() );

        // track-back path to root and add edges to subTree
        subTree.insert( newTree[e] );
        size_t pos = newTree[e].n1;
        for( ; e > 0; e-- )
            if( newTree[e-1].n2 == pos ) {
                subTree.insert( newTree[e-1] );
                pos = newTree[e-1].n1;
            }
    }
//  cout << "subTree: " << endl;
//  for( set<DEdge>::const_iterator sTi = subTree.begin(); sTi != subTree.end(); sTi++ )
//      cout << OR(sTi->n1).vars() << "->" << OR(sTi->n2).vars() << ",  ";
//  cout << endl;

    // Resulting Tree is a reordered copy of newTree
    // First add edges in subTree to Tree
    Tree.clear();
    for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        if( subTree.count( *e ) ) {
            Tree.push_back( *e );
//          cout << OR(e->n1).vars() << "->" << OR(e->n2).vars() << ",  ";
        }
//  cout << endl;
    // Then add edges pointing away from nsrem
    // FIXME
/*  for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        for( set<DEdge>::const_iterator sTi = subTree.begin(); sTi != subTree.end(); sTi++ )
            if( *e != *sTi ) {
                if( e->n1 == sTi->n1 || e->n1 == sTi->n2 ||
                    e->n2 == sTi->n1 || e->n2 == sTi->n2 ) {
                    Tree.push_back( *e );
//                  cout << OR(e->n1).vars() << "->" << OR(e->n2).vars() << ",  ";
                }
            }*/
    // FIXME
/*  for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        if( find( Tree.begin(), Tree.end(), *e) == Tree.end() ) {
            bool found = false;
            for( VarSet::const_iterator n = nsrem.begin(); n != nsrem.end(); n++ )
                if( (OR(e->n1).vars() && *n) ) {
                    found = true;
                    break;
                }
            if( found ) {
                Tree.push_back( *e );
                cout << OR(e->n1).vars() << "->" << OR(e->n2).vars() << ",  ";
            }
        }
    cout << endl;*/
    size_t subTreeSize = Tree.size();
    // Then add remaining edges
    for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        if( find( Tree.begin(), Tree.end(), *e ) == Tree.end() )
            Tree.push_back( *e );

    return subTreeSize;
}


// Cutset conditioning
// assumes that run() has been called already
Factor JTree::calcMarginal( const VarSet& ns ) {
    vector<Factor>::const_iterator beta;
    for( beta = _Qb.begin(); beta != _Qb.end(); beta++ )
        if( beta->vars() >> ns )
            break;
    if( beta != _Qb.end() )
        return( beta->marginal(ns) );
    else {
        vector<Factor>::const_iterator alpha;
        for( alpha = _Qa.begin(); alpha != _Qa.end(); alpha++ )
            if( alpha->vars() >> ns )
                break;
        if( alpha != _Qa.end() )
            return( alpha->marginal(ns) );
        else {
            // Find subtree to do efficient inference
            DEdgeVec T;
            size_t Tsize = findEfficientTree( ns, T );

            // Find remaining variables (which are not in the new root)
            VarSet nsrem = ns / OR(T.front().n1).vars();
            Factor Pns (ns, 0.0);
            
            multind mi( nsrem );

            // Save _Qa and _Qb on the subtree
            map<size_t,Factor> _Qa_old;
            map<size_t,Factor> _Qb_old;
            vector<size_t> b(Tsize, 0);
            for( size_t i = Tsize; (i--) != 0; ) {
                size_t alpha1 = T[i].n1;
                size_t alpha2 = T[i].n2;
                size_t beta;
                for( beta = 0; beta < nrIRs(); beta++ )
                    if( UEdge( _RTree[beta].n1, _RTree[beta].n2 ) == UEdge( alpha1, alpha2 ) )
                        break;
                assert( beta != nrIRs() );
                b[i] = beta;

                if( !_Qa_old.count( alpha1 ) )
                    _Qa_old[alpha1] = _Qa[alpha1];
                if( !_Qa_old.count( alpha2 ) )
                    _Qa_old[alpha2] = _Qa[alpha2];
                if( !_Qb_old.count( beta ) )
                    _Qb_old[beta] = _Qb[beta];
            }
                
            // For all states of nsrem
            for( size_t j = 0; j < mi.max(); j++ ) {
                vector<size_t> vi = mi.vi( j );
                
                // CollectEvidence
                double logZ = 0.0;
                for( size_t i = Tsize; (i--) != 0; ) {
            //      Make outer region T[i].n1 consistent with outer region T[i].n2
            //      IR(i) = seperator OR(T[i].n1) && OR(T[i].n2)

                    size_t k = 0;
                    for( VarSet::const_iterator n = nsrem.begin(); n != nsrem.end(); n++, k++ )
                        if( _Qa[T[i].n2].vars() >> *n ) {
                            Factor piet( *n, 0.0 );
                            piet[vi[k]] = 1.0;
                            _Qa[T[i].n2] *= piet; 
                        }

                    Factor new_Qb = _Qa[T[i].n2].part_sum( IR( b[i] ) );
                    logZ += log(new_Qb.normalize( Prob::NORMPROB ));
                    _Qa[T[i].n1] *= new_Qb.divided_by( _Qb[b[i]] ); 
                    _Qb[b[i]] = new_Qb;
                }
                logZ += log(_Qa[T[0].n1].normalize( Prob::NORMPROB ));

                Factor piet( nsrem, 0.0 );
                piet[j] = exp(logZ);
                Pns += piet * _Qa[T[0].n1].part_sum( ns / nsrem );      // OPTIMIZE ME

                // Restore clamped beliefs
                for( map<size_t,Factor>::const_iterator alpha = _Qa_old.begin(); alpha != _Qa_old.end(); alpha++ )
                    _Qa[alpha->first] = alpha->second;
                for( map<size_t,Factor>::const_iterator beta = _Qb_old.begin(); beta != _Qb_old.end(); beta++ )
                    _Qb[beta->first] = beta->second;
            }

            return( Pns.normalized(Prob::NORMPROB) );
        }
    }
}


} // end of namespace dai
