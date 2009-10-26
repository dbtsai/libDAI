/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <iostream>
#include <stack>
#include <dai/jtree.h>


namespace dai {


using namespace std;


const char *JTree::Name = "JTREE";


void JTree::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("verbose") );
    DAI_ASSERT( opts.hasKey("updates") );

    props.verbose = opts.getStringAs<size_t>("verbose");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");
    if( opts.hasKey("inference") )
        props.inference = opts.getStringAs<Properties::InfType>("inference");
    else
        props.inference = Properties::InfType::SUMPROD;
}


PropertySet JTree::getProperties() const {
    PropertySet opts;
    opts.Set( "verbose", props.verbose );
    opts.Set( "updates", props.updates );
    opts.Set( "inference", props.inference );
    return opts;
}


string JTree::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "verbose=" << props.verbose << ",";
    s << "updates=" << props.updates << ",";
    s << "inference=" << props.inference << "]";
    return s.str();
}


JTree::JTree( const FactorGraph &fg, const PropertySet &opts, bool automatic ) : DAIAlgRG(fg), _mes(), _logZ(), RTree(), Qa(), Qb(), props() {
    setProperties( opts );

    if( !isConnected() )
       DAI_THROW(FACTORGRAPH_NOT_CONNECTED);

    if( automatic ) {
        // Create ClusterGraph which contains factors as clusters
        vector<VarSet> cl;
        cl.reserve( fg.nrFactors() );
        for( size_t I = 0; I < nrFactors(); I++ )
            cl.push_back( factor(I).vars() );
        ClusterGraph _cg( cl );

        if( props.verbose >= 3 )
            cerr << "Initial clusters: " << _cg << endl;

        // Retain only maximal clusters
        _cg.eraseNonMaximal();
        if( props.verbose >= 3 )
            cerr << "Maximal clusters: " << _cg << endl;

        vector<VarSet> ElimVec = _cg.VarElim_MinFill().eraseNonMaximal().toVector();
        if( props.verbose >= 3 )
            cerr << "VarElim_MinFill result: " << ElimVec << endl;

        GenerateJT( ElimVec );
    }
}


void JTree::GenerateJT( const std::vector<VarSet> &Cliques ) {
    // Construct a weighted graph (each edge is weighted with the cardinality
    // of the intersection of the nodes, where the nodes are the elements of
    // Cliques).
    WeightedGraph<int> JuncGraph;
    for( size_t i = 0; i < Cliques.size(); i++ )
        for( size_t j = i+1; j < Cliques.size(); j++ ) {
            size_t w = (Cliques[i] & Cliques[j]).size();
            if( w )
                JuncGraph[UEdge(i,j)] = w;
        }

    // Construct maximal spanning tree using Prim's algorithm
    RTree = MaxSpanningTreePrims( JuncGraph );

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
                fac2OR.push_back( alpha );
                break;
            }
        DAI_ASSERT( alpha != nrORs() );
    }
    RecomputeORs();

    // Create inner regions and edges
    IRs.reserve( RTree.size() );
    vector<Edge> edges;
    edges.reserve( 2 * RTree.size() );
    for( size_t i = 0; i < RTree.size(); i++ ) {
        edges.push_back( Edge( RTree[i].n1, nrIRs() ) );
        edges.push_back( Edge( RTree[i].n2, nrIRs() ) );
        // inner clusters have counting number -1
        IRs.push_back( Region( Cliques[RTree[i].n1] & Cliques[RTree[i].n2], -1.0 ) );
    }

    // create bipartite graph
    G.construct( nrORs(), nrIRs(), edges.begin(), edges.end() );

    // Create messages and beliefs
    Qa.clear();
    Qa.reserve( nrORs() );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        Qa.push_back( OR(alpha) );

    Qb.clear();
    Qb.reserve( nrIRs() );
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        Qb.push_back( Factor( IR(beta), 1.0 ) );

    _mes.clear();
    _mes.reserve( nrORs() );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        _mes.push_back( vector<Factor>() );
        _mes[alpha].reserve( nbOR(alpha).size() );
        foreach( const Neighbor &beta, nbOR(alpha) )
            _mes[alpha].push_back( Factor( IR(beta), 1.0 ) );
    }

    // Check counting numbers
    checkCountingNumbers();

    if( props.verbose >= 3 ) {
        cerr << "Resulting regiongraph: " << *this << endl;
    }
}


string JTree::identify() const {
    return string(Name) + printProperties();
}


Factor JTree::belief( const VarSet &ns ) const {
    vector<Factor>::const_iterator beta;
    for( beta = Qb.begin(); beta != Qb.end(); beta++ )
        if( beta->vars() >> ns )
            break;
    if( beta != Qb.end() )
        return( beta->marginal(ns) );
    else {
        vector<Factor>::const_iterator alpha;
        for( alpha = Qa.begin(); alpha != Qa.end(); alpha++ )
            if( alpha->vars() >> ns )
                break;
        DAI_ASSERT( alpha != Qa.end() );
        return( alpha->marginal(ns) );
    }
}


vector<Factor> JTree::beliefs() const {
    vector<Factor> result;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        result.push_back( Qb[beta] );
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        result.push_back( Qa[alpha] );
    return result;
}


Factor JTree::belief( const Var &n ) const {
    return belief( (VarSet)n );
}


// Needs no init
void JTree::runHUGIN() {
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        Qa[alpha] = OR(alpha);

    for( size_t beta = 0; beta < nrIRs(); beta++ )
        Qb[beta].fill( 1.0 );

    // CollectEvidence
    _logZ = 0.0;
    for( size_t i = RTree.size(); (i--) != 0; ) {
//      Make outer region RTree[i].n1 consistent with outer region RTree[i].n2
//      IR(i) = seperator OR(RTree[i].n1) && OR(RTree[i].n2)
        Factor new_Qb;
        if( props.inference == Properties::InfType::SUMPROD )
            new_Qb = Qa[RTree[i].n2].marginal( IR( i ), false );
        else
            new_Qb = Qa[RTree[i].n2].maxMarginal( IR( i ), false );

        _logZ += log(new_Qb.normalize());
        Qa[RTree[i].n1] *= new_Qb / Qb[i];
        Qb[i] = new_Qb;
    }
    if( RTree.empty() )
        _logZ += log(Qa[0].normalize() );
    else
        _logZ += log(Qa[RTree[0].n1].normalize());

    // DistributeEvidence
    for( size_t i = 0; i < RTree.size(); i++ ) {
//      Make outer region RTree[i].n2 consistent with outer region RTree[i].n1
//      IR(i) = seperator OR(RTree[i].n1) && OR(RTree[i].n2)
        Factor new_Qb;
        if( props.inference == Properties::InfType::SUMPROD )
            new_Qb = Qa[RTree[i].n1].marginal( IR( i ) );
        else
            new_Qb = Qa[RTree[i].n1].maxMarginal( IR( i ) );

        Qa[RTree[i].n2] *= new_Qb / Qb[i];
        Qb[i] = new_Qb;
    }

    // Normalize
    for( size_t alpha = 0; alpha < nrORs(); alpha++ )
        Qa[alpha].normalize();
}


// Really needs no init! Initial messages can be anything.
void JTree::runShaferShenoy() {
    // First pass
    _logZ = 0.0;
    for( size_t e = nrIRs(); (e--) != 0; ) {
        // send a message from RTree[e].n2 to RTree[e].n1
        // or, actually, from the seperator IR(e) to RTree[e].n1

        size_t i = nbIR(e)[1].node; // = RTree[e].n2
        size_t j = nbIR(e)[0].node; // = RTree[e].n1
        size_t _e = nbIR(e)[0].dual;

        Factor msg = OR(i);
        foreach( const Neighbor &k, nbOR(i) )
            if( k != e )
                msg *= message( i, k.iter );
        if( props.inference == Properties::InfType::SUMPROD )
            message( j, _e ) = msg.marginal( IR(e), false );
        else
            message( j, _e ) = msg.maxMarginal( IR(e), false );
        _logZ += log( message(j,_e).normalize() );
    }

    // Second pass
    for( size_t e = 0; e < nrIRs(); e++ ) {
        size_t i = nbIR(e)[0].node; // = RTree[e].n1
        size_t j = nbIR(e)[1].node; // = RTree[e].n2
        size_t _e = nbIR(e)[1].dual;

        Factor msg = OR(i);
        foreach( const Neighbor &k, nbOR(i) )
            if( k != e )
                msg *= message( i, k.iter );
        if( props.inference == Properties::InfType::SUMPROD )
            message( j, _e ) = msg.marginal( IR(e) );
        else
            message( j, _e ) = msg.maxMarginal( IR(e) );
    }

    // Calculate beliefs
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        Factor piet = OR(alpha);
        foreach( const Neighbor &k, nbOR(alpha) )
            piet *= message( alpha, k.iter );
        if( nrIRs() == 0 ) {
            _logZ += log( piet.normalize() );
            Qa[alpha] = piet;
        } else if( alpha == nbIR(0)[0].node /*RTree[0].n1*/ ) {
            _logZ += log( piet.normalize() );
            Qa[alpha] = piet;
        } else
            Qa[alpha] = piet.normalized();
    }

    // Only for logZ (and for belief)...
    for( size_t beta = 0; beta < nrIRs(); beta++ ) {
        if( props.inference == Properties::InfType::SUMPROD )
            Qb[beta] = Qa[nbIR(beta)[0].node].marginal( IR(beta) );
        else
            Qb[beta] = Qa[nbIR(beta)[0].node].maxMarginal( IR(beta) );
    }
}


double JTree::run() {
    if( props.updates == Properties::UpdateType::HUGIN )
        runHUGIN();
    else if( props.updates == Properties::UpdateType::SHSH )
        runShaferShenoy();
    return 0.0;
}


Real JTree::logZ() const {
    Real s = 0.0;
    for( size_t beta = 0; beta < nrIRs(); beta++ )
        s += IR(beta).c() * Qb[beta].entropy();
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        s += OR(alpha).c() * Qa[alpha].entropy();
        s += (OR(alpha).log(true) * Qa[alpha]).sum();
    }
    return s;
}



size_t JTree::findEfficientTree( const VarSet& ns, DEdgeVec &Tree, size_t PreviousRoot ) const {
    // find new root clique (the one with maximal statespace overlap with ns)
    size_t maxval = 0, maxalpha = 0;
    for( size_t alpha = 0; alpha < nrORs(); alpha++ ) {
        size_t val = VarSet(ns & OR(alpha).vars()).nrStates();
        if( val > maxval ) {
            maxval = val;
            maxalpha = alpha;
        }
    }

    // grow new tree
    Graph oldTree;
    for( DEdgeVec::const_iterator e = RTree.begin(); e != RTree.end(); e++ )
        oldTree.insert( UEdge(e->n1, e->n2) );
    DEdgeVec newTree = GrowRootedTree( oldTree, maxalpha );

    // identify subtree that contains variables of ns which are not in the new root
    VarSet nsrem = ns / OR(maxalpha).vars();
    set<DEdge> subTree;
    // for each variable in ns that is not in the root clique
    for( VarSet::const_iterator n = nsrem.begin(); n != nsrem.end(); n++ ) {
        // find first occurence of *n in the tree, which is closest to the root
        size_t e = 0;
        for( ; e != newTree.size(); e++ ) {
            if( OR(newTree[e].n2).vars().contains( *n ) )
                break;
        }
        DAI_ASSERT( e != newTree.size() );

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
        DAI_ASSERT( e != newTree.size() );

        // track-back path to root and add edges to subTree
        subTree.insert( newTree[e] );
        size_t pos = newTree[e].n1;
        for( ; e > 0; e-- )
            if( newTree[e-1].n2 == pos ) {
                subTree.insert( newTree[e-1] );
                pos = newTree[e-1].n1;
            }
    }

    // Resulting Tree is a reordered copy of newTree
    // First add edges in subTree to Tree
    Tree.clear();
    for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        if( subTree.count( *e ) ) {
            Tree.push_back( *e );
        }
    // Then add edges pointing away from nsrem
    // FIXME
/*  for( DEdgeVec::const_iterator e = newTree.begin(); e != newTree.end(); e++ )
        for( set<DEdge>::const_iterator sTi = subTree.begin(); sTi != subTree.end(); sTi++ )
            if( *e != *sTi ) {
                if( e->n1 == sTi->n1 || e->n1 == sTi->n2 ||
                    e->n2 == sTi->n1 || e->n2 == sTi->n2 ) {
                    Tree.push_back( *e );
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
            }
        }*/
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
    for( beta = Qb.begin(); beta != Qb.end(); beta++ )
        if( beta->vars() >> ns )
            break;
    if( beta != Qb.end() )
        return( beta->marginal(ns) );
    else {
        vector<Factor>::const_iterator alpha;
        for( alpha = Qa.begin(); alpha != Qa.end(); alpha++ )
            if( alpha->vars() >> ns )
                break;
        if( alpha != Qa.end() )
            return( alpha->marginal(ns) );
        else {
            // Find subtree to do efficient inference
            DEdgeVec T;
            size_t Tsize = findEfficientTree( ns, T );

            // Find remaining variables (which are not in the new root)
            VarSet nsrem = ns / OR(T.front().n1).vars();
            Factor Pns (ns, 0.0);

            // Save Qa and Qb on the subtree
            map<size_t,Factor> Qa_old;
            map<size_t,Factor> Qb_old;
            vector<size_t> b(Tsize, 0);
            for( size_t i = Tsize; (i--) != 0; ) {
                size_t alpha1 = T[i].n1;
                size_t alpha2 = T[i].n2;
                size_t beta;
                for( beta = 0; beta < nrIRs(); beta++ )
                    if( UEdge( RTree[beta].n1, RTree[beta].n2 ) == UEdge( alpha1, alpha2 ) )
                        break;
                DAI_ASSERT( beta != nrIRs() );
                b[i] = beta;

                if( !Qa_old.count( alpha1 ) )
                    Qa_old[alpha1] = Qa[alpha1];
                if( !Qa_old.count( alpha2 ) )
                    Qa_old[alpha2] = Qa[alpha2];
                if( !Qb_old.count( beta ) )
                    Qb_old[beta] = Qb[beta];
            }

            // For all states of nsrem
            for( State s(nsrem); s.valid(); s++ ) {
                // CollectEvidence
                double logZ = 0.0;
                for( size_t i = Tsize; (i--) != 0; ) {
                // Make outer region T[i].n1 consistent with outer region T[i].n2
                // IR(i) = seperator OR(T[i].n1) && OR(T[i].n2)

                    for( VarSet::const_iterator n = nsrem.begin(); n != nsrem.end(); n++ )
                        if( Qa[T[i].n2].vars() >> *n ) {
                            Factor piet( *n, 0.0 );
                            piet[s(*n)] = 1.0;
                            Qa[T[i].n2] *= piet;
                        }

                    Factor new_Qb = Qa[T[i].n2].marginal( IR( b[i] ), false );
                    logZ += log(new_Qb.normalize());
                    Qa[T[i].n1] *= new_Qb / Qb[b[i]];
                    Qb[b[i]] = new_Qb;
                }
                logZ += log(Qa[T[0].n1].normalize());

                Factor piet( nsrem, 0.0 );
                piet[s] = exp(logZ);
                Pns += piet * Qa[T[0].n1].marginal( ns / nsrem, false );      // OPTIMIZE ME

                // Restore clamped beliefs
                for( map<size_t,Factor>::const_iterator alpha = Qa_old.begin(); alpha != Qa_old.end(); alpha++ )
                    Qa[alpha->first] = alpha->second;
                for( map<size_t,Factor>::const_iterator beta = Qb_old.begin(); beta != Qb_old.end(); beta++ )
                    Qb[beta->first] = beta->second;
            }

            return( Pns.normalized() );
        }
    }
}


/// Calculates upper bound to the treewidth of a FactorGraph
/** \relates JTree
 *  \return a pair (number of variables in largest clique, number of states in largest clique)
 */
std::pair<size_t,size_t> treewidth( const FactorGraph & fg ) {
    ClusterGraph _cg;

    // Copy factors
    for( size_t I = 0; I < fg.nrFactors(); I++ )
        _cg.insert( fg.factor(I).vars() );

    // Retain only maximal clusters
    _cg.eraseNonMaximal();

    // Obtain elimination sequence
    vector<VarSet> ElimVec = _cg.VarElim_MinFill().eraseNonMaximal().toVector();

    // Calculate treewidth
    size_t treewidth = 0;
    size_t nrstates = 0;
    for( size_t i = 0; i < ElimVec.size(); i++ ) {
        if( ElimVec[i].size() > treewidth )
            treewidth = ElimVec[i].size();
        size_t s = ElimVec[i].nrStates();
        if( s > nrstates )
            nrstates = s;
    }

    return pair<size_t,size_t>(treewidth, nrstates);
}


std::vector<size_t> JTree::findMaximum() const {
    vector<size_t> maximum( nrVars() );
    vector<bool> visitedVars( nrVars(), false );
    vector<bool> visitedFactors( nrFactors(), false );
    stack<size_t> scheduledFactors;
    for( size_t i = 0; i < nrVars(); ++i ) {
        if( visitedVars[i] )
            continue;
        visitedVars[i] = true;

        // Maximise with respect to variable i
        Prob prod = beliefV(i).p();
        maximum[i] = max_element( prod.begin(), prod.end() ) - prod.begin();

        foreach( const Neighbor &I, nbV(i) )
            if( !visitedFactors[I] )
                scheduledFactors.push(I);

        while( !scheduledFactors.empty() ){
            size_t I = scheduledFactors.top();
            scheduledFactors.pop();
            if( visitedFactors[I] )
                continue;
            visitedFactors[I] = true;

            // Evaluate if some neighboring variables still need to be fixed; if not, we're done
            bool allDetermined = true;
            foreach( const Neighbor &j, nbF(I) )
                if( !visitedVars[j.node] ) {
                    allDetermined = false;
                    break;
                }
            if( allDetermined )
                continue;

            // Calculate product of incoming messages on factor I
            Prob prod2 = beliefF(I).p();

            // The allowed configuration is restrained according to the variables assigned so far:
            // pick the argmax amongst the allowed states
            Real maxProb = numeric_limits<Real>::min();
            State maxState( factor(I).vars() );
            for( State s( factor(I).vars() ); s.valid(); ++s ){
                // First, calculate whether this state is consistent with variables that
                // have been assigned already
                bool allowedState = true;
                foreach( const Neighbor &j, nbF(I) )
                    if( visitedVars[j.node] && maximum[j.node] != s(var(j.node)) ) {
                        allowedState = false;
                        break;
                    }
                // If it is consistent, check if its probability is larger than what we have seen so far
                if( allowedState && prod2[s] > maxProb ) {
                    maxState = s;
                    maxProb = prod2[s];
                }
            }

            // Decode the argmax
            foreach( const Neighbor &j, nbF(I) ) {
                if( visitedVars[j.node] ) {
                    // We have already visited j earlier - hopefully our state is consistent
                    if( maximum[j.node] != maxState(var(j.node)) && props.verbose >= 1 )
                        cerr << "JTree::findMaximum - warning: maximum not consistent due to loops." << endl;
                } else {
                    // We found a consistent state for variable j
                    visitedVars[j.node] = true;
                    maximum[j.node] = maxState( var(j.node) );
                    foreach( const Neighbor &J, nbV(j) )
                        if( !visitedFactors[J] )
                            scheduledFactors.push(J);
                }
            }
        }
    }
    return maximum;
}


} // end of namespace dai
