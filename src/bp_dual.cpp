
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>

#include <dai/bbp.h>
//#include <dai/diffs.h>
#include <dai/util.h>
//#include "stlutil.h"
#include <dai/properties.h>

namespace dai {

using namespace std;

const char *BP_dual::Name = "BP_dual";

const char *BP_dual::PropertyList[] = {"tol","maxiter","updates","verbose"};

void BP_dual::setProperties( const PropertySet &opts ) {
//     DAI_DMSG("in BP_dual::setProperties");
//     DAI_PV(opts);

    bool die=false;
    foreach(const char *p, PropertyList) {
        if( !opts.hasKey(p) ) {
            cerr << "BP_dual: missing property " << p << endl;
            die=true;
        }
    }
    if(die) 
        DAI_THROW(NOT_ALL_PROPERTIES_SPECIFIED);
    
    props.tol = opts.getStringAs<double>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");
    props.verbose = opts.getStringAs<size_t>("verbose");

//     DAI_PV(printProperties());
}

PropertySet BP_dual::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "updates", props.updates );
    opts.Set( "verbose", props.verbose );
    return opts;
}

std::string BP_dual::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "updates=" << props.updates << ",";
    s << "verbose=" << props.verbose;
    s << "]";
    return s.str();
}

// void BP_dual::checkProperties() {
//     const char *props[] = {"updates","tol","maxiter","verbose"};
//     for(size_t i=0; i<sizeof(props)/sizeof(*props); i++) {
//         if(!HasProperty(props[i]))
//             die("BP_dual: Missing property \"%s\"", props[i]);
//     }
    
//     ConvertPropertyTo<double>("tol");
//     ConvertPropertyTo<size_t>("maxiter");
//     ConvertPropertyTo<size_t>("verbose");
//     ConvertPropertyTo<UpdateType>("updates");
// }

void BP_dual::RegenerateIndices() {
    _indices.clear();
    _indices.reserve(nr_edges());

    foreach(Edge iI, edges()) {
        vector<size_t> ind( factor(iI.second).states(), 0 );
        IndexFor i (var(iI.first), factor(iI.second).vars() );
        for( size_t j = 0; i >= 0; ++i,++j )
            ind[j] = i; 
        _indices.push_back( ind );
    }
}

void BP_dual::RegenerateMessages() {
    _msgs.Zn.resize(nr_edges(),1.0);
    _msgs.Zm.resize(nr_edges(),1.0);

    // clear messages
    _msgs.m.clear();
    _msgs.m.reserve(nr_edges());
    _msgs.n.clear();
    _msgs.n.reserve(nr_edges());

    // create messages and indices
    foreach(Edge iI, edges()) {
        // initialize to uniform distributions
        _msgs.m.push_back( Prob( var(iI.first).states() ) );
        _msgs.n.push_back( Prob( var(iI.first).states() ) );
    }

    // create new_messages
    _new_msgs = _msgs;
}

void BP_dual::RegenerateBeliefs() {
    _beliefs.b1.clear();
    _beliefs.b1.reserve(nrVars());
    _beliefs.Zb1.resize(nrVars(), 1.0);
    _beliefs.b2.clear();
    _beliefs.b2.reserve(nrFactors());
    _beliefs.Zb2.resize(nrFactors(), 1.0);

    for(size_t i=0; i<nrVars(); i++) {
        _beliefs.b1.push_back( Prob( var(i).states() ).setUniform() );
    }
    for(size_t I=0; I<nrFactors(); I++) {
        _beliefs.b2.push_back( Prob( factor(I).states() ).setUniform() );
    }
}

// called by constructor, called before 'init'
void BP_dual::Regenerate() {

    indexEdges(); // so we can use compatibility interface

//     DAIAlgFG::Regenerate(); // located in BipartiteGraph
  
    RegenerateIndices();
    RegenerateMessages();
    RegenerateBeliefs();

    _maxdiff = 0;
    _iters = 0;
}

void BP_dual::CalcBelief1(size_t i) {
    Prob prod( var(i).states(), 1.0 );
    foreach(size_t I, nbV(i)) {
        prod *= newMsgM(I,i);
    }
    _beliefs.Zb1[i] = prod.normalize();
    _beliefs.b1[i] = prod;
}

void BP_dual::CalcBelief2(size_t I) {
    Prob prod( factor(I).p() );
    foreach(size_t j, nbF(I)) {
        const _ind_t *ind = &(index(j, I));
        for(size_t r=0; r<prod.size(); r++) {
            Prob n(newMsgN(j,I));
            prod[r] *= n[(*ind)[r]];
        }
    }
    _beliefs.Zb2[I] = prod.normalize();
    _beliefs.b2[I] = prod;
}

// called after run()
void BP_dual::CalcBeliefs() {
    for(size_t i=0; i<nrVars(); i++) {
        // calculate b_i
        CalcBelief1(i);
    }
    for(size_t I=0; I<nrFactors(); I++) {
        // calculate b_I
        CalcBelief2(I);
    }
}

void BP_dual::calcNewM(size_t iI) {
    // calculate updated message I->i
    size_t i = edge(iI).first;
    size_t I = edge(iI).second;

    Prob prod( factor(I).p() );

    foreach(size_t j, nbF(I)) {
        if( j != i ) {     // for all j in I \ i
            _ind_t* ind = &(index(j,I));
            Prob n(msgN(j,I));
            for( size_t r = 0; r < prod.size(); r++ )
                prod[r] *= n[(*ind)[r]];
        }
    }

    // Marginalize onto i
    Prob marg( var(i).states(), 0.0 );
    // ind is the precalculated Index(i,I) i.e. to x_I == k corresponds x_i == ind[k]
    _ind_t* ind = &(index(i,I));
    for( size_t r = 0; r < prod.size(); r++ )
        marg[(*ind)[r]] += prod[r];
    
    _new_msgs.Zm[iI] = marg.normalize();
    _new_msgs.m[iI] = marg;
}

void BP_dual::calcNewN(size_t iI) {
    // XXX optimize
    // calculate updated message i->I
    size_t i = edge(iI).first;
    size_t I = edge(iI).second;

    Prob prod(var(i).states(), 1.0);
    foreach(size_t J, nbV(i)) {
        if(J != I) { // for all J in i \ I
            prod *= msgM(J,i);
        }
    }
    _new_msgs.Zn[iI] = prod.normalize();
    _new_msgs.n[iI] = prod;
}

void BP_dual::upMsgM(size_t iI) {
    _msgs.m[iI] = _new_msgs.m[iI];
    _msgs.Zm[iI] = _new_msgs.Zm[iI];
}

void BP_dual::upMsgN(size_t iI) {
    _msgs.n[iI] = _new_msgs.n[iI];
    _msgs.Zn[iI] = _new_msgs.Zn[iI];
}

double BP_dual::run() {
    DAI_IFVERB(1, "Starting " << identify() << "..." << endl);

    double tic = toc();
    // for some reason we need 2* here, where orig BP doesn't
    Diffs diffs(2*nrVars(), 1.0);
    
    vector<size_t> edge_seq;
    vector<double> residuals;

    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ ) {
        CalcBelief1(i);
        old_beliefs.push_back(belief1(i));
    }

    size_t iter = 0;

    if( Updates() == UpdateType::SEQMAX ) {
        // do the first pass
        for(size_t iI = 0; iI < nr_edges(); iI++ ) {
            calcNewM(iI);
            calcNewN(iI);
        }

        // calculate initial residuals
        residuals.reserve(nr_edges());
        for( size_t iI = 0; iI < nr_edges(); iI++ )
            residuals.push_back( dist( _new_msgs.m[iI], _msgs.m[iI], Prob::DISTLINF ) );
    } else {
        edge_seq.reserve( nr_edges() );
        for( size_t i = 0; i < nr_edges(); i++ )
            edge_seq.push_back( i );
    }

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( iter=0; iter < props.maxiter && diffs.maxDiff() > props.tol; iter++ ) {
        if( Updates() == UpdateType::SEQMAX ) {
            // Residuals-BP by Koller et al.
            for( size_t t = 0; t < nr_edges(); t++ ) {
                // update the message with the largest residual
                size_t iI = max_element(residuals.begin(), residuals.end()) - residuals.begin();
                upMsgM(iI);
                residuals[iI] = 0;

                // I->i has been updated, which means that residuals for all
                // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
                size_t i = edge(iI).first;
                size_t I = edge(iI).second;
                foreach(size_t J, nbV(i)) {
                    if(J != I) {
                        size_t iJ = VV2E(i,J);
                        calcNewN(iJ);
                        upMsgN(iJ);
                        foreach(size_t j, nbF(J)) {
                            if(j != i) {
                                size_t jJ = VV2E(j,J);
                                calcNewM(jJ);
                                residuals[jJ] = dist( _new_msgs.m[jJ], _msgs.m[jJ], Prob::DISTLINF );
                            }
                        }
                    }
                }
            }
        } else if( Updates() == UpdateType::PARALL ) {
            // Parallel updates 
            for( size_t t = 0; t < nr_edges(); t++ ) {
                calcNewM(t);
                calcNewN(t);
            }
            if(0) {
                for(size_t t=0; t<nr_edges(); t++) {
                    upMsgM(t); upMsgN(t);
                }
            } else {
                _msgs = _new_msgs;
            }
        } else {
            // Sequential updates
            if( Updates() == UpdateType::SEQRND )
                random_shuffle( edge_seq.begin(), edge_seq.end() );

            foreach(size_t k, edge_seq) {
                calcNewM(k);
                calcNewN(k);
                upMsgM(k);
                upMsgN(k);
            }
        }

        // calculate new beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            CalcBelief1(i);
            Factor nb( belief1(i) );
            diffs.push( dist( nb, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = nb;
        }

        DAI_IFVERB(3,"BP_dual::run:  maxdiff " << diffs.maxDiff() << " after " << iter+1 << " passes" << endl);

        _iters++;
    }

    updateMaxDiff( diffs.maxDiff() );

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            DAI_IFVERB(1, endl << "BP_dual::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.maxDiff() << endl);
        } else {
            DAI_IFVERB(3, "BP_dual::run:  converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl);
        }
    }

    CalcBeliefs();
    
    return diffs.maxDiff();
}

string BP_dual::identify() const { 
    return string(Name) + printProperties();
}

vector<Factor> BP_dual::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( belief1(i) );
    for( size_t I = 0; I < nrFactors(); I++ )
        result.push_back( belief2(I) );
    return result;
}

Factor BP_dual::belief( const VarSet &ns ) const {
    if( ns.size() == 1 )
        return belief( *(ns.begin()) );
    else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        assert( I != nrFactors() );
        return belief2(I).marginal(ns);
    }
}

Real BP_dual::logZ() const {
    Real sum = 0.0;
    for(size_t i = 0; i < nrVars(); i++ )
        sum += Real(1.0 - nbV(i).size()) * belief1(i).entropy();
    for( size_t I = 0; I < nrFactors(); I++ )
        sum -= dist( belief2(I), factor(I), Prob::DISTKL );
    return sum;
}

// reset only messages to/from certain variables
void BP_dual::init(const VarSet &ns) {
    _iters=0;
    foreach(Var n, ns) {
        size_t ni = findVar(n);
        size_t st = n.states();
        foreach(Neighbor I, nbV(ni)) {
            msgM(I.node,ni).fill(1.0/st);
            zM(I.node,ni) = 1.0;
            msgN(ni,I.node).fill(1.0/st);
            zN(ni,I.node) = 1.0;
        }
    }
}

void BP_dual::init() {
    _iters=0;
    for(size_t iI = 0; iI < nr_edges(); iI++ ) {
        _msgs.m[iI].setUniform();
        _msgs.Zm[iI] = 1;
        _msgs.n[iI].setUniform();
        _msgs.Zn[iI] = 1;
    }
    _new_msgs = _msgs;
}


void BP_dual::init(const vector<size_t>& state) {
    _iters=0;
    for(size_t iI = 0; iI < nr_edges(); iI++ ) {
        size_t i = edge(iI).first;
        _msgs.m[iI].fill(0.1);
        _msgs.m[iI][state[i]]=1;
        _msgs.Zm[iI] = _msgs.m[iI].normalize();
        _msgs.n[iI].fill(0.1);
        _msgs.n[iI][state[i]]=1;
        _msgs.Zn[iI] = _msgs.n[iI].normalize();
    }
    _new_msgs = _msgs;
}

void _clamp(FactorGraph &g, const Var & n, const vector<size_t> &is ) {
    Factor mask_n(n,0.0);

    foreach(size_t i, is) { assert( i <= n.states() ); mask_n[i] = 1.0; }

    for( size_t I = 0; I < g.nrFactors(); I++ ) 
        if( g.factor(I).vars().contains( n ) )
          g.factor(I) *= mask_n;
}

// clamp a factor to have one of a set of values
void _clampFactor(FactorGraph &g, size_t I, const vector<size_t> &is) {
    size_t st = g.factor(I).states();
    Prob mask_n(st,0.0);

    foreach(size_t i, is) { assert( i <= st ); mask_n[i] = 1.0; }

    g.factor(I).p() *= mask_n;
}

} // end of namespace dai
