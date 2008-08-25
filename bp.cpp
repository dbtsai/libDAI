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
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include "bp.h"
#include "diffs.h"
#include "util.h"
#include "properties.h"


namespace dai {


using namespace std;


const char *BP::Name = "BP";


bool BP::checkProperties() {
    if( !HasProperty("updates") )
        return false;
    if( !HasProperty("tol") )
        return false;
    if (!HasProperty("maxiter") )
        return false;
    if (!HasProperty("verbose") )
        return false;
    
    ConvertPropertyTo<double>("tol");
    ConvertPropertyTo<size_t>("maxiter");
    ConvertPropertyTo<size_t>("verbose");
    ConvertPropertyTo<UpdateType>("updates");

    return true;
}


void BP::Regenerate() {
    DAIAlgFG::Regenerate();
    
    // clear messages
    _messages.clear();
    _messages.reserve(nr_edges());

    // clear indices
    _indices.clear();
    _indices.reserve(nr_edges());

    // create messages and indices
    for( vector<_edge_t>::const_iterator iI=edges().begin(); iI!=edges().end(); iI++ ) {
        _messages.push_back( Prob( var(iI->first).states() ) );

        vector<size_t> ind( factor(iI->second).stateSpace(), 0 );
        Index i (var(iI->first), factor(iI->second).vars() );
        for( size_t j = 0; i >= 0; ++i,++j )
            ind[j] = i; 
        _indices.push_back( ind );
    }

    // create new_messages
    _newmessages = _messages;
}


void BP::init() {
    assert( checkProperties() );
    for( vector<Prob>::iterator mij = _messages.begin(); mij != _messages.end(); mij++ )
        mij->fill(1.0 / mij->size());
    _newmessages = _messages;
}


void BP::calcNewMessage (size_t iI) {
    // calculate updated message I->i
    size_t i = edge(iI).first;
    size_t I = edge(iI).second;

/*  UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION

    Factor prod( factor( I ) );
    for( _nb_cit j = nb2(I).begin(); j != nb2(I).end(); j++ )
        if( *j != i ) {     // for all j in I \ i
            for( _nb_cit J = nb1(*j).begin(); J != nb1(*j).end(); J++ ) 
                if( *J != I ) {     // for all J in nb(j) \ I 
                    prod *= Factor( *j, message(*j,*J) );
    Factor marg = prod.marginal(var(i));
*/
    
    Prob prod( factor(I).p() );

    // Calculate product of incoming messages and factor I
    for( _nb_cit j = nb2(I).begin(); j != nb2(I).end(); j++ )
        if( *j != i ) {     // for all j in I \ i
            // ind is the precalculated Index(j,I) i.e. to x_I == k corresponds x_j == ind[k]
            _ind_t* ind = &(index(*j,I));

            // prod_j will be the product of messages coming into j
            Prob prod_j( var(*j).states() ); 
            for( _nb_cit J = nb1(*j).begin(); J != nb1(*j).end(); J++ ) 
                if( *J != I )   // for all J in nb(j) \ I 
                    prod_j *= message(*j,*J);

            // multiply prod with prod_j
            for( size_t r = 0; r < prod.size(); r++ )
                prod[r] *= prod_j[(*ind)[r]];
        }

    // Marginalize onto i
    Prob marg( var(i).states(), 0.0 );
    // ind is the precalculated Index(i,I) i.e. to x_I == k corresponds x_i == ind[k]
    _ind_t* ind = &(index(i,I));
    for( size_t r = 0; r < prod.size(); r++ )
        marg[(*ind)[r]] += prod[r];
    marg.normalize( _normtype );
    
    // Store result
    _newmessages[iI] = marg;
}


// BP::run does not check for NANs for performance reasons
// Somehow NaNs do not often occur in BP...
double BP::run() {
    if( Verbose() >= 1 )
        cout << "Starting " << identify() << "...";
    if( Verbose() >= 3)
       cout << endl; 

    clock_t tic = toc();
    Diffs diffs(nrVars(), 1.0);
    
    vector<size_t> edge_seq;
    vector<double> residuals;

    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        old_beliefs.push_back(belief1(i));

    size_t iter = 0;

    if( Updates() == UpdateType::SEQMAX ) {
        // do the first pass
        for(size_t iI = 0; iI < nr_edges(); iI++ ) 
            calcNewMessage(iI);

        // calculate initial residuals
        residuals.reserve(nr_edges());
        for( size_t iI = 0; iI < nr_edges(); iI++ )
            residuals.push_back( dist( _newmessages[iI], _messages[iI], Prob::DISTLINF ) );
    } else {
        edge_seq.reserve( nr_edges() );
        for( size_t i = 0; i < nr_edges(); i++ )
            edge_seq.push_back( i );
    }

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( iter=0; iter < MaxIter() && diffs.max() > Tol(); iter++ ) {
        if( Updates() == UpdateType::SEQMAX ) {
            // Residuals-BP by Koller et al.
            for( size_t t = 0; t < nr_edges(); t++ ) {
                // update the message with the largest residual
                size_t iI = max_element(residuals.begin(), residuals.end()) - residuals.begin();
                _messages[iI] = _newmessages[iI];
                residuals[iI] = 0;

                // I->i has been updated, which means that residuals for all
                // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
                size_t i = edge(iI).first;
                size_t I = edge(iI).second;
                for( _nb_cit J = nb1(i).begin(); J != nb1(i).end(); J++ ) 
                    if( *J != I )
                        for( _nb_cit j = nb2(*J).begin(); j != nb2(*J).end(); j++ )
                            if( *j != i ) {
                                size_t jJ = VV2E(*j,*J);
                                calcNewMessage(jJ);
                                residuals[jJ] = dist( _newmessages[jJ], _messages[jJ], Prob::DISTLINF );
                            }
            }
        } else if( Updates() == UpdateType::PARALL ) {
            // Parallel updates 
            for( size_t t = 0; t < nr_edges(); t++ )
                calcNewMessage(t);

            for( size_t t = 0; t < nr_edges(); t++ )
                _messages[t] = _newmessages[t];
        } else {
            // Sequential updates
            if( Updates() == UpdateType::SEQRND )
                random_shuffle( edge_seq.begin(), edge_seq.end() );
            
            for( size_t t = 0; t < nr_edges(); t++ ) {
                size_t k = edge_seq[t];
                calcNewMessage(k);
                _messages[k] = _newmessages[k];
            }
        }

        // calculate new beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); i++ ) {
            Factor nb( belief1(i) );
            diffs.push( dist( nb, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = nb;
        }

        if( Verbose() >= 3 )
            cout << "BP::run:  maxdiff " << diffs.max() << " after " << iter+1 << " passes" << endl;
    }

    updateMaxDiff( diffs.max() );

    if( Verbose() >= 1 ) {
        if( diffs.max() > Tol() ) {
            if( Verbose() == 1 )
                cout << endl;
                cout << "BP::run:  WARNING: not converged within " << MaxIter() << " passes (" << toc() - tic << " clocks)...final maxdiff:" << diffs.max() << endl;
        } else {
            if( Verbose() >= 3 )
                cout << "BP::run:  ";
                cout << "converged in " << iter << " passes (" << toc() - tic << " clocks)." << endl;
        }
    }

    return diffs.max();
}


Factor BP::belief1( size_t i ) const {
    Prob prod( var(i).states() ); 
    for( _nb_cit I = nb1(i).begin(); I != nb1(i).end(); I++ ) 
        prod *= newMessage(i,*I);

    prod.normalize( Prob::NORMPROB );
    return( Factor( var(i), prod ) );
}


Factor BP::belief (const Var &n) const {
    return( belief1( findVar( n ) ) );
}


vector<Factor> BP::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( belief1(i) );
    for( size_t I = 0; I < nrFactors(); I++ )
        result.push_back( belief2(I) );
    return result;
}


Factor BP::belief( const VarSet &ns ) const {
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


Factor BP::belief2 (size_t I) const {
    Prob prod( factor(I).p() );

    for( _nb_cit j = nb2(I).begin(); j != nb2(I).end(); j++ ) {
        // ind is the precalculated Index(j,I) i.e. to x_I == k corresponds x_j == ind[k]
        const _ind_t *ind = &(index(*j, I));

        // prod_j will be the product of messages coming into j
        Prob prod_j( var(*j).states() ); 
        for( _nb_cit J = nb1(*j).begin(); J != nb1(*j).end(); J++ ) 
            if( *J != I )   // for all J in nb(j) \ I 
                prod_j *= newMessage(*j,*J);

        // multiply prod with prod_j
        for( size_t r = 0; r < prod.size(); r++ )
            prod[r] *= prod_j[(*ind)[r]];
    }

    Factor result( factor(I).vars(), prod );
    result.normalize( Prob::NORMPROB );
    
    return( result );

/*  UNOPTIMIZED VERSION
 
    Factor prod( factor(I) );
    for( _nb_cit i = nb2(I).begin(); i != nb2(I).end(); i++ ) {
        for( _nb_cit J = nb1(*i).begin(); J != nb1(*i).end(); J++ )
            if( *J != I )
                prod *= Factor( var(*i), newMessage(*i,*J)) );
    }
    return prod.normalize( Prob::NORMPROB );*/
}


Complex BP::logZ() const {
    Complex sum = 0.0;
    for(size_t i = 0; i < nrVars(); i++ )
        sum += Complex(1.0 - nb1(i).size()) * belief1(i).entropy();
    for( size_t I = 0; I < nrFactors(); I++ )
        sum -= KL_dist( belief2(I), factor(I) );
    return sum;
}


string BP::identify() const { 
    stringstream result (stringstream::out);
    result << Name << GetProperties();
    return result.str();
}


void BP::init( const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ ) {
        size_t ni = findVar( *n );
        for( _nb_cit I = nb1(ni).begin(); I != nb1(ni).end(); I++ )
            message(ni,*I).fill( 1.0 );
    }
}


}
