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


#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <dai/bp.h>
#include <dai/util.h>
#include <dai/properties.h>


namespace dai {


using namespace std;


const char *BP::Name = "BP";


void BP::setProperties( const PropertySet &opts ) {
    assert( opts.hasKey("tol") );
    assert( opts.hasKey("maxiter") );
    assert( opts.hasKey("logdomain") );
    assert( opts.hasKey("updates") );
    
    props.tol = opts.getStringAs<double>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.logdomain = opts.getStringAs<bool>("logdomain");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");

    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<double>("damping");
    else
        props.damping = 0.0;
    if( opts.hasKey("inference") )
        props.inference = opts.getStringAs<Properties::InfType>("inference");
    else
        props.inference = Properties::InfType::SUMPROD;
}


PropertySet BP::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "logdomain", props.logdomain );
    opts.Set( "updates", props.updates );
    opts.Set( "damping", props.damping );
    opts.Set( "inference", props.inference );
    return opts;
}


string BP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "logdomain=" << props.logdomain << ",";
    s << "updates=" << props.updates << ",";
    s << "damping=" << props.damping << ",";
    s << "inference=" << props.inference << "]";
    return s.str();
}


void BP::construct() {
    // create edge properties
    _edges.clear();
    _edges.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i ) {
        _edges.push_back( vector<EdgeProp>() );
        _edges[i].reserve( nbV(i).size() ); 
        foreach( const Neighbor &I, nbV(i) ) {
            EdgeProp newEP;
            newEP.message = Prob( var(i).states() );
            newEP.newMessage = Prob( var(i).states() );

            newEP.index.reserve( factor(I).states() );
            for( IndexFor k( var(i), factor(I).vars() ); k >= 0; ++k )
                newEP.index.push_back( k );

            newEP.residual = 0.0;
            _edges[i].push_back( newEP );
        }
    }
}


void BP::init() {
    double c = props.logdomain ? 0.0 : 1.0;
    for( size_t i = 0; i < nrVars(); ++i ) {
        foreach( const Neighbor &I, nbV(i) ) {
            message( i, I.iter ).fill( c );
            newMessage( i, I.iter ).fill( c );
        }
    }
}


void BP::findMaxResidual( size_t &i, size_t &_I ) {
    i = 0;
    _I = 0;
    double maxres = residual( i, _I );
    for( size_t j = 0; j < nrVars(); ++j )
        foreach( const Neighbor &I, nbV(j) )
            if( residual( j, I.iter ) > maxres ) {
                i = j;
                _I = I.iter;
                maxres = residual( i, _I );
            }
}


void BP::calcNewMessage( size_t i, size_t _I ) {
    // calculate updated message I->i
    size_t I = nbV(i,_I);

    if( 0 == 1 ) {
        /* UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION */
        Factor prod( factor( I ) );
        foreach( const Neighbor &j, nbF(I) )
            if( j != i ) {     // for all j in I \ i
                foreach( const Neighbor &J, nbV(j) )
                    if( J != I ) {     // for all J in nb(j) \ I 
                        prod *= Factor( var(j), message(j, J.iter) );
                    }
            }
        newMessage(i,_I) = prod.marginal( var(i) ).p();
    } else {
        /* OPTIMIZED VERSION */
        Prob prod( factor(I).p() );
        if( props.logdomain ) 
            prod.takeLog();

        // Calculate product of incoming messages and factor I
        foreach( const Neighbor &j, nbF(I) ) {
            if( j != i ) {     // for all j in I \ i
                size_t _I = j.dual;
                // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
                const ind_t &ind = index(j, _I);

                // prod_j will be the product of messages coming into j
                Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 ); 
                foreach( const Neighbor &J, nbV(j) )
                    if( J != I ) { // for all J in nb(j) \ I 
                        if( props.logdomain )
                            prod_j += message( j, J.iter );
                        else
                            prod_j *= message( j, J.iter );
                    }

                // multiply prod with prod_j
                for( size_t r = 0; r < prod.size(); ++r )
                    if( props.logdomain )
                        prod[r] += prod_j[ind[r]];
                    else
                        prod[r] *= prod_j[ind[r]];
            }
        }
        if( props.logdomain ) {
            prod -= prod.maxVal();
            prod.takeExp();
        }

        // Marginalize onto i
        Prob marg( var(i).states(), 0.0 );
        // ind is the precalculated IndexFor(i,I) i.e. to x_I == k corresponds x_i == ind[k]
        const ind_t ind = index(i,_I);
        if( props.inference == Properties::InfType::SUMPROD ) 
            for( size_t r = 0; r < prod.size(); ++r )
                marg[ind[r]] += prod[r];
        else
            for( size_t r = 0; r < prod.size(); ++r )
                if( prod[r] > marg[ind[r]] ) 
                    marg[ind[r]] = prod[r];
        marg.normalize();

        // Store result
        if( props.logdomain )
            newMessage(i,_I) = marg.log();
        else
            newMessage(i,_I) = marg;
    }
}


// BP::run does not check for NANs for performance reasons
// Somehow NaNs do not often occur in BP...
double BP::run() {
    if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";
    if( props.verbose >= 3)
       cout << endl; 

    double tic = toc();
    Diffs diffs(nrVars(), 1.0);
    
    vector<Edge> update_seq;

    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        old_beliefs.push_back( beliefV(i) );

    size_t nredges = nrEdges();

    if( props.updates == Properties::UpdateType::SEQMAX ) {
        // do the first pass
        for( size_t i = 0; i < nrVars(); ++i )
            foreach( const Neighbor &I, nbV(i) ) {
                calcNewMessage( i, I.iter );
                // calculate initial residuals
                residual( i, I.iter ) = dist( newMessage( i, I.iter ), message( i, I.iter ), Prob::DISTLINF );
            }
    } else {
        update_seq.reserve( nredges );
        for( size_t i = 0; i < nrVars(); ++i )
            foreach( const Neighbor &I, nbV(i) )
                update_seq.push_back( Edge( i, I.iter ) );
    }

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    for( _iters=0; _iters < props.maxiter && diffs.maxDiff() > props.tol; ++_iters ) {
        if( props.updates == Properties::UpdateType::SEQMAX ) {
            // Residuals-BP by Koller et al.
            for( size_t t = 0; t < nredges; ++t ) {
                // update the message with the largest residual
                size_t i, _I;
                findMaxResidual( i, _I );
                updateMessage( i, _I );

                // I->i has been updated, which means that residuals for all
                // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
                foreach( const Neighbor &J, nbV(i) ) {
                    if( J.iter != _I ) {
                        foreach( const Neighbor &j, nbF(J) ) {
                            size_t _J = j.dual;
                            if( j != i ) {
                                calcNewMessage( j, _J );
                                residual( j, _J ) = dist( newMessage( j, _J ), message( j, _J ), Prob::DISTLINF );
                            }
                        }
                    }
                }
            }
        } else if( props.updates == Properties::UpdateType::PARALL ) {
            // Parallel updates 
            for( size_t i = 0; i < nrVars(); ++i )
                foreach( const Neighbor &I, nbV(i) )
                    calcNewMessage( i, I.iter );

            for( size_t i = 0; i < nrVars(); ++i )
                foreach( const Neighbor &I, nbV(i) )
                    updateMessage( i, I.iter );
        } else {
            // Sequential updates
            if( props.updates == Properties::UpdateType::SEQRND )
                random_shuffle( update_seq.begin(), update_seq.end() );
            
            foreach( const Edge &e, update_seq ) {
                calcNewMessage( e.first, e.second );
                updateMessage( e.first, e.second );
            }
        }

        // calculate new beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor nb( beliefV(i) );
            diffs.push( dist( nb, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = nb;
        }

        if( props.verbose >= 3 )
            cout << Name << "::run:  maxdiff " << diffs.maxDiff() << " after " << _iters+1 << " passes" << endl;
    }

    if( diffs.maxDiff() > _maxdiff )
        _maxdiff = diffs.maxDiff();

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
                cout << Name << "::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( props.verbose >= 3 )
                cout << Name << "::run:  ";
                cout << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    return diffs.maxDiff();
}


Factor BP::beliefV( size_t i ) const {
    Prob prod( var(i).states(), props.logdomain ? 0.0 : 1.0 ); 
    foreach( const Neighbor &I, nbV(i) )
        if( props.logdomain )
            prod += newMessage( i, I.iter );
        else
            prod *= newMessage( i, I.iter );
    if( props.logdomain ) {
        prod -= prod.maxVal();
        prod.takeExp();
    }

    prod.normalize();
    return( Factor( var(i), prod ) );
}


Factor BP::belief (const Var &n) const {
    return( beliefV( findVar( n ) ) );
}


vector<Factor> BP::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); ++i )
        result.push_back( beliefV(i) );
    for( size_t I = 0; I < nrFactors(); ++I )
        result.push_back( beliefF(I) );
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
        return beliefF(I).marginal(ns);
    }
}


Factor BP::beliefF (size_t I) const {
    if( 0 == 1 ) {
        /*  UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION */

        Factor prod( factor(I) );
        foreach( const Neighbor &j, nbF(I) ) {
            foreach( const Neighbor &J, nbV(j) ) {
                if( J != I )  // for all J in nb(j) \ I
                    prod *= Factor( var(j), newMessage(j, J.iter) );
            }
        }
        return prod.normalized();
    } else {
        /* OPTIMIZED VERSION */
        Prob prod( factor(I).p() );
        if( props.logdomain )
            prod.takeLog();

        foreach( const Neighbor &j, nbF(I) ) {
            size_t _I = j.dual;
            // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
            const ind_t & ind = index(j, _I);

            // prod_j will be the product of messages coming into j
            Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 ); 
            foreach( const Neighbor &J, nbV(j) ) {
                if( J != I ) { // for all J in nb(j) \ I 
                    if( props.logdomain )
                        prod_j += newMessage( j, J.iter );
                    else
                        prod_j *= newMessage( j, J.iter );
                }
            }

            // multiply prod with prod_j
            for( size_t r = 0; r < prod.size(); ++r ) {
                if( props.logdomain )
                    prod[r] += prod_j[ind[r]];
                else
                    prod[r] *= prod_j[ind[r]];
            }
        }

        if( props.logdomain ) {
            prod -= prod.maxVal();
            prod.takeExp();
        }

        Factor result( factor(I).vars(), prod );
        result.normalize();

        return( result );
    }
}


Real BP::logZ() const {
    Real sum = 0.0;
    for(size_t i = 0; i < nrVars(); ++i )
        sum += (1.0 - nbV(i).size()) * beliefV(i).entropy();
    for( size_t I = 0; I < nrFactors(); ++I )
        sum -= dist( beliefF(I), factor(I), Prob::DISTKL );
    return sum;
}


string BP::identify() const { 
    return string(Name) + printProperties();
}


void BP::init( const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); ++n ) {
        size_t ni = findVar( *n );
        foreach( const Neighbor &I, nbV( ni ) )
            message( ni, I.iter ).fill( props.logdomain ? 0.0 : 1.0 );
    }
}


} // end of namespace dai
