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


#include <vector>
#include <dai/daialg.h>


namespace dai {


using namespace std;


Factor calcMarginal( const InfAlg & obj, const VarSet & ns, bool reInit ) {
    Factor Pns (ns);
    
    InfAlg *clamped = obj.clone();
    if( !reInit )
        clamped->init();

    Real logZ0 = 0.0;
    for( State s(ns); s.valid(); s++ ) {
        // save unclamped factors connected to ns
        clamped->backupFactors( ns );

        // set clamping Factors to delta functions
        for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
            clamped->clamp( *n, s(*n) );
        
        // run DAIAlg, calc logZ, store in Pns
        if( reInit )
            clamped->init();
        else
            clamped->init(ns);
        clamped->run();

        Real Z;
        if( s == 0 ) {
            logZ0 = clamped->logZ();
            Z = 1.0;
        } else {
            // subtract logZ0 to avoid very large numbers
            Z = exp(clamped->logZ() - logZ0);
        }

        Pns[s] = Z;
        
        // restore clamped factors
        clamped->restoreFactors( ns );
    }

    delete clamped;

    return( Pns.normalized() );
}


vector<Factor> calcPairBeliefs( const InfAlg & obj, const VarSet& ns, bool reInit ) {
    // convert ns to vector<VarSet>
    size_t N = ns.size();
    vector<Var> vns;
    vns.reserve( N );
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); n++ )
        vns.push_back( *n );

    vector<Factor> pairbeliefs;
    pairbeliefs.reserve( N * N );
    for( size_t j = 0; j < N; j++ )
        for( size_t k = 0; k < N; k++ )
            if( j == k )
                pairbeliefs.push_back( Factor() );
            else
                pairbeliefs.push_back( Factor( vns[j] | vns[k] ) );

    InfAlg *clamped = obj.clone();
    if( !reInit )
        clamped->init();

    Real logZ0 = 0.0;
    for( size_t j = 0; j < N; j++ ) {
        // clamp Var j to its possible values
        for( size_t j_val = 0; j_val < vns[j].states(); j_val++ ) {
            clamped->clamp( vns[j], j_val, true );
            if( reInit )
                clamped->init();
            else
                clamped->init(ns);
            clamped->run();

            //if( j == 0 )
            //  logZ0 = obj.logZ();
            double Z_xj = 1.0;
            if( j == 0 && j_val == 0 ) {
                logZ0 = clamped->logZ();
            } else {
                // subtract logZ0 to avoid very large numbers
                Z_xj = exp(clamped->logZ() - logZ0);
            }

            for( size_t k = 0; k < N; k++ ) 
                if( k != j ) {
                    Factor b_k = clamped->belief(vns[k]);
                    for( size_t k_val = 0; k_val < vns[k].states(); k_val++ ) 
                        if( vns[j].label() < vns[k].label() )
                            pairbeliefs[j * N + k][j_val + (k_val * vns[j].states())] = Z_xj * b_k[k_val];
                        else
                            pairbeliefs[j * N + k][k_val + (j_val * vns[k].states())] = Z_xj * b_k[k_val];
                }

            // restore clamped factors
            clamped->restoreFactors( ns );
        }
    }
    
    delete clamped;

    // Calculate result by taking the geometric average
    vector<Factor> result;
    result.reserve( N * (N - 1) / 2 );
    for( size_t j = 0; j < N; j++ )
        for( size_t k = j+1; k < N; k++ )
            result.push_back( (pairbeliefs[j * N + k] * pairbeliefs[k * N + j]) ^ 0.5 );

    return result;
}


Factor calcMarginal2ndO( const InfAlg & obj, const VarSet& ns, bool reInit ) {
    // returns a a probability distribution whose 1st order interactions
    // are unspecified, whose 2nd order interactions approximate those of 
    // the marginal on ns, and whose higher order interactions are absent.

    vector<Factor> pairbeliefs = calcPairBeliefs( obj, ns, reInit );

    Factor Pns (ns);
    for( size_t ij = 0; ij < pairbeliefs.size(); ij++ )
        Pns *= pairbeliefs[ij];
    
    return( Pns.normalized() );
}


vector<Factor> calcPairBeliefsNew( const InfAlg & obj, const VarSet& ns, bool reInit ) {
    vector<Factor> result;
    result.reserve( ns.size() * (ns.size() - 1) / 2 );

    InfAlg *clamped = obj.clone();
    if( !reInit )
        clamped->init();

    Real logZ0 = 0.0;
    VarSet::const_iterator nj = ns.begin();
    for( long j = 0; j < (long)ns.size() - 1; j++, nj++ ) {
        size_t k = 0;
        for( VarSet::const_iterator nk = nj; (++nk) != ns.end(); k++ ) {
            Factor pairbelief( *nj | *nk );

            // clamp Vars j and k to their possible values
            for( size_t j_val = 0; j_val < nj->states(); j_val++ ) 
                for( size_t k_val = 0; k_val < nk->states(); k_val++ ) {
                    // save unclamped factors connected to ns
                    clamped->backupFactors( ns );
            
                    clamped->clamp( *nj, j_val );
                    clamped->clamp( *nk, k_val );
                    if( reInit )
                        clamped->init();
                    else
                        clamped->init(ns);
                    clamped->run();

                    double Z_xj = 1.0;
                    if( j_val == 0 && k_val == 0 ) {
                        logZ0 = clamped->logZ();
                    } else {
                        // subtract logZ0 to avoid very large numbers
                        Z_xj = exp(clamped->logZ() - logZ0);
                    }

                    // we assume that j.label() < k.label()
                    // i.e. we make an assumption here about the indexing
                    pairbelief[j_val + (k_val * nj->states())] = Z_xj;

                    // restore clamped factors
                    clamped->restoreFactors( ns );
                }
        
            result.push_back( pairbelief );
        }
    }
    
    delete clamped;

    assert( result.size() == (ns.size() * (ns.size() - 1) / 2) );

    return result;
}


} // end of namespace dai
