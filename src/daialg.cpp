/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <vector>
#include <dai/daialg.h>


namespace dai {


using namespace std;


Factor calcMarginal( const InfAlg &obj, const VarSet &vs, bool reInit ) {
    Factor Pvs (vs);

    InfAlg *clamped = obj.clone();
    if( !reInit )
        clamped->init();

    map<Var,size_t> varindices;
    for( VarSet::const_iterator n = vs.begin(); n != vs.end(); n++ )
        varindices[*n] = obj.fg().findVar( *n );

    Real logZ0 = -INFINITY;
    for( State s(vs); s.valid(); s++ ) {
        // save unclamped factors connected to vs
        clamped->backupFactors( vs );

        // set clamping Factors to delta functions
        for( VarSet::const_iterator n = vs.begin(); n != vs.end(); n++ )
            clamped->clamp( varindices[*n], s(*n) );

        // run DAIAlg, calc logZ, store in Pvs
        if( reInit )
            clamped->init();
        else
            clamped->init(vs);

        Real logZ;
        try {
            clamped->run();
            logZ = clamped->logZ();
        } catch( Exception &e ) {
            if( e.code() == Exception::NOT_NORMALIZABLE )
                logZ = -INFINITY;
            else
                throw;
        }

        if( logZ0 == -INFINITY )
            if( logZ != -INFINITY )
                logZ0 = logZ;

        if( logZ == -INFINITY )
            Pvs.set( s, 0 );
        else
            Pvs.set( s, exp(logZ - logZ0) ); // subtract logZ0 to avoid very large numbers

        // restore clamped factors
        clamped->restoreFactors( vs );
    }

    delete clamped;

    return( Pvs.normalized() );
}


vector<Factor> calcPairBeliefs( const InfAlg & obj, const VarSet& vs, bool reInit, bool accurate ) {
    vector<Factor> result;
    size_t N = vs.size();
    result.reserve( N * (N - 1) / 2 );

    InfAlg *clamped = obj.clone();
    if( !reInit )
        clamped->init();

    map<Var,size_t> varindices;
    for( VarSet::const_iterator v = vs.begin(); v != vs.end(); v++ )
        varindices[*v] = obj.fg().findVar( *v );

    if( accurate ) {
        Real logZ0 = 0.0;
        VarSet::const_iterator nj = vs.begin();
        for( long j = 0; j < (long)N - 1; j++, nj++ ) {
            size_t k = 0;
            for( VarSet::const_iterator nk = nj; (++nk) != vs.end(); k++ ) {
                Factor pairbelief( VarSet(*nj, *nk) );

                // clamp Vars j and k to their possible values
                for( size_t j_val = 0; j_val < nj->states(); j_val++ )
                    for( size_t k_val = 0; k_val < nk->states(); k_val++ ) {
                        // save unclamped factors connected to vs
                        clamped->backupFactors( vs );

                        clamped->clamp( varindices[*nj], j_val );
                        clamped->clamp( varindices[*nk], k_val );
                        if( reInit )
                            clamped->init();
                        else
                            clamped->init(vs);

                        Real logZ;
                        try {
                            clamped->run();
                            logZ = clamped->logZ();
                        } catch( Exception &e ) {
                            if( e.code() == Exception::NOT_NORMALIZABLE )
                                logZ = -INFINITY;
                            else
                                throw;
                        }

                        if( logZ0 == -INFINITY )
                            if( logZ != -INFINITY )
                                logZ0 = logZ;

                        Real Z_xj;
                        if( logZ == -INFINITY )
                            Z_xj = 0;
                        else
                            Z_xj = exp(logZ - logZ0); // subtract logZ0 to avoid very large numbers

                        // we assume that j.label() < k.label()
                        // i.e. we make an assumption here about the indexing
                        pairbelief.set( j_val + (k_val * nj->states()), Z_xj );

                        // restore clamped factors
                        clamped->restoreFactors( vs );
                    }

                result.push_back( pairbelief.normalized() );
            }
        }
    } else {
        // convert vs to vector<VarSet>
        vector<Var> vvs( vs.begin(), vs.end() );

        vector<Factor> pairbeliefs;
        pairbeliefs.reserve( N * N );
        for( size_t j = 0; j < N; j++ )
            for( size_t k = 0; k < N; k++ )
                if( j == k )
                    pairbeliefs.push_back( Factor() );
                else
                    pairbeliefs.push_back( Factor( VarSet(vvs[j], vvs[k]) ) );

        Real logZ0 = -INFINITY;
        for( size_t j = 0; j < N; j++ ) {
            // clamp Var j to its possible values
            for( size_t j_val = 0; j_val < vvs[j].states(); j_val++ ) {
                clamped->clamp( varindices[vvs[j]], j_val, true );
                if( reInit )
                    clamped->init();
                else
                    clamped->init(vs);

                Real logZ;
                try {
                    clamped->run();
                    logZ = clamped->logZ();
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_NORMALIZABLE )
                        logZ = -INFINITY;
                    else
                        throw;
                }

                if( logZ0 == -INFINITY )
                    if( logZ != -INFINITY )
                        logZ0 = logZ;

                Real Z_xj;
                if( logZ == -INFINITY )
                    Z_xj = 0;
                else
                    Z_xj = exp(logZ - logZ0); // subtract logZ0 to avoid very large numbers

                for( size_t k = 0; k < N; k++ )
                    if( k != j ) {
                        Factor b_k = clamped->belief(vvs[k]);
                        for( size_t k_val = 0; k_val < vvs[k].states(); k_val++ )
                            if( vvs[j].label() < vvs[k].label() )
                                pairbeliefs[j * N + k].set( j_val + (k_val * vvs[j].states()), Z_xj * b_k[k_val] );
                            else
                                pairbeliefs[j * N + k].set( k_val + (j_val * vvs[k].states()), Z_xj * b_k[k_val] );
                    }

                // restore clamped factors
                clamped->restoreFactors( vs );
            }
        }

        // Calculate result by taking the geometric average
        for( size_t j = 0; j < N; j++ )
            for( size_t k = j+1; k < N; k++ )
                result.push_back( ((pairbeliefs[j * N + k] * pairbeliefs[k * N + j]) ^ 0.5).normalized() );
    }
    delete clamped;
    return result;
}


} // end of namespace dai
