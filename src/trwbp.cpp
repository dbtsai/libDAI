/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij  [joris dot mooij at libdai dot org]
 */


#include <dai/trwbp.h>


#define DAI_TRWBP_FAST 1


namespace dai {


using namespace std;


const char *TRWBP::Name = "TRWBP";


string TRWBP::identify() const {
    return string(Name) + printProperties();
}


/* This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour */
Real TRWBP::logZ() const {
    Real sum = 0.0;
    for( size_t I = 0; I < nrFactors(); I++ ) {
        sum += (beliefF(I) * factor(I).log(true)).sum();  // TRWBP
        if( factor(I).vars().size() == 2 )
            sum -= edgeWeight(I) * MutualInfo( beliefF(I) );  // TRWBP
    }
    for( size_t i = 0; i < nrVars(); ++i )
        sum += beliefV(i).entropy();  // TRWBP
    return sum;
}


/* This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour */
void TRWBP::calcNewMessage( size_t i, size_t _I ) {
    // calculate updated message I->i
    size_t I = nbV(i,_I);
    const Var &v_i = var(i);
    const VarSet &v_I = factor(I).vars();
    Real c_I = edgeWeight(I); // TRWBP: c_I (\mu_I in the paper)

    Prob marg;
    if( v_I.size() == 1 ) { // optimization
        marg = factor(I).p();
    } else
        Factor Fprod( factor(I) );
        Prob &prod = Fprod.p();
        if( props.logdomain ) {
            prod.takeLog();
            prod /= c_I;         // TRWBP
        } else
            prod ^= (1.0 / c_I); // TRWBP
    
        // Calculate product of incoming messages and factor I
        foreach( const Neighbor &j, nbF(I) )
            if( j != i ) { // for all j in I \ i
                const Var &v_j = var(j);

                // TRWBP: corresponds to messages n_jI
                // prod_j will be the product of messages coming into j
                Prob prod_j( v_j.states(), props.logdomain ? 0.0 : 1.0 );
                foreach( const Neighbor &J, nbV(j) ) {
                    Real c_J = edgeWeight(J);
                    if( J != I ) { // for all J in nb(j) \ I
                        if( props.logdomain )
                            prod_j += message( j, J.iter ) * c_J;
                        else
                            prod_j *= message( j, J.iter ) ^ c_J;
                    } else { // TRWBP: multiply by m_Ij^(c_I-1)
                        if( props.logdomain )
                            prod_j += message( j, J.iter ) * (c_J - 1.0);
                        else
                            prod_j *= message( j, J.iter ) ^ (c_J - 1.0);
                    }
                }

                // multiply prod with prod_j
                if( !DAI_TRWBP_FAST ) {
                    /* UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION */
                    if( props.logdomain )
                        Fprod += Factor( v_j, prod_j );
                    else
                        Fprod *= Factor( v_j, prod_j );
                } else {
                    /* OPTIMIZED VERSION */
                    // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
                    const ind_t &ind = index(j, _I);
                    for( size_t r = 0; r < prod.size(); ++r )
                        if( props.logdomain )
                            prod[r] += prod_j[ind[r]];
                        else
                            prod[r] *= prod_j[ind[r]];
                }
            }

        if( props.logdomain ) {
            prod -= prod.max();
            prod.takeExp();
        }

        // Marginalize onto i
        if( !DAI_TRWBP_FAST ) {
            /* UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION */
            if( props.inference == Properties::InfType::SUMPROD )
                marg = Fprod.marginal( v_i ).p();
            else
                marg = Fprod.maxMarginal( v_i ).p();
        } else {
            /* OPTIMIZED VERSION */
            marg = Prob( v_i.states(), 0.0 );
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
        }

        // Store result
        if( props.logdomain )
            newMessage(i,_I) = marg.log();
        else
            newMessage(i,_I) = marg;

        // Update the residual if necessary
        if( props.updates == Properties::UpdateType::SEQMAX )
            updateResidual( i, _I , dist( newMessage( i, _I ), message( i, _I ), Prob::DISTLINF ) );
    }
}


/* This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour */
void TRWBP::calcBeliefV( size_t i, Prob &p ) const {
    p = Prob( var(i).states(), props.logdomain ? 0.0 : 1.0 );
    foreach( const Neighbor &I, nbV(i) ) {
        Real c_I = edgeWeight(I);
        if( props.logdomain )
            p += newMessage( i, I.iter ) * c_I;
        else
            p *= newMessage( i, I.iter ) ^ c_I;
    }
}


/* This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour */
void TRWBP::calcBeliefF( size_t I, Prob &p ) const {
    Real c_I = edgeWeight(I); // TRWBP: c_I
    const VarSet &v_I = factor(I).vars();

    Factor Fprod( factor(I) );
    Prob &prod = Fprod.p();

    if( props.logdomain ) {
        prod.takeLog();
        prod /= c_I; // TRWBP
    } else
        prod ^= (1.0 / c_I); // TRWBP

    // Calculate product of incoming messages and factor I
    foreach( const Neighbor &j, nbF(I) ) {
        const Var &v_j = var(j);

        // TRWBP: corresponds to messages n_jI
        // prod_j will be the product of messages coming into j
        Prob prod_j( v_j.states(), props.logdomain ? 0.0 : 1.0 );
        foreach( const Neighbor &J, nbV(j) ) {
            Real c_J = edgeWeight(J);
            if( J != I ) { // for all J in nb(j) \ I
                if( props.logdomain )
                    prod_j += newMessage( j, J.iter ) * c_J;
                else
                    prod_j *= newMessage( j, J.iter ) ^ c_J;
            } else { // TRWBP: multiply by m_Ij^(c_I-1)
                if( props.logdomain )
                    prod_j += newMessage( j, J.iter ) * (c_J - 1.0);
                else
                    prod_j *= newMessage( j, J.iter ) ^ (c_J - 1.0);
            }
        }

        // multiply prod with prod_j
        if( !DAI_TRWBP_FAST ) {
            /* UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION */
            if( props.logdomain )
                Fprod += Factor( v_j, prod_j );
            else
                Fprod *= Factor( v_j, prod_j );
        } else {
            /* OPTIMIZED VERSION */
            size_t _I = j.dual;
            // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
            const ind_t &ind = index(j, _I);

            for( size_t r = 0; r < prod.size(); ++r ) {
                if( props.logdomain )
                    prod[r] += prod_j[ind[r]];
                else
                    prod[r] *= prod_j[ind[r]];
            }
        }
    }
    
    p = prod;
}


void TRWBP::construct() {
    BP::construct();
    _edge_weight.resize( nrFactors(), 1.0 );
}


} // end of namespace dai
