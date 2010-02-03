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


// This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour
Real TRWBP::logZ() const {
    Real sum = 0.0;
    for( size_t I = 0; I < nrFactors(); I++ ) {
        sum += (beliefF(I) * factor(I).log(true)).sum();  // TRWBP/FBP
        sum += Weight(I) * beliefF(I).entropy();  // TRWBP/FBP
    }
    for( size_t i = 0; i < nrVars(); ++i ) {
        Real c_i = 0.0;
        foreach( const Neighbor &I, nbV(i) )
            c_i += Weight(I);
        sum += (1.0 - c_i) * beliefV(i).entropy();  // TRWBP/FBP
    }
    return sum;
}


// This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour
Prob TRWBP::calcIncomingMessageProduct( size_t I, bool without_i, size_t i ) const {
    Real c_I = Weight(I); // TRWBP: c_I

    Factor Fprod( factor(I) );
    Prob &prod = Fprod.p();
    if( props.logdomain ) {
        prod.takeLog();
        prod /= c_I; // TRWBP
    } else
        prod ^= (1.0 / c_I); // TRWBP

    // Calculate product of incoming messages and factor I
    foreach( const Neighbor &j, nbF(I) )
        if( !(without_i && (j == i)) ) {
            const Var &v_j = var(j);
            // prod_j will be the product of messages coming into j
            // TRWBP: corresponds to messages n_jI
            Prob prod_j( v_j.states(), props.logdomain ? 0.0 : 1.0 );
            foreach( const Neighbor &J, nbV(j) ) {
                Real c_J = Weight(J);  // TRWBP
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
                // UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION
                if( props.logdomain )
                    Fprod += Factor( v_j, prod_j );
                else
                    Fprod *= Factor( v_j, prod_j );
            } else {
                // OPTIMIZED VERSION
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
    
    return prod;
}


// This code has been copied from bp.cpp, except where comments indicate TRWBP-specific behaviour
void TRWBP::calcBeliefV( size_t i, Prob &p ) const {
    p = Prob( var(i).states(), props.logdomain ? 0.0 : 1.0 );
    foreach( const Neighbor &I, nbV(i) ) {
        Real c_I = Weight(I);
        if( props.logdomain )
            p += newMessage( i, I.iter ) * c_I;
        else
            p *= newMessage( i, I.iter ) ^ c_I;
    }
}


void TRWBP::construct() {
    BP::construct();
    _weight.resize( nrFactors(), 1.0 );
}


} // end of namespace dai
