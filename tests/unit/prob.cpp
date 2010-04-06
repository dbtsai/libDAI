/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/prob.h>
#include <strstream>


using namespace dai;


const double tol = 1e-8;


#define BOOST_TEST_MODULE ProbTest


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    Prob x1;
    BOOST_CHECK_EQUAL( x1.size(), 0 );
    BOOST_CHECK( x1.p() == Prob::container_type() );

    Prob x2( 3 );
    BOOST_CHECK_EQUAL( x2.size(), 3 );
    BOOST_CHECK( x2.p() == Prob::container_type( 3, 1.0 / 3.0 ) );
    BOOST_CHECK_EQUAL( x2[0], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x2[1], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x2[2], 1.0 / 3.0 );

    Prob x3( 4, 1.0 );
    BOOST_CHECK_EQUAL( x3.size(), 4 );
    BOOST_CHECK( x3.p() == Prob::container_type( 4, 1.0 ) );
    BOOST_CHECK_EQUAL( x3[0], 1.0 );
    BOOST_CHECK_EQUAL( x3[1], 1.0 );
    BOOST_CHECK_EQUAL( x3[2], 1.0 );
    BOOST_CHECK_EQUAL( x3[3], 1.0 );
    x3.set( 0, 0.5 );
    x3.set( 1, 1.0 );
    x3.set( 2, 2.0 );
    x3.set( 3, 4.0 );

    std::vector<Real> v;
    v.push_back( 0.5 );
    v.push_back( 1.0 );
    v.push_back( 2.0 );
    v.push_back( 4.0 );
    Prob x4( v.begin(), v.end(), 0 );
    BOOST_CHECK_EQUAL( x4.size(), 4 );
    BOOST_CHECK( x4.p() == x3.p() );
    BOOST_CHECK( x4 == x3 );
    BOOST_CHECK_EQUAL( x4[0], 0.5 );
    BOOST_CHECK_EQUAL( x4[1], 1.0 );
    BOOST_CHECK_EQUAL( x4[2], 2.0 );
    BOOST_CHECK_EQUAL( x4[3], 4.0 );

    Prob x5( v.begin(), v.end(), v.size() );
    BOOST_CHECK_EQUAL( x5.size(), 4 );
    BOOST_CHECK( x5.p() == x3.p() );
    BOOST_CHECK( x5 == x3 );
    BOOST_CHECK_EQUAL( x5[0], 0.5 );
    BOOST_CHECK_EQUAL( x5[1], 1.0 );
    BOOST_CHECK_EQUAL( x5[2], 2.0 );
    BOOST_CHECK_EQUAL( x5[3], 4.0 );

    std::vector<int> y( 3, 2 );
    Prob x6( y );
    BOOST_CHECK_EQUAL( x6.size(), 3 );
    BOOST_CHECK_EQUAL( x6[0], 2.0 );
    BOOST_CHECK_EQUAL( x6[1], 2.0 );
    BOOST_CHECK_EQUAL( x6[2], 2.0 );

    Prob x7( x6 );
    BOOST_CHECK( x7 == x6 );
    
    Prob x8 = x6;
    BOOST_CHECK( x8 == x6 );
}


#ifndef DAI_SPARSE
BOOST_AUTO_TEST_CASE( IteratorTest ) {
    Prob x( 5, 0.0 );
    size_t i;
    for( i = 0; i < x.size(); i++ )
        x.set( i, i );

    i = 0;
    for( Prob::const_iterator cit = x.begin(); cit != x.end(); cit++, i++ )
        BOOST_CHECK_EQUAL( *cit, i );
    
    i = 0;
    for( Prob::iterator it = x.begin(); it != x.end(); it++, i++ )
        *it = 4 - i;
    
    i = 0;
    for( Prob::const_iterator it = x.begin(); it != x.end(); it++, i++ )
        BOOST_CHECK_EQUAL( *it, 4 - i );

    i = 0;
    for( Prob::const_reverse_iterator crit = x.rbegin(); crit != x.rend(); crit++, i++ )
        BOOST_CHECK_EQUAL( *crit, i );

    i = 0;
    for( Prob::reverse_iterator rit = x.rbegin(); rit != x.rend(); rit++, i++ )
        *rit = 2 * i;
    
    i = 0;
    for( Prob::const_reverse_iterator crit = x.rbegin(); crit != x.rend(); crit++, i++ )
        BOOST_CHECK_EQUAL( *crit, 2 * i );
}
#endif


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    Prob x( 5, 0.0 );
    for( size_t i = 0; i < x.size(); i++ )
        x.set( i, 2.0 - i );

    // test accumulate, min, max, sum, sumAbs, maxAbs
    BOOST_CHECK_EQUAL( x.sum(), 0.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( 0.0, fo_id<Real>() ), 0.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( 1.0, fo_id<Real>() ), 1.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( -1.0, fo_id<Real>() ), -1.0 );
    BOOST_CHECK_EQUAL( x.max(), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -INFINITY, fo_id<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 3.0, fo_id<Real>(), false ), 3.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -5.0, fo_id<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.min(), -2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( INFINITY, fo_id<Real>(), true ), -2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -3.0, fo_id<Real>(), true ), -3.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 5.0, fo_id<Real>(), true ), -2.0 );
    BOOST_CHECK_EQUAL( x.sumAbs(), 6.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( 0.0, fo_abs<Real>() ), 6.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( 1.0, fo_abs<Real>() ), 7.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( -1.0, fo_abs<Real>() ), 7.0 );
    BOOST_CHECK_EQUAL( x.maxAbs(), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 0.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 1.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -1.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 3.0, fo_abs<Real>(), false ), 3.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -3.0, fo_abs<Real>(), false ), 3.0 );
    x.set( 1, 1.0 );
    BOOST_CHECK_EQUAL( x.maxAbs(), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 0.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 1.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -1.0, fo_abs<Real>(), false ), 2.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( 3.0, fo_abs<Real>(), false ), 3.0 );
    BOOST_CHECK_EQUAL( x.accumulateMax( -3.0, fo_abs<Real>(), false ), 3.0 );
    for( size_t i = 0; i < x.size(); i++ )
        x.set( i, i ? (1.0 / i) : 0.0 );
    BOOST_CHECK_EQUAL( x.accumulateSum( 0.0, fo_inv0<Real>() ), 10.0 );
    x /= x.sum();

    // test entropy
    BOOST_CHECK( x.entropy() < Prob(5).entropy() );
    for( size_t i = 1; i < 100; i++ )
        BOOST_CHECK_CLOSE( Prob(i).entropy(), std::log(i), tol );

    // test hasNaNs and hasNegatives
    BOOST_CHECK( !Prob( 3, 0.0 ).hasNaNs() );
    Real c = 0.0;
    BOOST_CHECK( Prob( 3, c / c ).hasNaNs() );
    BOOST_CHECK( !Prob( 3, 0.0 ).hasNegatives() );
    BOOST_CHECK( !Prob( 3, 1.0 ).hasNegatives() );
    BOOST_CHECK( Prob( 3, -1.0 ).hasNegatives() );
    x.set( 0, 0.0 ); x.set( 1, 0.0 ); x.set( 2, -1.0 ); x.set( 3, 1.0 ); x.set( 4, 100.0 );
    BOOST_CHECK( x.hasNegatives() );
    x.set( 2, -INFINITY );
    BOOST_CHECK( x.hasNegatives() );
    x.set( 2, INFINITY );
    BOOST_CHECK( !x.hasNegatives() );
    x.set( 2, -1.0 );

    // test argmax
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)4, (Real)100.0 ) );
    x.set( 4, 0.5 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)3, (Real)1.0 ) );
    x.set( 3, -2.0 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)4, (Real)0.5 ) );
    x.set( 4, -1.0 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)0, (Real)0.0 ) );
    x.set( 0, -2.0 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)1, (Real)0.0 ) );
    x.set( 1, -3.0 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)2, (Real)-1.0 ) );
    x.set( 2, -2.0 );
    BOOST_CHECK( x.argmax() == std::make_pair( (size_t)4, (Real)-1.0 ) );

    // test draw
    for( size_t i = 0; i < x.size(); i++ )
        x.set( i, i ? (1.0 / i) : 0.0 );
    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        BOOST_CHECK( x.draw() < x.size() );
        BOOST_CHECK( x.draw() != 0 );
    }
    x.set( 2, 0.0 );
    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        BOOST_CHECK( x.draw() < x.size() );
        BOOST_CHECK( x.draw() != 0 );
        BOOST_CHECK( x.draw() != 2 );
    }
    x.set( 4, 0.0 );
    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        BOOST_CHECK( x.draw() < x.size() );
        BOOST_CHECK( x.draw() != 0 );
        BOOST_CHECK( x.draw() != 2 );
        BOOST_CHECK( x.draw() != 4 );
    }
    x.set( 1, 0.0 );
    for( size_t repeat = 0; repeat < 10000; repeat++ )
        BOOST_CHECK( x.draw() == 3 );

    // test <, ==
    Prob a(3, 1.0), b(3, 1.0);
    BOOST_CHECK( !(a < b) );
    BOOST_CHECK( !(b < a) );
    BOOST_CHECK( a == b );
    a.set( 0, 0.0 );
    BOOST_CHECK( a < b );
    BOOST_CHECK( !(b < a) );
    BOOST_CHECK( !(a == b) );
    b.set( 2, 0.0 );
    BOOST_CHECK( a < b );
    BOOST_CHECK( !(b < a) );
    BOOST_CHECK( !(a == b) );
    b.set( 0, 0.0 );
    BOOST_CHECK( !(a < b) );
    BOOST_CHECK( b < a );
    BOOST_CHECK( !(a == b) );
    a.set( 1, 0.0 );
    BOOST_CHECK( a < b );
    BOOST_CHECK( !(b < a) );
    BOOST_CHECK( !(a == b) );
    b.set( 1, 0.0 );
    BOOST_CHECK( !(a < b) );
    BOOST_CHECK( b < a );
    BOOST_CHECK( !(a == b) );
    a.set( 2, 0.0 );
    BOOST_CHECK( !(a < b) );
    BOOST_CHECK( !(b < a) );
    BOOST_CHECK( a == b );
}


BOOST_AUTO_TEST_CASE( UnaryTransformationsTest ) {
    Prob x( 3 );
    x.set( 0, -2.0 );
    x.set( 1, 0.0 );
    x.set( 2, 2.0 );

    Prob y = -x;
    Prob z = x.pwUnaryTr( std::negate<Real>() );
    BOOST_CHECK_EQUAL( y[0], 2.0 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], -2.0 );
    BOOST_CHECK( y == z );

    y = x.abs();
    z = x.pwUnaryTr( fo_abs<Real>() );
    BOOST_CHECK_EQUAL( y[0], 2.0 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 2.0 );
    BOOST_CHECK( y == z );

    y = x.exp();
    z = x.pwUnaryTr( fo_exp<Real>() );
    BOOST_CHECK_CLOSE( y[0], std::exp(-2.0), tol );
    BOOST_CHECK_EQUAL( y[1], 1.0 );
    BOOST_CHECK_CLOSE( y[2], 1.0 / y[0], tol );
    BOOST_CHECK( y == z );

    y = x.log(false);
    z = x.pwUnaryTr( fo_log<Real>() );
    BOOST_CHECK( isnan( y[0] ) );
    BOOST_CHECK_EQUAL( y[1], -INFINITY );
    BOOST_CHECK_CLOSE( y[2], std::log(2.0), tol );
    BOOST_CHECK( !(y == z) );
    y.set( 0, 0.0 );
    z.set( 0, 0.0 );
    BOOST_CHECK( y == z );

    y = x.log(true);
    z = x.pwUnaryTr( fo_log0<Real>() );
    BOOST_CHECK( isnan( y[0] ) );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], std::log(2.0) );
    BOOST_CHECK( !(y == z) );
    y.set( 0, 0.0 );
    z.set( 0, 0.0 );
    BOOST_CHECK( y == z );

    y = x.inverse(false);
    z = x.pwUnaryTr( fo_inv<Real>() );
    BOOST_CHECK_EQUAL( y[0], -0.5 );
    BOOST_CHECK_EQUAL( y[1], INFINITY );
    BOOST_CHECK_EQUAL( y[2], 0.5 );
    BOOST_CHECK( y == z );

    y = x.inverse(true);
    z = x.pwUnaryTr( fo_inv0<Real>() );
    BOOST_CHECK_EQUAL( y[0], -0.5 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 0.5 );
    BOOST_CHECK( y == z );

    x.set( 0, 2.0 );
    y = x.normalized();
    BOOST_CHECK_EQUAL( y[0], 0.5 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 0.5 );

    y = x.normalized( NORMPROB );
    BOOST_CHECK_EQUAL( y[0], 0.5 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 0.5 );

    x.set( 0, -2.0 );
    y = x.normalized( NORMLINF );
    BOOST_CHECK_EQUAL( y[0], -1.0 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 1.0 );
}


BOOST_AUTO_TEST_CASE( UnaryOperationsTest ) {
    Prob xorg(3);
    xorg.set( 0, 2.0 );
    xorg.set( 1, 0.0 );
    xorg.set( 2, 1.0 );
    Prob y(3);

    Prob x = xorg;
    BOOST_CHECK( x.setUniform() == Prob(3) );
    BOOST_CHECK( x == Prob(3) );

    y.set( 0, std::exp(2.0) );
    y.set( 1, 1.0 );
    y.set( 2, std::exp(1.0) );
    x = xorg;
    BOOST_CHECK( x.takeExp() == y );
    BOOST_CHECK( x == y );
    x = xorg;
    BOOST_CHECK( x.pwUnaryOp( fo_exp<Real>() ) == y );
    BOOST_CHECK( x == y );

    y.set( 0, std::log(2.0) );
    y.set( 1, -INFINITY );
    y.set( 2, 0.0 );
    x = xorg;
    BOOST_CHECK( x.takeLog() == y );
    BOOST_CHECK( x == y );
    x = xorg;
    BOOST_CHECK( x.takeLog(false) == y );
    BOOST_CHECK( x == y );
    x = xorg;
    BOOST_CHECK( x.pwUnaryOp( fo_log<Real>() ) == y );
    BOOST_CHECK( x == y );

    y.set( 1, 0.0 );
    x = xorg;
    BOOST_CHECK( x.takeLog(true) == y );
    BOOST_CHECK( x == y );
    x = xorg;
    BOOST_CHECK( x.pwUnaryOp( fo_log0<Real>() ) == y );
    BOOST_CHECK( x == y );

    y.set( 0, 2.0 / 3.0 );
    y.set( 1, 0.0 / 3.0 );
    y.set( 2, 1.0 / 3.0 );
    x = xorg;
    BOOST_CHECK_EQUAL( x.normalize(), 3.0 );
    BOOST_CHECK( x == y );

    x = xorg;
    BOOST_CHECK_EQUAL( x.normalize( NORMPROB ), 3.0 );
    BOOST_CHECK( x == y );

    y.set( 0, 2.0 / 2.0 );
    y.set( 1, 0.0 / 2.0 );
    y.set( 2, 1.0 / 2.0 );
    x = xorg;
    BOOST_CHECK_EQUAL( x.normalize( NORMLINF ), 2.0 );
    BOOST_CHECK( x == y );

    xorg.set( 0, -2.0 );
    y.set( 0, 2.0 );
    y.set( 1, 0.0 );
    y.set( 2, 1.0 );
    x = xorg;
    BOOST_CHECK( x.takeAbs() == y );
    BOOST_CHECK( x == y );

    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        x.randomize();
        for( size_t i = 0; i < x.size(); i++ ) {
            BOOST_CHECK( x[i] < 1.0 );
            BOOST_CHECK( x[i] >= 0.0 );
        }
    }
}


BOOST_AUTO_TEST_CASE( ScalarTransformationsTest ) {
    Prob x(3);
    x.set( 0, 2.0 );
    x.set( 1, 0.0 );
    x.set( 2, 1.0 );
    Prob y(3);

    y.set( 0, 3.0 ); y.set( 1, 1.0 ); y.set( 2, 2.0 );
    BOOST_CHECK( (x + 1.0) == y );
    y.set( 0, 0.0 ); y.set( 1, -2.0 ); y.set( 2, -1.0 );
    BOOST_CHECK( (x + (-2.0)) == y );

    y.set( 0, 1.0 ); y.set( 1, -1.0 ); y.set( 2, 0.0 );
    BOOST_CHECK( (x - 1.0) == y );
    y.set( 0, 4.0 ); y.set( 1, 2.0 ); y.set( 2, 3.0 );
    BOOST_CHECK( (x - (-2.0)) == y );

    BOOST_CHECK( (x * 1.0) == x );
    y.set( 0, 4.0 ); y.set( 1, 0.0 ); y.set( 2, 2.0 );
    BOOST_CHECK( (x * 2.0) == y );
    y.set( 0, -1.0 ); y.set( 1, 0.0 ); y.set( 2, -0.5 );
    BOOST_CHECK( (x * -0.5) == y );

    BOOST_CHECK( (x / 1.0) == x );
    y.set( 0, 1.0 ); y.set( 1, 0.0 ); y.set( 2, 0.5 );
    BOOST_CHECK( (x / 2.0) == y );
    y.set( 0, -4.0 ); y.set( 1, 0.0 ); y.set( 2, -2.0 );
    BOOST_CHECK( (x / -0.5) == y );
    BOOST_CHECK( (x / 0.0) == Prob(3, 0.0) );

    BOOST_CHECK( (x ^ 1.0) == x );
    BOOST_CHECK( (x ^ 0.0) == Prob(3, 1.0) );
    y.set( 0, 4.0 ); y.set( 1, 0.0 ); y.set( 2, 1.0 );
    BOOST_CHECK( (x ^ 2.0) == y );
    y.set( 0, 1.0 / std::sqrt(2.0) ); y.set( 1, INFINITY ); y.set( 2, 1.0 );
    Prob z = (x ^ -0.5);
    BOOST_CHECK_CLOSE( z[0], y[0], tol );
    BOOST_CHECK_EQUAL( z[1], y[1] );
    BOOST_CHECK_CLOSE( z[2], y[2], tol );
}


BOOST_AUTO_TEST_CASE( ScalarOperationsTest ) {
    Prob xorg(3), x(3);
    xorg.set( 0, 2.0 );
    xorg.set( 1, 0.0 );
    xorg.set( 2, 1.0 );
    Prob y(3);

    x = xorg;
    BOOST_CHECK( x.fill( 1.0 ) == Prob(3, 1.0) );
    BOOST_CHECK( x == Prob(3, 1.0) );
    BOOST_CHECK( x.fill( 2.0 ) == Prob(3, 2.0) );
    BOOST_CHECK( x == Prob(3, 2.0) );
    BOOST_CHECK( x.fill( 0.0 ) == Prob(3, 0.0) );
    BOOST_CHECK( x == Prob(3, 0.0) );

    x = xorg;
    y.set( 0, 3.0 ); y.set( 1, 1.0 ); y.set( 2, 2.0 );
    BOOST_CHECK( (x += 1.0) == y );
    BOOST_CHECK( x == y );
    y.set( 0, 1.0 ); y.set( 1, -1.0 ); y.set( 2, 0.0 );
    BOOST_CHECK( (x += -2.0) == y );
    BOOST_CHECK( x == y );

    x = xorg;
    y.set( 0, 1.0 ); y.set( 1, -1.0 ); y.set( 2, 0.0 );
    BOOST_CHECK( (x -= 1.0) == y );
    BOOST_CHECK( x == y );
    y.set( 0, 3.0 ); y.set( 1, 1.0 ); y.set( 2, 2.0 );
    BOOST_CHECK( (x -= -2.0) == y );
    BOOST_CHECK( x == y );

    x = xorg;
    BOOST_CHECK( (x *= 1.0) == x );
    BOOST_CHECK( x == x );
    y.set( 0, 4.0 ); y.set( 1, 0.0 ); y.set( 2, 2.0 );
    BOOST_CHECK( (x *= 2.0) == y );
    BOOST_CHECK( x == y );
    y.set( 0, -1.0 ); y.set( 1, 0.0 ); y.set( 2, -0.5 );
    BOOST_CHECK( (x *= -0.25) == y );
    BOOST_CHECK( x == y );

    x = xorg;
    BOOST_CHECK( (x /= 1.0) == x );
    BOOST_CHECK( x == x );
    y.set( 0, 1.0 ); y.set( 1, 0.0 ); y.set( 2, 0.5 );
    BOOST_CHECK( (x /= 2.0) == y );
    BOOST_CHECK( x == y );
    y.set( 0, -4.0 ); y.set( 1, 0.0 ); y.set( 2, -2.0 );
    BOOST_CHECK( (x /= -0.25) == y );
    BOOST_CHECK( x == y );
    BOOST_CHECK( (x /= 0.0) == Prob(3, 0.0) );
    BOOST_CHECK( x == Prob(3, 0.0) );

    x = xorg;
    BOOST_CHECK( (x ^= 1.0) == x );
    BOOST_CHECK( x == x );
    BOOST_CHECK( (x ^= 0.0) == Prob(3, 1.0) );
    BOOST_CHECK( x == Prob(3, 1.0) );
    x = xorg;
    y.set( 0, 4.0 ); y.set( 1, 0.0 ); y.set( 2, 1.0 );
    BOOST_CHECK( (x ^= 2.0) == y );
    BOOST_CHECK( x == y );
    y.set( 0, 0.5 ); y.set( 1, INFINITY ); y.set( 2, 1.0 );
    BOOST_CHECK( (x ^= -0.5) == y );
    BOOST_CHECK( x == y );
}


BOOST_AUTO_TEST_CASE( VectorOperationsTest ) {
    size_t N = 6;
    Prob xorg(N), x(N);
    xorg.set( 0, 2.0 ); xorg.set( 1, 0.0 ); xorg.set( 2, 1.0 ); xorg.set( 3, 0.0 ); xorg.set( 4, 2.0 ); xorg.set( 5, 3.0 );
    Prob y(N);
    y.set( 0, 0.5 ); y.set( 1, -1.0 ); y.set( 2, 0.0 ); y.set( 3, 0.0 ); y.set( 4, -2.0 ); y.set( 5, 3.0 );
    Prob z(N), r(N);

    z.set( 0, 2.5 ); z.set( 1, -1.0 ); z.set( 2, 1.0 ); z.set( 3, 0.0 ); z.set( 4, 0.0 ); z.set( 5, 6.0 );
    x = xorg;
    r = (x += y);
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    BOOST_CHECK( x == z );
    x = xorg;
    BOOST_CHECK( x.pwBinaryOp( y, std::plus<Real>() ) == z );
    BOOST_CHECK( x == z );

    z.set( 0, 1.5 ); z.set( 1, 1.0 ); z.set( 2, 1.0 ); z.set( 3, 0.0 ); z.set( 4, 4.0 ); z.set( 5, 0.0 );
    x = xorg;
    r = (x -= y);
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    BOOST_CHECK( x == z );
    x = xorg;
    BOOST_CHECK( x.pwBinaryOp( y, std::minus<Real>() ) == z );
    BOOST_CHECK( x == z );

    z.set( 0, 1.0 ); z.set( 1, 0.0 ); z.set( 2, 0.0 ); z.set( 3, 0.0 ); z.set( 4, -4.0 ); z.set( 5, 9.0 );
    x = xorg;
    r = (x *= y);
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    BOOST_CHECK( x == z );
    x = xorg;
    BOOST_CHECK( x.pwBinaryOp( y, std::multiplies<Real>() ) == z );
    BOOST_CHECK( x == z );

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, 0.0 ); z.set( 3, 0.0 ); z.set( 4, -1.0 ); z.set( 5, 1.0 );
    x = xorg;
    r = (x /= y);
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    BOOST_CHECK( x == z );
    x = xorg;
    BOOST_CHECK( x.pwBinaryOp( y, fo_divides0<Real>() ) == z );
    BOOST_CHECK( x == z );

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, INFINITY ); /*z.set( 3, INFINITY );*/ z.set( 4, -1.0 ); z.set( 5, 1.0 );
    x = xorg;
    r = (x.divide( y ));
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_CLOSE( r[1], z[1], tol );
    BOOST_CHECK_EQUAL( r[2], z[2] );
    BOOST_CHECK( isnan(r[3]) );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );
    x.set( 3, 0.0 ); r.set( 3, 0.0 );
    BOOST_CHECK( x == r );
    x = xorg;
    r = x.pwBinaryOp( y, std::divides<Real>() );
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_CLOSE( r[1], z[1], tol );
    BOOST_CHECK_EQUAL( r[2], z[2] );
    BOOST_CHECK( isnan(r[3]) );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );
    x.set( 3, 0.0 ); r.set( 3, 0.0 );
    BOOST_CHECK( x == r );

    z.set( 0, std::sqrt(2.0) ); z.set( 1, INFINITY ); z.set( 2, 1.0 ); z.set( 3, 1.0 ); z.set( 4, 0.25 ); z.set( 5, 27.0 );
    x = xorg;
    r = (x ^= y);
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_EQUAL( r[1], z[1] );
    BOOST_CHECK_CLOSE( r[2], z[2], tol );
    BOOST_CHECK_CLOSE( r[3], z[3], tol );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );
    BOOST_CHECK( x == z );
    x = xorg;
    BOOST_CHECK( x.pwBinaryOp( y, fo_pow<Real>() ) == z );
    BOOST_CHECK( x == z );
}


BOOST_AUTO_TEST_CASE( VectorTransformationsTest ) {
    size_t N = 6;
    Prob x(N);
    x.set( 0, 2.0 ); x.set( 1, 0.0 ); x.set( 2, 1.0 ); x.set( 3, 0.0 ); x.set( 4, 2.0 ); x.set( 5, 3.0 );
    Prob y(N);
    y.set( 0, 0.5 ); y.set( 1, -1.0 ); y.set( 2, 0.0 ); y.set( 3, 0.0 ); y.set( 4, -2.0 ); y.set( 5, 3.0 );
    Prob z(N), r(N);

    z.set( 0, 2.5 ); z.set( 1, -1.0 ); z.set( 2, 1.0 ); z.set( 3, 0.0 ); z.set( 4, 0.0 ); z.set( 5, 6.0 );
    r = x + y;
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    z = x.pwBinaryTr( y, std::plus<Real>() );
    BOOST_CHECK( r == z );

    z.set( 0, 1.5 ); z.set( 1, 1.0 ); z.set( 2, 1.0 ); z.set( 3, 0.0 ); z.set( 4, 4.0 ); z.set( 5, 0.0 );
    r = x - y;
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    z = x.pwBinaryTr( y, std::minus<Real>() );
    BOOST_CHECK( r == z );

    z.set( 0, 1.0 ); z.set( 1, 0.0 ); z.set( 2, 0.0 ); z.set( 3, 0.0 ); z.set( 4, -4.0 ); z.set( 5, 9.0 );
    r = x * y;
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    z = x.pwBinaryTr( y, std::multiplies<Real>() );
    BOOST_CHECK( r == z );

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, 0.0 ); z.set( 3, 0.0 ); z.set( 4, -1.0 ); z.set( 5, 1.0 );
    r = x / y;
    for( size_t i = 0; i < N; i++ )
        BOOST_CHECK_CLOSE( r[i], z[i], tol );
    z = x.pwBinaryTr( y, fo_divides0<Real>() );
    BOOST_CHECK( r == z );

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, INFINITY ); /*z.set( 3, INFINITY );*/ z.set( 4, -1.0 ); z.set( 5, 1.0 );
    r = x.divided_by( y );
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_CLOSE( r[1], z[1], tol );
    BOOST_CHECK_EQUAL( r[2], z[2] );
    BOOST_CHECK( isnan(r[3]) );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );
    z = x.pwBinaryTr( y, std::divides<Real>() );
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_CLOSE( r[1], z[1], tol );
    BOOST_CHECK_EQUAL( r[2], z[2] );
    BOOST_CHECK( isnan(r[3]) );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );

    z.set( 0, std::sqrt(2.0) ); z.set( 1, INFINITY ); z.set( 2, 1.0 ); z.set( 3, 1.0 ); z.set( 4, 0.25 ); z.set( 5, 27.0 );
    r = x ^ y;
    BOOST_CHECK_CLOSE( r[0], z[0], tol );
    BOOST_CHECK_EQUAL( r[1], z[1] );
    BOOST_CHECK_CLOSE( r[2], z[2], tol );
    BOOST_CHECK_CLOSE( r[3], z[3], tol );
    BOOST_CHECK_CLOSE( r[4], z[4], tol );
    BOOST_CHECK_CLOSE( r[5], z[5], tol );
    z = x.pwBinaryTr( y, fo_pow<Real>() );
    BOOST_CHECK( r == z );
}


BOOST_AUTO_TEST_CASE( RelatedFunctionsTest ) {
    Prob x(3), y(3), z(3);
    x.set( 0, 0.2 );
    x.set( 1, 0.8 );
    x.set( 2, 0.0 );
    y.set( 0, 0.0 );
    y.set( 1, 0.6 );
    y.set( 2, 0.4 );

    z = min( x, y );
    BOOST_CHECK_EQUAL( z[0], 0.0 );
    BOOST_CHECK_EQUAL( z[1], 0.6 );
    BOOST_CHECK_EQUAL( z[2], 0.0 );
    z = max( x, y );
    BOOST_CHECK_EQUAL( z[0], 0.2 );
    BOOST_CHECK_EQUAL( z[1], 0.8 );
    BOOST_CHECK_EQUAL( z[2], 0.4 );

    BOOST_CHECK_EQUAL( dist( x, x, DISTL1 ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, DISTL1 ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTL1 ), 0.2 + 0.2 + 0.4 );
    BOOST_CHECK_EQUAL( dist( y, x, DISTL1 ), 0.2 + 0.2 + 0.4 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTL1 ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTL1 ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, DISTLINF ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, DISTLINF ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTLINF ), 0.4 );
    BOOST_CHECK_EQUAL( dist( y, x, DISTLINF ), 0.4 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTLINF ), x.innerProduct( y, 0.0, fo_max<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTLINF ), y.innerProduct( x, 0.0, fo_max<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, DISTTV ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, DISTTV ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTTV ), 0.5 * (0.2 + 0.2 + 0.4) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTTV ), 0.5 * (0.2 + 0.2 + 0.4) );
    BOOST_CHECK_EQUAL( dist( x, y, DISTTV ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( y, x, DISTTV ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( x, x, DISTKL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, DISTKL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTKL ), INFINITY );
    BOOST_CHECK_EQUAL( dist( y, x, DISTKL ), INFINITY );
    BOOST_CHECK_EQUAL( dist( x, y, DISTKL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTKL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, DISTHEL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, DISTHEL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, DISTHEL ), 0.5 * (0.2 + std::pow(std::sqrt(0.8) - std::sqrt(0.6), 2.0) + 0.4) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTHEL ), 0.5 * (0.2 + std::pow(std::sqrt(0.8) - std::sqrt(0.6), 2.0) + 0.4) );
    BOOST_CHECK_EQUAL( dist( x, y, DISTHEL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_Hellinger<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( y, x, DISTHEL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_Hellinger<Real>() ) / 2.0 );
    x.set( 1, 0.7 ); x.set( 2, 0.1 );
    y.set( 0, 0.1 ); y.set( 1, 0.5 );
    BOOST_CHECK_CLOSE( dist( x, y, DISTKL ), 0.2 * std::log(0.2 / 0.1) + 0.7 * std::log(0.7 / 0.5) + 0.1 * std::log(0.1 / 0.4), tol );
    BOOST_CHECK_CLOSE( dist( y, x, DISTKL ), 0.1 * std::log(0.1 / 0.2) + 0.5 * std::log(0.5 / 0.7) + 0.4 * std::log(0.4 / 0.1), tol );
    BOOST_CHECK_EQUAL( dist( x, y, DISTKL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, DISTKL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_KL<Real>() ) );

    Prob xx(4), yy(4);
    for( size_t i = 0; i < 3; i++ ) {
        xx.set( i, x[i] );
        yy.set( i, y[i] );
    }
    std::stringstream ss;
    ss << xx;
    std::string s;
    std::getline( ss, s );
#ifdef DAI_SPARSE
    BOOST_CHECK_EQUAL( s, std::string("(size:4, def:0.25, 0:0.2, 1:0.7, 2:0.1)") );
#else
    BOOST_CHECK_EQUAL( s, std::string("(0.2, 0.7, 0.1, 0.25)") );
#endif
    std::stringstream ss2;
    ss2 << yy;
    std::getline( ss2, s );
#ifdef DAI_SPARSE
    BOOST_CHECK_EQUAL( s, std::string("(size:4, def:0.25, 0:0.1, 1:0.5, 2:0.4)") );
#else
    BOOST_CHECK_EQUAL( s, std::string("(0.1, 0.5, 0.4, 0.25)") );
#endif

    z = min( x, y );
    BOOST_CHECK_EQUAL( z[0], 0.1 );
    BOOST_CHECK_EQUAL( z[1], 0.5 );
    BOOST_CHECK_EQUAL( z[2], 0.1 );
    z = max( x, y );
    BOOST_CHECK_EQUAL( z[0], 0.2 );
    BOOST_CHECK_EQUAL( z[1], 0.7 );
    BOOST_CHECK_EQUAL( z[2], 0.4 );

    BOOST_CHECK_CLOSE( x.innerProduct( y, 0.0, std::plus<Real>(), std::multiplies<Real>() ), 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
    BOOST_CHECK_CLOSE( y.innerProduct( x, 0.0, std::plus<Real>(), std::multiplies<Real>() ), 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
    BOOST_CHECK_CLOSE( x.innerProduct( y, 1.0, std::plus<Real>(), std::multiplies<Real>() ), 1.0 + 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
    BOOST_CHECK_CLOSE( y.innerProduct( x, 1.0, std::plus<Real>(), std::multiplies<Real>() ), 1.0 + 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
    BOOST_CHECK_CLOSE( x.innerProduct( y, -1.0, std::plus<Real>(), std::multiplies<Real>() ), -1.0 + 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
    BOOST_CHECK_CLOSE( y.innerProduct( x, -1.0, std::plus<Real>(), std::multiplies<Real>() ), -1.0 + 0.2*0.1 + 0.7*0.5 + 0.1*0.4, tol );
}
