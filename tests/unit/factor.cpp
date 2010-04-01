/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/factor.h>
#include <strstream>


using namespace dai;


const double tol = 1e-8;


#define BOOST_TEST_MODULE FactorTest


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    Factor x1;
    BOOST_CHECK_EQUAL( x1.states(), 1 );
    BOOST_CHECK( x1.p() == Prob( 1, 1.0 ) );
    BOOST_CHECK( x1.vars() == VarSet() );

    Factor x2( 5.0 );
    BOOST_CHECK_EQUAL( x2.states(), 1 );
    BOOST_CHECK( x2.p() == Prob( 1, 5.0 ) );
    BOOST_CHECK( x2.vars() == VarSet() );

    Var v1( 0, 3 );
    Factor x3( v1 );
    BOOST_CHECK_EQUAL( x3.states(), 3 );
    BOOST_CHECK( x3.p() == Prob( 3, 1.0 / 3.0 ) );
    BOOST_CHECK( x3.vars() == VarSet( v1 ) );
    BOOST_CHECK_EQUAL( x3[0], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x3[1], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x3[2], 1.0 / 3.0 );

    Var v2( 1, 2 );
    Factor x4( VarSet( v1, v2 ) );
    BOOST_CHECK_EQUAL( x4.states(), 6 );
    BOOST_CHECK( x4.p() == Prob( 6, 1.0 / 6.0 ) );
    BOOST_CHECK( x4.vars() == VarSet( v1, v2 ) );
    for( size_t i = 0; i < 6; i++ )
        BOOST_CHECK_EQUAL( x4[i], 1.0 / 6.0 );

    Factor x5( VarSet( v1, v2 ), 1.0 );
    BOOST_CHECK_EQUAL( x5.states(), 6 );
    BOOST_CHECK( x5.p() == Prob( 6, 1.0 ) );
    BOOST_CHECK( x5.vars() == VarSet( v1, v2 ) );
    for( size_t i = 0; i < 6; i++ )
        BOOST_CHECK_EQUAL( x5[i], 1.0 );

    std::vector<Real> x( 6, 1.0 );
    for( size_t i = 0; i < 6; i++ )
        x[i] = 10.0 - i;
    Factor x6( VarSet( v1, v2 ), x );
    BOOST_CHECK_EQUAL( x6.states(), 6 );
    BOOST_CHECK( x6.vars() == VarSet( v1, v2 ) );
    for( size_t i = 0; i < 6; i++ )
        BOOST_CHECK_EQUAL( x6[i], x[i] );

    x.resize( 4 );
    BOOST_CHECK_THROW( Factor x7( VarSet( v1, v2 ), x ), Exception );

    x.resize( 6 );
    x[4] = 10.0 - 4; x[5] = 10.0 - 5;
    Factor x8( VarSet( v2, v1 ), &(x[0]) );
    BOOST_CHECK_EQUAL( x8.states(), 6 );
    BOOST_CHECK( x8.vars() == VarSet( v1, v2 ) );
    for( size_t i = 0; i < 6; i++ )
        BOOST_CHECK_EQUAL( x8[i], x[i] );

    Prob xx( x );
    Factor x9( VarSet( v2, v1 ), xx );
    BOOST_CHECK_EQUAL( x9.states(), 6 );
    BOOST_CHECK( x9.vars() == VarSet( v1, v2 ) );
    for( size_t i = 0; i < 6; i++ )
        BOOST_CHECK_EQUAL( x9[i], x[i] );

    xx.resize( 4 );
    BOOST_CHECK_THROW( Factor x10( VarSet( v2, v1 ), xx ), Exception );

    std::vector<Real> w;
    w.push_back( 0.1 );
    w.push_back( 3.5 );
    w.push_back( 2.8 );
    w.push_back( 6.3 );
    w.push_back( 8.4 );
    w.push_back( 0.0 );
    w.push_back( 7.4 );
    w.push_back( 2.4 );
    w.push_back( 8.9 );
    w.push_back( 1.3 );
    w.push_back( 1.6 );
    w.push_back( 2.6 );
    Var v4( 4, 3 );
    Var v8( 8, 2 );
    Var v7( 7, 2 );
    std::vector<Var> vars;
    vars.push_back( v4 );
    vars.push_back( v8 );
    vars.push_back( v7 );
    Factor x11( vars, w );
    BOOST_CHECK_EQUAL( x11.states(), 12 );
    BOOST_CHECK( x11.vars() == VarSet( vars.begin(), vars.end() ) );
    BOOST_CHECK_EQUAL( x11[0], 0.1 );
    BOOST_CHECK_EQUAL( x11[1], 3.5 );
    BOOST_CHECK_EQUAL( x11[2], 2.8 );
    BOOST_CHECK_EQUAL( x11[3], 7.4 );
    BOOST_CHECK_EQUAL( x11[4], 2.4 );
    BOOST_CHECK_EQUAL( x11[5], 8.9 );
    BOOST_CHECK_EQUAL( x11[6], 6.3 );
    BOOST_CHECK_EQUAL( x11[7], 8.4 );
    BOOST_CHECK_EQUAL( x11[8], 0.0 );
    BOOST_CHECK_EQUAL( x11[9], 1.3 );
    BOOST_CHECK_EQUAL( x11[10], 1.6 );
    BOOST_CHECK_EQUAL( x11[11], 2.6 );

    Factor x12( x11 );
    BOOST_CHECK( x12 == x11 );
    
    Factor x13 = x12;
    BOOST_CHECK( x13 == x11 );
}


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    Factor x( Var( 5, 5 ), 0.0 );
    for( size_t i = 0; i < x.states(); i++ )
        x.set( i, 2.0 - i );

    // test min, max, sum, sumAbs, maxAbs
    BOOST_CHECK_EQUAL( x.sum(), 0.0 );
    BOOST_CHECK_EQUAL( x.max(), 2.0 );
    BOOST_CHECK_EQUAL( x.min(), -2.0 );
    BOOST_CHECK_EQUAL( x.sumAbs(), 6.0 );
    BOOST_CHECK_EQUAL( x.maxAbs(), 2.0 );
    x.set( 1, 1.0 );
    BOOST_CHECK_EQUAL( x.maxAbs(), 2.0 );
    x /= x.sum();

    // test entropy
    BOOST_CHECK( x.entropy() < Prob(5).entropy() );
    for( size_t i = 1; i < 100; i++ )
        BOOST_CHECK_CLOSE( Factor( Var(0,i) ).entropy(), std::log(i), tol );

    // test hasNaNs and hasNegatives
    BOOST_CHECK( !Factor( 0.0 ).hasNaNs() );
    Real c = 0.0;
    BOOST_CHECK( Factor( c / c ).hasNaNs() );
    BOOST_CHECK( !Factor( 0.0 ).hasNegatives() );
    BOOST_CHECK( !Factor( 1.0 ).hasNegatives() );
    BOOST_CHECK( Factor( -1.0 ).hasNegatives() );
    x.set( 0, 0.0 ); x.set( 1, 0.0 ); x.set( 2, -1.0 ); x.set( 3, 1.0 ); x.set( 4, 100.0 );
    BOOST_CHECK( x.hasNegatives() );
    x.set( 2, -INFINITY );
    BOOST_CHECK( x.hasNegatives() );
    x.set( 2, INFINITY );
    BOOST_CHECK( !x.hasNegatives() );
    x.set( 2, -1.0 );

    // test strength
    Var x0(0,2);
    Var x1(1,2);
    BOOST_CHECK_CLOSE( createFactorIsing( x0, x1, 1.0 ).strength( x0, x1 ), std::tanh( 1.0 ), tol );
    BOOST_CHECK_CLOSE( createFactorIsing( x0, x1, -1.0 ).strength( x0, x1 ), std::tanh( 1.0 ), tol );
    BOOST_CHECK_CLOSE( createFactorIsing( x0, x1, 0.5 ).strength( x0, x1 ), std::tanh( 0.5 ), tol );

    // test ==
    Factor a(Var(0,3)), b(Var(0,3));
    Factor d(Var(1,3));
    BOOST_CHECK( !(a == d) );
    BOOST_CHECK( !(b == d) );
    BOOST_CHECK( a == b );
    a.set( 0, 0.0 );
    BOOST_CHECK( !(a == b) );
    b.set( 2, 0.0 );
    BOOST_CHECK( !(a == b) );
    b.set( 0, 0.0 );
    BOOST_CHECK( !(a == b) );
    a.set( 1, 0.0 );
    BOOST_CHECK( !(a == b) );
    b.set( 1, 0.0 );
    BOOST_CHECK( !(a == b) );
    a.set( 2, 0.0 );
    BOOST_CHECK( a == b );
}

/*
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

    y = x.normalized( Prob::NORMPROB );
    BOOST_CHECK_EQUAL( y[0], 0.5 );
    BOOST_CHECK_EQUAL( y[1], 0.0 );
    BOOST_CHECK_EQUAL( y[2], 0.5 );

    x.set( 0, -2.0 );
    y = x.normalized( Prob::NORMLINF );
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
    BOOST_CHECK_EQUAL( x.normalize( Prob::NORMPROB ), 3.0 );
    BOOST_CHECK( x == y );

    y.set( 0, 2.0 / 2.0 );
    y.set( 1, 0.0 / 2.0 );
    y.set( 2, 1.0 / 2.0 );
    x = xorg;
    BOOST_CHECK_EQUAL( x.normalize( Prob::NORMLINF ), 2.0 );
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

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, INFINITY ); 
    // z.set( 3, INFINITY );
    z.set( 4, -1.0 ); z.set( 5, 1.0 );
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

    z.set( 0, 4.0 ); z.set( 1, 0.0 ); z.set( 2, INFINITY ); 
    // z.set( 3, INFINITY );
    z.set( 4, -1.0 ); z.set( 5, 1.0 );
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

    BOOST_CHECK_EQUAL( dist( x, x, Prob::DISTL1 ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, Prob::DISTL1 ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTL1 ), 0.2 + 0.2 + 0.4 );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTL1 ), 0.2 + 0.2 + 0.4 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTL1 ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTL1 ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, Prob::DISTLINF ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, Prob::DISTLINF ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTLINF ), 0.4 );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTLINF ), 0.4 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTLINF ), x.innerProduct( y, 0.0, fo_max<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTLINF ), y.innerProduct( x, 0.0, fo_max<Real>(), fo_absdiff<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, Prob::DISTTV ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, Prob::DISTTV ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTTV ), 0.5 * (0.2 + 0.2 + 0.4) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTTV ), 0.5 * (0.2 + 0.2 + 0.4) );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTTV ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTTV ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_absdiff<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( x, x, Prob::DISTKL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, Prob::DISTKL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTKL ), INFINITY );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTKL ), INFINITY );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTKL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTKL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( x, x, Prob::DISTHEL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( y, y, Prob::DISTHEL ), 0.0 );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTHEL ), 0.5 * (0.2 + std::pow(std::sqrt(0.8) - std::sqrt(0.6), 2.0) + 0.4) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTHEL ), 0.5 * (0.2 + std::pow(std::sqrt(0.8) - std::sqrt(0.6), 2.0) + 0.4) );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTHEL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_Hellinger<Real>() ) / 2.0 );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTHEL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_Hellinger<Real>() ) / 2.0 );
    x.set( 1, 0.7 ); x.set( 2, 0.1 );
    y.set( 0, 0.1 ); y.set( 1, 0.5 );
    BOOST_CHECK_CLOSE( dist( x, y, Prob::DISTKL ), 0.2 * std::log(0.2 / 0.1) + 0.7 * std::log(0.7 / 0.5) + 0.1 * std::log(0.1 / 0.4), tol );
    BOOST_CHECK_CLOSE( dist( y, x, Prob::DISTKL ), 0.1 * std::log(0.1 / 0.2) + 0.5 * std::log(0.5 / 0.7) + 0.4 * std::log(0.4 / 0.1), tol );
    BOOST_CHECK_EQUAL( dist( x, y, Prob::DISTKL ), x.innerProduct( y, 0.0, std::plus<Real>(), fo_KL<Real>() ) );
    BOOST_CHECK_EQUAL( dist( y, x, Prob::DISTKL ), y.innerProduct( x, 0.0, std::plus<Real>(), fo_KL<Real>() ) );

    std::stringstream ss;
    ss << x;
    std::string s;
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, std::string("(0.2, 0.7, 0.1)") );
    std::stringstream ss2;
    ss2 << y;
    std::getline( ss2, s );
    BOOST_CHECK_EQUAL( s, std::string("(0.1, 0.5, 0.4)") );

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
*/
