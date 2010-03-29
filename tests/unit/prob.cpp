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


#define BOOST_TEST_MODULE ProbTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    Prob x1;
    BOOST_CHECK_EQUAL( x1.size(), 0 );
    BOOST_CHECK( x1.p() == std::vector<Real>() );

    Prob x2( 3 );
    BOOST_CHECK_EQUAL( x2.size(), 3 );
    BOOST_CHECK( x2.p() == std::vector<Real>( 3, 1.0 / 3.0 ) );
    BOOST_CHECK_EQUAL( x2[0], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x2[1], 1.0 / 3.0 );
    BOOST_CHECK_EQUAL( x2[2], 1.0 / 3.0 );

    Prob x3( 4, 1.0 );
    BOOST_CHECK_EQUAL( x3.size(), 4 );
    BOOST_CHECK( x3.p() == std::vector<Real>( 4, 1.0 ) );
    BOOST_CHECK_EQUAL( x3[0], 1.0 );
    BOOST_CHECK_EQUAL( x3[1], 1.0 );
    BOOST_CHECK_EQUAL( x3[2], 1.0 );
    BOOST_CHECK_EQUAL( x3[3], 1.0 );
    x3[0] = 0.5;
    x3[1] = 1.0;
    x3[2] = 2.0;
    x3[3] = 4.0;

    Prob x4( x3.begin(), x3.end() );
    BOOST_CHECK_EQUAL( x4.size(), 4 );
    BOOST_CHECK( x4.p() == x3.p() );
    BOOST_CHECK_EQUAL( x4[0], 0.5 );
    BOOST_CHECK_EQUAL( x4[1], 1.0 );
    BOOST_CHECK_EQUAL( x4[2], 2.0 );
    BOOST_CHECK_EQUAL( x4[3], 4.0 );

    x3.p() = std::vector<Real>( 4, 2.5 );
    Prob x5( x3.begin(), x3.end(), x3.size() );
    BOOST_CHECK_EQUAL( x5.size(), 4 );
    BOOST_CHECK( x5.p() == x3.p() );
    BOOST_CHECK_EQUAL( x5[0], 2.5 );
    BOOST_CHECK_EQUAL( x5[1], 2.5 );
    BOOST_CHECK_EQUAL( x5[2], 2.5 );
    BOOST_CHECK_EQUAL( x5[3], 2.5 );

    std::vector<int> y( 3, 2 );
    Prob x6( y );
    BOOST_CHECK_EQUAL( x6.size(), 3 );
    BOOST_CHECK( x6.p() == std::vector<Real>( 3, 2.0 ) );
    BOOST_CHECK_EQUAL( x6[0], 2.0 );
    BOOST_CHECK_EQUAL( x6[1], 2.0 );
    BOOST_CHECK_EQUAL( x6[2], 2.0 );

    Prob x7( x6 );
    BOOST_CHECK( x7 == x6 );
    
    Prob x8 = x6;
    BOOST_CHECK( x8 == x6 );
}


BOOST_AUTO_TEST_CASE( IteratorTest ) {
    Prob x( 5, 0.0 );
    size_t i;
    for( i = 0; i < x.size(); i++ )
        x[i] = i;

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


BOOST_AUTO_TEST_CASE( QueriesTest ) {
}


BOOST_AUTO_TEST_CASE( UnaryTransformationsTest ) {
}


BOOST_AUTO_TEST_CASE( UnaryOperationsTest ) {
}


BOOST_AUTO_TEST_CASE( ScalarOperationsTest ) {
}


BOOST_AUTO_TEST_CASE( VectorTransformationsTest ) {
}


BOOST_AUTO_TEST_CASE( VectorOperationsTest ) {
}


BOOST_AUTO_TEST_CASE( RelatedFunctionsTest ) {
}
