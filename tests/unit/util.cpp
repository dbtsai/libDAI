/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/util.h>
#include <strstream>
#include <string>
#include <vector>
#include <map>
#include <set>


using namespace dai;


#define BOOST_TEST_MODULE UtilTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( IsNanTest ) {
    double x = 0.0, y = 0.0;
    BOOST_CHECK( isnan( x / y ) );
    BOOST_CHECK( !isnan( x ) );
}


BOOST_AUTO_TEST_CASE( RndTest ) {
    rnd_seed( 123 );
    Real a1 = rnd_uniform();
    Real b1 = rnd_stdnormal();
    int c1 = rnd_int(0,1);
    int d1 = rnd(2);

    rnd_seed( 123 );
    Real a2 = rnd_uniform();
    Real b2 = rnd_stdnormal();
    int c2 = rnd_int(0,1);
    int d2 = rnd(2);

    BOOST_CHECK_EQUAL( a1, a2 );
    BOOST_CHECK_EQUAL( b1, b2 );
    BOOST_CHECK_EQUAL( c1, c2 );
    BOOST_CHECK_EQUAL( d1, d2 );

    for( size_t i = 0; i < 10000; i++ ) {
        Real x = rnd_uniform();
        BOOST_CHECK( x >= 0.0 );
        BOOST_CHECK( x < 1 );
    }

    for( int min = -5; min <= 5; min++ )
        for( int max = min; max <= min + 10; max++ )
            for( size_t i = 0; i < 1000; i++ ) {
                int j = rnd_int( min, max );
                BOOST_CHECK( j >= min );
                BOOST_CHECK( j <= max );
            }

    for( int max = 1; max <= 100; max++ )
        for( size_t i = 0; i < 100; i++ ) {
            int j = rnd( max );
            BOOST_CHECK( j >= 0 );
            BOOST_CHECK( j < max );
        }
}


BOOST_AUTO_TEST_CASE( StreamTest ) {
    std::vector<int> v;
    v.push_back( 1 );
    v.push_back( 3 );
    v.push_back( 2 );
    std::stringstream ss;
    ss << v;
    std::string s;
    std::getline( ss, s );
    BOOST_CHECK_EQUAL( s, "(1, 3, 2)" );

    std::set<int> x;
    x.insert( 6 );
    x.insert( 5 );
    x.insert( 4 );
    std::stringstream ss2;
    ss2 << x;
    std::getline( ss2, s );
    BOOST_CHECK_EQUAL( s, "{4, 5, 6}" );

    std::map<int,int> y;
    y[1] = 6;
    y[3] = 5;
    y[2] = 4;
    std::stringstream ss3;
    ss3 << y;
    std::getline( ss3, s );
    BOOST_CHECK_EQUAL( s, "{1->6, 2->4, 3->5}" );

    std::pair<int,double> z;
    z.first = 5;
    z.second = 1.2345;
    std::stringstream ss4;
    ss4 << z;
    std::getline( ss4, s );
    BOOST_CHECK_EQUAL( s, "(5, 1.2345)" );
}


BOOST_AUTO_TEST_CASE( concatTest ) {
    std::vector<int> a;
    a.push_back( 0 );
    a.push_back( 1 );
    std::vector<int> b;
    b.push_back( 2 );
    b.push_back( 3 );
    b.push_back( 4 );
    std::vector<int> c;
    c.push_back( 0 );
    c.push_back( 1 );
    c.push_back( 2 );
    c.push_back( 3 );
    c.push_back( 4 );
    BOOST_CHECK( concat( a, b ) == c );
}


BOOST_AUTO_TEST_CASE( tokenizeStringTest ) {
    std::string s("Hello\tworld.\nThis is it.");
    std::vector<std::string> words;
    tokenizeString( s, words );
    BOOST_CHECK_EQUAL( words.size(), 3 );
    BOOST_CHECK_EQUAL( words[0], "Hello" );
    BOOST_CHECK_EQUAL( words[1], "world." );
    BOOST_CHECK_EQUAL( words[2], "This is it." );
    words.clear();
    tokenizeString( s, words, " " );
    BOOST_CHECK_EQUAL( words.size(), 3 );
    BOOST_CHECK_EQUAL( words[0], "Hello\tworld.\nThis" );
    BOOST_CHECK_EQUAL( words[1], "is" );
    BOOST_CHECK_EQUAL( words[2], "it." );
}
