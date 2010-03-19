/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/var.h>
#include <strstream>


using namespace dai;


#define BOOST_TEST_MODULE VarTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    // check constructors
    Var x;
    BOOST_CHECK_EQUAL( x.label(), 0 );
    BOOST_CHECK_EQUAL( x.states(), 0 );

    x = Var( 0, 2 );
    BOOST_CHECK_EQUAL( x.label(), 0 );
    BOOST_CHECK_EQUAL( x.states(), 2 );
}

BOOST_AUTO_TEST_CASE( AccMutTest ) {
    // check states and labels mutators
    Var x;
    x.states() = 3;
    BOOST_CHECK_EQUAL( x.states(), 3 );

    x.label() = 5;
    BOOST_CHECK_EQUAL( x.label(), 5 );
}

BOOST_AUTO_TEST_CASE( ComparisonTest ) {
    // check comparison operators
    Var x( 5, 3 );
    Var y( 6, 3 );
    Var z( 5, 3 );
    BOOST_CHECK( x < y );
    BOOST_CHECK( !(x < z) );
    BOOST_CHECK( y > x );
    BOOST_CHECK( !(z > x) );
    BOOST_CHECK( x <= y );
    BOOST_CHECK( x <= z );
    BOOST_CHECK( !(x >= y) );
    BOOST_CHECK( x >= z );
    BOOST_CHECK( !(x == y) );
    BOOST_CHECK( x == z );
    BOOST_CHECK( x != y );
    BOOST_CHECK( !(x != z) );
}

BOOST_AUTO_TEST_CASE( StreamTest ) {
    // check stream output
    Var x( 5, 3 );
    std::stringstream ss;
    ss << x;
    std::string s;
    ss >> s;
    BOOST_CHECK_EQUAL( s, "x5" );
}
