/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/exceptions.h>
#include <strstream>


using namespace dai;


#define BOOST_TEST_MODULE ExceptionsTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( ExceptionsTest ) {
    BOOST_CHECK_THROW( DAI_THROW(NOT_IMPLEMENTED), Exception );
    BOOST_CHECK_THROW( DAI_THROW(NOT_IMPLEMENTED), std::runtime_error );
    BOOST_CHECK_THROW( DAI_THROWE(NOT_IMPLEMENTED,"Detailed error message"), Exception );
    BOOST_CHECK_THROW( DAI_THROWE(NOT_IMPLEMENTED,"Detailed error messgae"), std::runtime_error );
    BOOST_CHECK_THROW( DAI_ASSERT( 0 ), Exception );
    BOOST_CHECK_THROW( DAI_ASSERT( 0 == 1 ), std::runtime_error );

    try {
        DAI_THROW(NOT_IMPLEMENTED);
    } catch( Exception& e ) {
        BOOST_CHECK_EQUAL( e.code(), Exception::NOT_IMPLEMENTED );
        BOOST_CHECK_EQUAL( e.message(e.code()), std::string("Feature not implemented") );
    }

    try {
        DAI_THROWE(NOT_IMPLEMENTED,"Detailed error message");
    } catch( Exception& e ) {
        BOOST_CHECK_EQUAL( e.code(), Exception::NOT_IMPLEMENTED );
        BOOST_CHECK_EQUAL( e.message(e.code()), std::string("Feature not implemented") );
    }

    try {
        DAI_THROW(NOT_IMPLEMENTED);
    } catch( std::runtime_error& e ) {
        BOOST_CHECK_EQUAL( e.what(), std::string("Feature not implemented [tests/unit/exceptions_test.cpp, line 45]") );
    }

    try {
        DAI_THROWE(NOT_IMPLEMENTED,"Detailed error message");
    } catch( std::runtime_error& e ) {
        BOOST_CHECK_EQUAL( e.what(), std::string("Feature not implemented [tests/unit/exceptions_test.cpp, line 51]") );
    }
}
