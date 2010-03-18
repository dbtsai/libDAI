#define BOOST_TEST_DYN_LINK

#include <dai/var.h>

using namespace dai;

#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( my_test ) {
    Var x( 0, 2 );
    BOOST_CHECK( x.states() == 2 );
}
