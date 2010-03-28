/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2010  Joris Mooij      [joris dot mooij at libdai dot org]
 */


#define BOOST_TEST_DYN_LINK


#include <dai/index.h>
#include <strstream>
#include <map>


using namespace dai;


#define BOOST_TEST_MODULE IndexTest


#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( IndexForTest ) {
    IndexFor x;
    BOOST_CHECK( !x.valid() );
    x.reset();
    BOOST_CHECK( x.valid() );

    size_t nrVars = 5;
    std::vector<Var> vars;
    for( size_t i = 0; i < nrVars; i++ )
        vars.push_back( Var( i, i+2 ) );

    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        VarSet indexVars;
        VarSet forVars;
        for( size_t i = 0; i < 5; i++ ) {
            if( rnd(2) == 0 )
                indexVars |= vars[i];
            if( rnd(2) == 0 )
                forVars |= vars[i];
        }
        IndexFor ind( indexVars, forVars );
        size_t iter = 0;
        for( ; ind.valid(); ind++, iter++ )
            BOOST_CHECK_EQUAL( calcLinearState( indexVars, calcState( forVars, iter ) ), (size_t)ind );
        BOOST_CHECK_EQUAL( iter, forVars.nrStates() );
        iter = 0;
        ind.reset();
        for( ; ind.valid(); ++ind, iter++ )
            BOOST_CHECK_EQUAL( calcLinearState( indexVars, calcState( forVars, iter ) ), (size_t)ind );
        BOOST_CHECK_EQUAL( iter, forVars.nrStates() );
    }
}


BOOST_AUTO_TEST_CASE( PermuteTest ) {
    Permute x;

    Var x0(0, 2);
    Var x1(1, 3);
    Var x2(2, 2);
    std::vector<Var> V;
    V.push_back( x1 );
    V.push_back( x2 );
    V.push_back( x0 );
    VarSet X( V.begin(), V.end() );
    Permute sigma(V);
    BOOST_CHECK_EQUAL( sigma.sigma()[0], 2 );
    BOOST_CHECK_EQUAL( sigma.sigma()[1], 0 );
    BOOST_CHECK_EQUAL( sigma.sigma()[2], 1 );
    BOOST_CHECK_EQUAL( sigma[0], 2 );
    BOOST_CHECK_EQUAL( sigma[1], 0 );
    BOOST_CHECK_EQUAL( sigma[2], 1 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 0 ), 0 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 1 ), 2 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 2 ), 4 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 3 ), 6 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 4 ), 8 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 5 ), 10 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 6 ), 1 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 7 ), 3 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 8 ), 5 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 9 ), 7 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 10 ), 9 );
    BOOST_CHECK_EQUAL( sigma.convertLinearIndex( 11 ), 11 );

    std::vector<size_t> rs, sig;
    rs.push_back(2);
    rs.push_back(3);
    rs.push_back(2);
    sig.push_back(2);
    sig.push_back(0);
    sig.push_back(1);
    Permute tau( rs, sig );
    BOOST_CHECK( tau.sigma() == sig );
    BOOST_CHECK_EQUAL( tau[0], 2 );
    BOOST_CHECK_EQUAL( tau[1], 0 );
    BOOST_CHECK_EQUAL( tau[2], 1 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 0 ), 0 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 1 ), 2 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 2 ), 4 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 3 ), 6 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 4 ), 8 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 5 ), 10 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 6 ), 1 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 7 ), 3 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 8 ), 5 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 9 ), 7 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 10 ), 9 );
    BOOST_CHECK_EQUAL( tau.convertLinearIndex( 11 ), 11 );
}


BOOST_AUTO_TEST_CASE( multiforTest ) {
    multifor x;
    BOOST_CHECK( x.valid() );

    std::vector<size_t> ranges;
    ranges.push_back( 3 );
    ranges.push_back( 4 );
    ranges.push_back( 5 );
    multifor S(ranges);
    size_t s = 0;
    for( size_t s2 = 0; s2 < 5; s2++ )
        for( size_t s1 = 0; s1 < 4; s1++ )
            for( size_t s0 = 0; s0 < 3; s0++, s++, S++ ) {
                BOOST_CHECK( S.valid() );
                BOOST_CHECK_EQUAL( s, (size_t)S );
                BOOST_CHECK_EQUAL( S[0], s0 );
                BOOST_CHECK_EQUAL( S[1], s1 );
                BOOST_CHECK_EQUAL( S[2], s2 );
            }
    BOOST_CHECK( !S.valid() );

    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        std::vector<size_t> dims;
        size_t total = 1;
        for( size_t i = 0; i < 4; i++ ) {
            dims.push_back( rnd(3) + 1 );
            total *= dims.back();
        }
        multifor ind( dims );
        size_t iter = 0;
        for( ; ind.valid(); ind++, iter++ ) {
            BOOST_CHECK_EQUAL( (size_t)ind, iter );
            BOOST_CHECK_EQUAL( ind[0], iter % dims[0] );
            BOOST_CHECK_EQUAL( ind[1], (iter / dims[0]) % dims[1] );
            BOOST_CHECK_EQUAL( ind[2], (iter / (dims[0] * dims[1])) % dims[2] );
            BOOST_CHECK_EQUAL( ind[3], (iter / (dims[0] * dims[1] * dims[2])) % dims[3] );
        }
        BOOST_CHECK_EQUAL( iter, total );
        iter = 0;
        ind.reset();
        for( ; ind.valid(); ++ind, iter++ ) {
            BOOST_CHECK_EQUAL( (size_t)ind, iter );
            BOOST_CHECK_EQUAL( ind[0], iter % dims[0] );
            BOOST_CHECK_EQUAL( ind[1], (iter / dims[0]) % dims[1] );
            BOOST_CHECK_EQUAL( ind[2], (iter / (dims[0] * dims[1])) % dims[2] );
            BOOST_CHECK_EQUAL( ind[3], (iter / (dims[0] * dims[1] * dims[2])) % dims[3] );
        }
        BOOST_CHECK_EQUAL( iter, total );
    }
}


BOOST_AUTO_TEST_CASE( StateTest ) {
    State x;
    BOOST_CHECK( x.valid() );

    Var v0( 0, 3 );
    Var v1( 1, 4 );
    Var v2( 3, 5 );
    VarSet vars;
    vars |= v2;
    vars |= v1;
    vars |= v0;
    State S( vars );
    size_t s = 0;
    for( size_t s2 = 0; s2 < 5; s2++ )
        for( size_t s1 = 0; s1 < 4; s1++ )
            for( size_t s0 = 0; s0 < 3; s0++, s++, S++ ) {
                BOOST_CHECK( S.valid() );
                BOOST_CHECK_EQUAL( s, (size_t)S );
                BOOST_CHECK_EQUAL( S(v0), s0 );
                BOOST_CHECK_EQUAL( S(v1), s1 );
                BOOST_CHECK_EQUAL( S(v2), s2 );
                BOOST_CHECK_EQUAL( S( Var( 2, 2 ) ), 0 );
            }
    BOOST_CHECK( !S.valid() );
    S.reset();
    std::vector<std::pair<Var, size_t> > ps;
    ps.push_back( std::make_pair( Var( 2, 2 ), 1 ) );
    ps.push_back( std::make_pair( Var( 4, 2 ), 1 ) );
    S.insert( ps.begin(), ps.end() );
    BOOST_CHECK( S.valid() );
    BOOST_CHECK_EQUAL( (size_t)S, 132 );

    for( size_t repeat = 0; repeat < 10000; repeat++ ) {
        std::vector<size_t> dims;
        size_t total = 1;
        for( size_t i = 0; i < 4; i++ ) {
            dims.push_back( rnd(3) + 1 );
            total *= dims.back();
        }
        std::vector<Var> vs;
        for( size_t i = 0; i < 4; i++ )
            vs.push_back( Var( i, dims[i] ) );
        State ind( VarSet( vs.begin(), vs.end() ) );
        size_t iter = 0;
        for( ; ind.valid(); ind++, iter++ ) {
            BOOST_CHECK_EQUAL( (size_t)ind, iter );
            BOOST_CHECK_EQUAL( ind(vs[0]), iter % dims[0] );
            BOOST_CHECK_EQUAL( ind(vs[1]), (iter / dims[0]) % dims[1] );
            BOOST_CHECK_EQUAL( ind(vs[2]), (iter / (dims[0] * dims[1])) % dims[2] );
            BOOST_CHECK_EQUAL( ind(vs[3]), (iter / (dims[0] * dims[1] * dims[2])) % dims[3] );
            BOOST_CHECK_EQUAL( ind(VarSet(vs[0], vs[1])), iter % (dims[0] * dims[1]) );
            BOOST_CHECK_EQUAL( ind(VarSet(vs[1], vs[2])), (iter / dims[0]) % (dims[1] * dims[2]) );
            BOOST_CHECK_EQUAL( ind(VarSet(vs[2], vs[3])), (iter / (dims[0] * dims[1])) % (dims[2] * dims[3]) );
            BOOST_CHECK_EQUAL( ind(VarSet(vs.begin(), vs.end())), iter );
            State indcopy( VarSet(vs.begin(), vs.end()), (size_t)ind );
            BOOST_CHECK_EQUAL( ind(vs[0]), indcopy(vs[0]) );
            BOOST_CHECK_EQUAL( ind(vs[1]), indcopy(vs[1]) );
            BOOST_CHECK_EQUAL( ind(vs[2]), indcopy(vs[2]) );
            BOOST_CHECK_EQUAL( ind(vs[3]), indcopy(vs[3]) );
            State indcopy2( indcopy.get() );
            BOOST_CHECK_EQUAL( ind(vs[0]), indcopy2(vs[0]) );
            BOOST_CHECK_EQUAL( ind(vs[1]), indcopy2(vs[1]) );
            BOOST_CHECK_EQUAL( ind(vs[2]), indcopy2(vs[2]) );
            BOOST_CHECK_EQUAL( ind(vs[3]), indcopy2(vs[3]) );
            std::map<Var,size_t> indmap( ind );
            State indcopy3( indmap );
            BOOST_CHECK_EQUAL( ind(vs[0]), indcopy3(vs[0]) );
            BOOST_CHECK_EQUAL( ind(vs[1]), indcopy3(vs[1]) );
            BOOST_CHECK_EQUAL( ind(vs[2]), indcopy3(vs[2]) );
            BOOST_CHECK_EQUAL( ind(vs[3]), indcopy3(vs[3]) );
        }
        BOOST_CHECK_EQUAL( iter, total );
        iter = 0;
        ind.reset();
        for( ; ind.valid(); ++ind, iter++ ) {
            BOOST_CHECK_EQUAL( (size_t)ind, iter );
            BOOST_CHECK_EQUAL( ind(vs[0]), iter % dims[0] );
            BOOST_CHECK_EQUAL( ind(vs[1]), (iter / dims[0]) % dims[1] );
            BOOST_CHECK_EQUAL( ind(vs[2]), (iter / (dims[0] * dims[1])) % dims[2] );
            BOOST_CHECK_EQUAL( ind(vs[3]), (iter / (dims[0] * dims[1] * dims[2])) % dims[3] );
            State::const_iterator ci = ind.begin();
            BOOST_CHECK_EQUAL( (ci++)->second, iter % dims[0] );
            BOOST_CHECK_EQUAL( (ci++)->second, (iter / dims[0]) % dims[1] );
            BOOST_CHECK_EQUAL( (ci++)->second, (iter / (dims[0] * dims[1])) % dims[2] );
            BOOST_CHECK_EQUAL( (ci++)->second, (iter / (dims[0] * dims[1] * dims[2])) % dims[3] );
            BOOST_CHECK( ci == ind.end() );
        }
        BOOST_CHECK_EQUAL( iter, total );
        State::const_iterator ci = ind.begin();
        BOOST_CHECK_EQUAL( (ci++)->first, vs[0] );
        BOOST_CHECK_EQUAL( (ci++)->first, vs[1] );
        BOOST_CHECK_EQUAL( (ci++)->first, vs[2] );
        BOOST_CHECK_EQUAL( (ci++)->first, vs[3] );
        BOOST_CHECK( ci == ind.end() );
    }
}
