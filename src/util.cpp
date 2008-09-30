/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

    This file is part of libDAI.

    libDAI is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    libDAI is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libDAI; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include <dai/util.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#ifdef WINDOWS
    #include <windows.h>
    #include <boost/math/special_functions/atanh.hpp>  // for atanh
    #include <boost/math/special_functions/log1p.hpp>  // for log1p
    #include <float.h>  // for _isnan
#else
    // Assume POSIX compliant system. We need the following for querying the CPU time for this process
    #include <sys/times.h>
    #include <sys/param.h>
#endif


#ifdef WINDOWS
bool isnan( double x ) {
    return _isnan( x );
}
double atanh( double x ) {
    return boost::math::atanh( x );
}
double log1p( double x ) {
    return boost::math::log1p( x );
}
#endif


namespace dai {


// Returns user+system time in seconds
double toc() {
#ifdef WINDOWS
    SYSTEMTIME  tbuf;
    GetSystemTime(&tbuf);
    return( (double)(tbuf.wSecond + (double)tbuf.wMilliseconds / 1000.0) );
#else
    tms tbuf;
    times(&tbuf);
    return( (double)(tbuf.tms_utime + tbuf.tms_stime) / HZ );
#endif
}

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand _rnd_gen_type;

_rnd_gen_type _rnd_gen(42U);

// Define a uniform random number distribution which produces
// values between 0 and 1 (0 inclusive, 1 exclusive).
boost::uniform_real<> _uni_dist(0,1);
boost::variate_generator<_rnd_gen_type&, boost::uniform_real<> > _uni_rnd(_rnd_gen, _uni_dist);

// Define a normal distribution with mean 0 and standard deviation 1.
boost::normal_distribution<> _normal_dist;
boost::variate_generator<_rnd_gen_type&, boost::normal_distribution<> > _normal_rnd(_rnd_gen, _normal_dist);


void rnd_seed( size_t seed ) {
    _rnd_gen.seed(seed);
}

double rnd_uniform() {
    return _uni_rnd();
}

double rnd_stdnormal() {
    return _normal_rnd();
}

int rnd_int( int min, int max ) {
    return (int)floor(_uni_rnd() * (max + 1 - min) + min);
}


} // end of namespace dai
