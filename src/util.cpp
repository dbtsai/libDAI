/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
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
    // Assume POSIX compliant system. We need the following for querying the system time
    #include <sys/time.h>
#endif


#ifdef CYGWIN
bool isnan( double x ) {
    return __isnand( x );  // isnan() is a macro in Cygwin (as required by C99)
}
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


Real max( const std::vector<Real> &v ) {
    if( v.size() == 0 )
        return INFINITY;
    else
        return *std::max_element( v.begin(), v.end() );
}

// Returns user+system time in seconds
double toc() {
#ifdef WINDOWS
    SYSTEMTIME tbuf;
    GetSystemTime(&tbuf);
    return( (double)(tbuf.wSecond + (double)tbuf.wMilliseconds / 1000.0) );
#else
    struct timeval tv;
    struct timezone tz;
    gettimeofday( &tv, &tz );
    return( (double)(tv.tv_sec + (double)tv.tv_usec / 1000000.0) );
#endif
}

/// Type of global random number generator
typedef boost::minstd_rand _rnd_gen_type;  // Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand

/// Global random number generator
_rnd_gen_type _rnd_gen(42U);

/// Uniform distribution with values between 0 and 1 (0 inclusive, 1 exclusive).
boost::uniform_real<Real> _uni_dist(0,1);

/// Global uniform random random number
boost::variate_generator<_rnd_gen_type&, boost::uniform_real<Real> > _uni_rnd(_rnd_gen, _uni_dist);

/// Normal distribution with mean 0 and standard deviation 1.
boost::normal_distribution<Real> _normal_dist;

/// Global random number generator with standard normal distribution
boost::variate_generator<_rnd_gen_type&, boost::normal_distribution<Real> > _normal_rnd(_rnd_gen, _normal_dist);


void rnd_seed( size_t seed ) {
    _rnd_gen.seed(seed);
}

Real rnd_uniform() {
    return _uni_rnd();
}

Real rnd_stdnormal() {
    return _normal_rnd();
}

int rnd_int( int min, int max ) {
    return (int)floor(_uni_rnd() * (max + 1 - min) + min);
}

void tokenizeString(const std::string& s, std::vector<std::string>& outTokens, const std::string& delim) {
    size_t start = 0;
    while (start < s.size()) {
        size_t end = s.find_first_of(delim, start);
        if (end > s.size())
            end = s.size();
        outTokens.push_back(s.substr(start, end - start));
        start = end + 1;
    }
}

} // end of namespace dai
