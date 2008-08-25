/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
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


#ifndef __defined_libdai_x2x_h
#define __defined_libdai_x2x_h


#include <cmath>
#include <cstring>


namespace x2x {

    // Probability tables store -1 first, then +1
 
    /// Convert moments to cumulants upto order k
    void m2c (int N, double *x, int k);

    /// Convert cumulants to moments upto order k
    void c2m (int N, double *x, int k);

    /// Convert (generalized) weights to log probability or energy
    void w2logp (int N, double *x);

    /// Convert log probability or energy to (generalized) weights
    void logp2w (int N, double *x);

    /// Convert probability to moments
    void p2m (int N, double *x);

    /// Convert moments to probability
    void m2p (int N, double *x);

    /// Convert log probability to probability
    void logp2p (int N, double *x);

    /// Convert probability to log probability
    void p2logp (int N, double *x);

    /// Normalize a log probability table
    void logpnorm (int N, double *x);

    /// Normalize a probability table, use logpnorm whenever possible
    void pnorm (int N, double *x);

    /// Fills table with v for all entries with more than k indices
    /// used for example when cumulants or moments are converted upto some order
    void fill (int N, double *x, int k, double v);

} // end of namespace x2x


#endif
