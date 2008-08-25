/*
Copyright (C) 2005  Martijn Leisink  m.leisink at science.ru.nl

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

CHANGES:
  2006-11-20  Joris Mooij
    * removed MATLAB interface
    * put into namespace "x2x"
*/


#include <cmath>
#include <cstring>



namespace x2x {

    // helper functions to compute the sum over all partitions
    double psum (double *x, long s, int n=0);
    double psumx (double *x, long a, long s, int n) {
        // recursively process all one bits remaining after psum
        // put some of them into the same subset and call psum again for the rest
        if((s>>n)>1) {
            while(!((s>>++n)&1));
            return(psumx(x,a,s,n)+psumx(x,a^(1l<<n),s,n));
        } else {
            return(x[a]*psum(x,s^a,0));
        };
    }
    double psum (double *x, long s, int n) {
        // take the first one bit and put it in the first subset, then call psumx
        if(s>>n) {
            while(!((s>>n)&1)) ++n;
            return(psumx(x,1l<<n,s,n));
        } else return(1);
    }

    // convert moments to cumulants upto order k
    void m2c (int N, double *x, int k) {
        int *c=new int[k+1];
        long s;
        int z;
        x[0]=log(x[0]); // to get the correct answer if not normalized
        for(int b=1;b<=k;++b) {
            // start with marginals, then two-point correlations and so on
            c[z=b]=-1;
            s=0; // index into x
            do {
                while(--z>=0) s^=(1l<<(c[z]=c[z+1]+1));
                // c_ijk = m_ijk - c_i*c_jk - c_j*c_ik - c_k*c_ij - c_i*c_j*c_k
                x[s]=2*x[s]-psum(x,s);
                // increment b indices
                for(z=0;z<b&&(s^=3l<<c[z],++c[z]==N-z);++z) s^=1l<<c[z];
            } while(z<b);
        };
        delete[] c;
    }

    // convert cumulants to moments upto order k
    void c2m (int N, double *x, int k) {
        int *c=new int[k+1];
        long s;
        int z;
        for(int b=k;b>=1;--b) {
            // start with k-point cumulants, then k-1-point cumulants and so on
            c[z=b]=-1;
            s=0; // index into x
            do {
                while(--z>=0) s^=(1l<<(c[z]=c[z+1]+1));
                // m_ijk = c_ijk + c_i*c_jk + c_j*c_ik + c_k*c_ij + c_i*c_j*c_k
                x[s]=psum(x,s);
                // increment b indices
                for(z=0;z<b&&(s^=3l<<c[z],++c[z]==N-z);++z) s^=1l<<c[z];
            } while(z<b);
        };
        x[0]=exp(x[0]); // to get the correct answer if not normalized
        delete[] c;
    }

    // convert (generalized) weights to log probability or energy
    void w2logp (int N, double *x) {
        for(long s=1l<<N;s>>=1;) for(long j=1l<<N;(j-=(s<<1))>=0;)
            for(long k=j+s;--k>=j;) { x[k+s]+=x[k]; x[k]=2*x[k]-x[k+s]; };
    }

    // convert log probability or energy to (generalized) weights
    void logp2w (int N, double *x) {
        for(long s=1l<<N;s>>=1;) for(long j=1l<<N;(j-=(s<<1))>=0;)
            for(long k=j+s;--k>=j;) { x[k]=(x[k]+x[k+s])/2; x[k+s]-=x[k]; };
    }

    // convert probability to moments
    void p2m (int N, double *x) {
        for(long s=1l<<N;s>>=1;) for(long j=1l<<N;(j-=(s<<1))>=0;)
            for(long k=j+s;--k>=j;) { x[k+s]-=x[k]; x[k]=2*x[k]+x[k+s]; };
    }

    // convert moments to probability
    void m2p (int N, double *x) {
        for(long s=1l<<N;s>>=1;) for(long j=1l<<N;(j-=(s<<1))>=0;)
            for(long k=j+s;--k>=j;) { x[k]=(x[k]-x[k+s])/2; x[k+s]+=x[k]; };
    }

    // convert log probability to probability
    void logp2p (int N, double *x) {
        for(long s=1l<<N;s--;) x[s]=exp(x[s]);
    }

    // convert probability to log probability
    void p2logp (int N, double *x) {
        for(long s=1l<<N;s--;) x[s]=log(x[s]);
    }

    // normalize a log probability table
    void logpnorm (int N, double *x) {
        double f=x[0];
        for(long s=1l<<N;--s;)
            if(f>x[s]) f+=log1p(exp(x[s]-f));
            else f=x[s]+log1p(exp(f-x[s]));
        for(long s=1l<<N;s--;) x[s]-=f;
    }

    // normalize a probability table, use logpnorm whenever possible
    void pnorm (int N, double *x) {
        double z=0;
        for(long s=1l<<N;s--;) z+=x[s];
        for(long s=1l<<N;s--;) x[s]/=z;
    }

    // fills table with v for all entries with more than k indices
    // used for example when cumulants or moments are converted upto some order
    void fill (int N, double *x, int k, double v) {
        if(k>N) return;
        long ss=0,s;
        int n=0;
        for(long i=1l<<N;i--;) { // make use of gray code to count number of bits
            if(n>k) x[ss]=v;
            s=i^(i>>1);
            if(s>ss) ++n; else --n;
            ss=s;
        };
    }

}
