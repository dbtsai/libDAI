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


#ifndef __defined_libdai_mr_h
#define __defined_libdai_mr_h


#include <vector>
#include <string>
#include <dai/factorgraph.h>
#include <dai/daialg.h>
#include <dai/enum.h>
#include <dai/properties.h>


namespace dai {


class sub_nb;


class MR : public DAIAlgFG {
    private:
        bool supported;                                            // is the underlying factor graph supported?

        std::vector<size_t>                             con;       // con[i] = connectivity of spin i
        std::vector<std::vector<size_t> >               nb;        // nb[i] are the neighbours of spin i
        std::vector<std::vector<double> >               tJ;        // tJ[i][_j] is the tanh of the interaction between spin i and its neighbour nb[i][_j]
        std::vector<double>                             theta;     // theta[i] is the local field on spin i
        std::vector<std::vector<double> >               M;         // M[i][_j] is M^{(i)}_j
        std::vector<std::vector<size_t> >               kindex;    // the _j'th neighbour of spin i has spin i as its kindex[i][_j]'th neighbour
        std::vector<std::vector<std::vector<double> > > cors;
    
        static const size_t kmax = 31;
        
        size_t N;

        std::vector<double> Mag;

    public:
        struct Properties {
            size_t verbose;
            double tol;
            ENUM2(UpdateType,FULL,LINEAR)
            ENUM3(InitType,RESPPROP,CLAMPING,EXACT)
            UpdateType updates;
            InitType inits;
        } props;
        double maxdiff;

    public:
        MR( const FactorGraph & fg, const PropertySet &opts );
        void init(size_t Nin, double *_w, double *_th);
        void makekindex();
        void read_files();
        void init_cor();
        double init_cor_resp();
        void solvemcav();
        void solveM();
        double run();
        Factor belief( const Var &n ) const;
        Factor belief( const VarSet &/*ns*/ ) const { assert( 0 == 1 ); }
        std::vector<Factor> beliefs() const;
        Complex logZ() const { return NAN; }
        void init() {}
        static const char *Name;
        std::string identify() const;
        double _tJ(size_t i, sub_nb A);

        double Omega(size_t i, size_t _j, size_t _l);
        double T(size_t i, sub_nb A);
        double T(size_t i, size_t _j);
        double Gamma(size_t i, size_t _j, size_t _l1, size_t _l2);
        double Gamma(size_t i, size_t _l1, size_t _l2);

        double appM(size_t i, sub_nb A);
        void sum_subs(size_t j, sub_nb A, double *sum_even, double *sum_odd);

        double sign(double a) { return (a >= 0) ? 1.0 : -1.0; }
        MR* clone() const { assert( 0 == 1 ); }

        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        double maxDiff() const { return maxdiff; }
}; 


// represents a subset of nb[i] as a binary number
// the elements of a subset should be thought of as indices in nb[i]
class sub_nb {
    private:
        size_t s;
        size_t bits;
    
    public:
        // construct full subset containing nr_elmt elements
        sub_nb(size_t nr_elmt) {
#ifdef DAI_DEBUG
            assert( nr_elmt < sizeof(size_t) / sizeof(char) * 8 );
#endif
            bits = nr_elmt;
            s = (1U << bits) - 1;
        }

        // copy constructor
        sub_nb( const sub_nb & x ) : s(x.s), bits(x.bits) {}

        // assignment operator 
        sub_nb & operator=( const sub_nb &x ) {
            if( this != &x ) {
                s = x.s;
                bits = x.bits;
            }
            return *this;
        }

        // returns number of elements
        size_t size() {
            size_t size = 0;
            for(size_t bit = 0; bit < bits; bit++)
                if( s & (1U << bit) )
                    size++;
            return size;
        }

        // increases s by one (for enumeration in lexicographical order)
        sub_nb operator++() { 
            s++; 
            if( s >= (1U << bits) )
                s = 0;
            return *this; 
        }
        
        // return i'th element of this subset
        size_t operator[](size_t i) { 
            size_t bit;
            for(bit = 0; bit < bits; bit++ )
                if( s & (1U << bit) ) {
                    if( i == 0 )
                        break;
                    else
                        i--;
                }
#ifdef DAI_DEBUG
            assert( bit < bits );
#endif
            return bit;
        }

        // add index _j to this subset
        sub_nb &operator +=(size_t _j) {
            s |= (1U << _j); 
            return *this;
        }

        // return copy with index _j
        sub_nb operator+(size_t _j) {
            sub_nb x = *this;
            x += _j;
            return x;
        }

        // delete index _j from this subset
        sub_nb &operator -=(size_t _j) {
            s &= ~(1U << _j); 
            return *this;
        }

        // return copy without index _j
        sub_nb operator-(size_t _j) {
            sub_nb x = *this;
            x -= _j;
            return x;
        }

        // empty this subset
        sub_nb & clear() {
            s = 0;
            return *this;
        }

        // returns true if subset is empty
        bool empty() { return (s == 0); }

        // return 1 if _j is contained, 0 otherwise ("is element of")
        size_t operator&(size_t _j) { return s & (1U << _j); }

        friend std::ostream& operator<< (std::ostream& os, const sub_nb x) {
            if( x.bits == 0 )
                os << "empty";
            else {
                for(size_t bit = x.bits; bit > 0; bit-- )
                    if( x.s & (1U << (bit-1)) )
                        os << "1";
                    else
                        os << "0";
            }
            return os;
        }
};


} // end of namespace dai


#endif
