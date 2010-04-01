/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2009  Patrick Pletscher  [pletscher at inf dot ethz dot ch]
 *                2009  Joris Mooij        [joris dot mooij at libdai dot org]
 */


%module dai

        struct Neighbor {
            size_t iter;
            size_t node;
            size_t dual;

            Neighbor() {}
            Neighbor( size_t iter, size_t node, size_t dual ) : iter(iter), node(node), dual(dual) {}

            operator size_t () const { return node; }
        };

%{
#include "../include/dai/var.h"
#include "../include/dai/smallset.h"
#include "../include/dai/varset.h"
#include "../include/dai/prob.h"
#include "../include/dai/factor.h"
#include "../include/dai/bipgraph.h"
#include "../include/dai/factorgraph.h"
#include "../include/dai/util.h"
%}

%ignore dai::TProb::operator[];
%ignore dai::TFactor::operator[];

%ignore dai::Var::label() const;
%ignore dai::Var::states() const;

%include "../include/dai/util.h"
%include "../include/dai/var.h"
%include "../include/dai/smallset.h"
%template(SmallSetVar) dai::SmallSet< dai::Var >;
%include "../include/dai/varset.h"
%extend dai::VarSet {
        inline void append(const dai::Var &v) { (*self) |= v; }   /* for python, octave */
};

%include "../include/dai/prob.h"
%template(Prob) dai::TProb<dai::Real>;
%extend dai::TProb<dai::Real> {
        inline dai::Real __getitem__(int i) const {return (*self).get(i);} /* for python */
        inline void __setitem__(int i,dai::Real d) {(*self).set(i,d);}   /* for python */
        inline dai::Real __paren(int i) const {return (*self).get(i);}     /* for octave */
        inline void __paren_asgn(int i,dai::Real d) {(*self).set(i,d);}  /* for octave */
};
%include "../include/dai/factor.h"
%extend dai::TFactor<dai::Real> {
        inline dai::Real __getitem__(int i) const {return (*self).get(i);} /* for python */
        inline void __setitem__(int i,dai::Real d) {(*self).set(i,d);}   /* for python */
        inline dai::Real __paren__(int i) const {return (*self).get(i);}     /* for octave */
        inline void __paren_asgn__(int i,dai::Real d) {(*self).set(i,d);}  /* for octave */
};

%template(Factor) dai::TFactor<dai::Real>;
%include "../include/dai/bipgraph.h"
%include "../include/dai/factorgraph.h"
%include "std_vector.i"
// TODO: typemaps for the vectors (input/output python arrays)
%inline{
typedef std::vector<dai::Factor> VecFactor;
typedef std::vector< VecFactor > VecVecFactor;
}
%template(VecFactor) std::vector< dai::Factor >;
%template(VecVecFactor) std::vector< VecFactor >;

%{
typedef dai::BipartiteGraph::Neighbor Neighbor;
%}

%include "../include/dai/index.h"
%extend dai::multifor {
    inline size_t __getitem__(int i) const {
        return (*self)[i];
    }
    inline void next() {
        return (*self)++;
    }
};
