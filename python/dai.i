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

%include "../include/dai/util.h"
%include "../include/dai/var.h"
%include "../include/dai/smallset.h"
%template(SmallSetVar) dai::SmallSet< dai::Var >;
%include "../include/dai/varset.h"
%include "../include/dai/prob.h"
%template(Prob) dai::TProb<Real>;
%extend dai::TProb<Real> {
        inline Real __getitem__(int i) const {return (*self)[i];}
        inline void __setitem__(int i,Real d) {(*self)[i] = d;}
        %template(TProbRealConstructor) TProb<double *>;
        %template(TProbIntConstructor)  TProb<size_t *>;
};
%include "../include/dai/factor.h"
%extend dai::TFactor<Real> {
        inline Real __getitem__(int i) const {return (*self)[i];}
        inline void __setitem__(int i,Real d) {(*self)[i] = d;}
};

%template(Factor) dai::TFactor<Real>;
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
%extend dai::MultiFor {
    inline size_t __getitem__(int i) const {
        return (*self)[i];
    }
    inline void next() {
        return (*self)++;
    }
};
