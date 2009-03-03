#ifndef ___defined_libdai_bbp_h
#define ___defined_libdai_bbp_h

#include <vector>
#include <utility>
#include <ext/hash_map>

#include <boost/tuple/tuple.hpp>

#include <dai/prob.h>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>

#include <dai/bp_dual.h>

namespace dai {

using namespace std;
using namespace __gnu_cxx;
using boost::tuple;
using boost::get;

// utility function to subtract 1 from a size_t, or return 0 if the
// argument is 0
inline size_t oneLess(size_t v) { return v==0?v:v-1; }

/// function to compute adj_w_unnorm from w, Z_w, adj_w
Prob unnormAdjoint(const Prob &w, Real Z_w, const Prob &adj_w);

vector<size_t> getGibbsState(const BP_dual& fg, size_t iters);

vector<Prob> get_zero_adj_2(const BP_dual& bp_dual);
vector<Prob> get_zero_adj_1(const BP_dual& bp_dual);

class BBP {
  // ----------------------------------------------------------------
  // inputs
  const BP_dual* _bp_dual;
  const FactorGraph* _fg;
  vector<Prob> _adj_b_1, _adj_b_2;
  vector<Prob> _adj_b_1_unnorm, _adj_b_2_unnorm;
  // in case the objective function depends on factors as well:
  vector<Prob> _init_adj_psi_1;
  vector<Prob> _init_adj_psi_2;
  // ----------------------------------------------------------------
  // ouptuts
  vector<Prob> _adj_psi_1, _adj_psi_2;
  // the following vectors are length _fg->nr_edges()
  vector<Prob> _adj_n, _adj_m; 
  vector<Prob> _adj_n_unnorm, _adj_m_unnorm;
  vector<Prob> _new_adj_n, _new_adj_m;
  // ----------------------------------------------------------------
  // auxiliary data
  typedef vector<hash_map<size_t, vector<size_t> > > _trip_t;
  typedef vector<size_t> _trip_end_t;

  // members to store indices of VFV and FVF triples
   // XXX need to have vector going in other direction so that we can
   // do faster looping on iIj and IiJ
  _trip_t _VFV;
  _trip_t _FVF;

  typedef vector<tuple<size_t,size_t,size_t> > _trip_ind_t;
  _trip_ind_t _VFV_ind;
  _trip_ind_t _FVF_ind;
  
  // helper quantities computed from the BP messages
  vector<Prob> _T,_U,_S,_R;
  size_t _iters;
  size_t _max_iter;
  
 public:
  void RegenerateInds();
  void RegenerateT();
  void RegenerateU();
  void RegenerateS();
  void RegenerateR();
  void RegenerateInputs();
  /// initialise members for messages and factor adjoints
  void RegeneratePsiAdjoints();
  void RegenerateMessageAdjoints();
  void Regenerate();

  size_t VFV(size_t i, size_t I, size_t j) const { return (*const_cast<_trip_t*>(&_VFV))[i][I][j]; }
  size_t FVF(size_t I, size_t i, size_t J) const { return (*const_cast<_trip_t*>(&_FVF))[I][i][J]; }
  DAI_ACCMUT(Prob & T(size_t i, size_t I), { return _T[_fg->VV2E(i,I)]; });
  DAI_ACCMUT(Prob & U(size_t I, size_t i), { return _U[_fg->VV2E(i,I)]; });
  DAI_ACCMUT(Prob & S(size_t i, size_t I, size_t j), { return _S[VFV(i,I,j)]; });
  DAI_ACCMUT(Prob & R(size_t I, size_t i, size_t J), { return _R[FVF(I,i,J)]; });

  DAI_ACCMUT(Prob & T(size_t iI), { return _T[iI]; });
  DAI_ACCMUT(Prob & U(size_t iI), { return _U[iI]; });
  DAI_ACCMUT(Prob & S(size_t iIj), { return _S[iIj]; });
  DAI_ACCMUT(Prob & R(size_t IiJ), { return _R[IiJ]; });
  
  void calcNewN(size_t iI);
  void calcNewM(size_t iI);
  void upMsgM(size_t iI);
  void upMsgN(size_t iI);
  void doParUpdate();
  Real getUnMsgMag();
  void getMsgMags(Real &s, Real &new_s);

  void zero_adj_b_2() {
    _adj_b_2.clear();
    _adj_b_2.reserve(_bp_dual->nrFactors());
    for(size_t I=0; I<_bp_dual->nrFactors(); I++) {
      _adj_b_2.push_back(Prob(_bp_dual->factor(I).states(),Real(0.0)));
    }
  }
 public:
  BBP(const BP_dual &bp_dual) :
    _bp_dual(&bp_dual), _fg(&bp_dual), _max_iter(50) {}

  DAI_ACCMUT(size_t& maxIter(), { return _max_iter; });

  void init(const vector<Prob> &adj_b_1, const vector<Prob> &adj_b_2,
            const vector<Prob> &adj_psi_1, const vector<Prob> &adj_psi_2) {
    _adj_b_1 = adj_b_1;
    _adj_b_2 = adj_b_2;
    _init_adj_psi_1 = adj_psi_1;
    _init_adj_psi_2 = adj_psi_2;
    Regenerate(); 
  }
  void init(const vector<Prob> &adj_b_1, const vector<Prob> &adj_b_2) {
    init(adj_b_1, adj_b_2, get_zero_adj_1(*_bp_dual), get_zero_adj_2(*_bp_dual));
  }
  void init(const vector<Prob> &adj_b_1) {
    init(adj_b_1, get_zero_adj_2(*_bp_dual));
  }

  /// run until change is less than given tolerance
  void run(Real tol);

  size_t doneIters() { return _iters; }

  DAI_ACCMUT(Prob& adj_psi_1(size_t i), { return _adj_psi_1[i]; });
  DAI_ACCMUT(Prob& adj_psi_2(size_t I), { return _adj_psi_2[I]; });
  DAI_ACCMUT(Prob& adj_b_1(size_t i), { return _adj_b_1[i]; });
  DAI_ACCMUT(Prob& adj_b_2(size_t I), { return _adj_b_2[I]; });
  DAI_ACCMUT(Prob& adj_n(size_t i, size_t I), { return _adj_n[_fg->VV2E(i,I)]; });
  DAI_ACCMUT(Prob& adj_m(size_t I, size_t i), { return _adj_m[_fg->VV2E(i,I)]; });
};

/// Cost functions. Not used by BBP class, only used by following functions
DAI_ENUM(bbp_cfn_t, cfn_b, cfn_b2, cfn_exp, cfn_var_ent, cfn_factor_ent, cfn_bethe_ent);

/// initialise BBP using BP_dual and (cfn_type, stateP)
void gibbsInitBBPCostFnAdj(BBP& bbp, const BP_dual& fg, bbp_cfn_t cfn_type, const vector<size_t>* stateP);

/// calculate actual value of cost function (cfn_type, stateP)
Real gibbsCostFn(const BP_dual& fg, bbp_cfn_t cfn_type, const vector<size_t> *stateP);

/// function to test the validity of adjoints computed by BBP
double numericBBPTest(const BP_dual& bp_dual, bbp_cfn_t cfn, double bbpTol, double h);

/// output to "./bbp-data/" a set of files which evaluate the effect
/// on the given cost function (for a random Gibbs state) of
/// perturbing each graph variable
void makeBBPGraph(const BP_dual& bp_dual, bbp_cfn_t cfn);

}

#endif
