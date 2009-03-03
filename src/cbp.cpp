#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>

#include <dai/util.h>
#include <dai/properties.h>

#include <dai/bp.h>
#include <dai/cbp.h>
#include <dai/bbp.h>

using namespace std;

namespace dai {

const char *CBP::Name = "CBP";

const char *CBP::PropertyList[] = {"updates","tol","rec_tol","maxiter","verbose","max_levels","min_max_adj","choose","clamp","recursion","bbp_cfn","rand_seed"};

#define rnd_multi(x) rnd_int(0,(x)-1)

void CBP::setProperties(const PropertySet &opts) {
//     DAI_DMSG("in CBP::setProperties");
//     DAI_PV(opts);
    foreach(const char* p, PropertyList) {
        if(!opts.hasKey(p)) {
            // XXX probably leaks pointer?
            throw (string("CBP: Missing property ")+p).c_str();
        }
    }
    
    props.tol = opts.getStringAs<double>("tol");
    props.rec_tol = opts.getStringAs<double>("rec_tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.verbose = opts.getStringAs<size_t>("verbose");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");
    props.max_levels = opts.getStringAs<size_t>("max_levels");
    props.min_max_adj = opts.getStringAs<double>("min_max_adj");
    props.choose = opts.getStringAs<Properties::ChooseMethodType>("choose");
    props.recursion = opts.getStringAs<Properties::RecurseType>("recursion");
    props.clamp = opts.getStringAs<Properties::ClampType>("clamp");
    props.bbp_cfn = opts.getStringAs<bbp_cfn_t>("bbp_cfn");
    props.rand_seed = opts.getStringAs<size_t>("rand_seed");
}

PropertySet CBP::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "rec_tol", props.rec_tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "updates", props.updates );
    opts.Set( "max_levels", props.max_levels );
    opts.Set( "min_max_adj", props.min_max_adj );
    opts.Set( "choose", props.choose );
    opts.Set( "recursion", props.recursion );
    opts.Set( "clamp", props.clamp );
    opts.Set( "bbp_cfn", props.bbp_cfn );
    opts.Set( "rand_seed", props.rand_seed );
    return opts;
}

std::string CBP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "rec_tol=" << props.rec_tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "updates=" << props.updates << ",";
    s << "max_levels=" << props.max_levels << ",";
    s << "min_max_adj=" << props.min_max_adj << ",";
    s << "choose=" << props.choose << ",";
    s << "recursion=" << props.recursion << ",";
    s << "clamp=" << props.clamp << ",";
    s << "bbp_cfn=" << props.bbp_cfn << ",";
    s << "rand_seed=" << props.rand_seed << ",";
    s << "]";
    return s.str();
}

void CBP::construct() {
//     DAIAlgFG::Regenerate();
    indexEdges();

    _beliefs1.clear(); _beliefs1.reserve(nrVars());
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefs1.push_back( Factor(var(i)).normalized() );

    _beliefs2.clear(); _beliefs2.reserve(nrFactors());
    for( size_t I = 0; I < nrFactors(); I++ ) {
        Factor f = factor(I);
        f.fill(1); f.normalize();
        _beliefs2.push_back(f);
    }

    // to compute average level
    _sum_level = 0;
    _num_leaves = 0;

    _maxdiff = 0;
    _iters = 0;
}

static
vector<Factor> mixBeliefs(Real p, vector<Factor> b, vector<Factor> c) {
  vector<Factor> out;
  assert(b.size()==c.size());
  out.reserve(b.size());
  Real pc = 1-p;
  for(size_t i=0; i<b.size(); i++) {
      // XXX probably already normalized
    out.push_back(b[i].normalized()*p+
                  c[i].normalized()*pc);
  }
  return out;
}

double CBP::run() {
//     BP bp(getBP());
//     InfAlg *bp = newInfAlg( GetPropertyAs<string>("bp_alg"), *this, GetPropertyAs<Properties>("bp_opts") );
    size_t seed = props.rand_seed;
    if(seed>0) rnd_seed(seed);

    BP_dual bp_dual(getBP_dual());
    bp_dual.init();
    bp_dual.run();
    _iters += bp_dual.Iterations();

    vector<Factor> beliefs_out;
    Real lz_out;
    runRecurse(bp_dual, bp_dual.logZ(), 
               vector<size_t>(0), set<size_t>(),
               _num_leaves, _sum_level,
               lz_out, beliefs_out);
    if(props.verbose>=1) {
      cerr << "CBP average levels = " << (_sum_level/_num_leaves) << ", leaves = " << _num_leaves << endl;
    }
    setBeliefs(beliefs_out, lz_out);
    return 0;
}

BP CBP::getBP() {
    PropertySet bpProps;
    bpProps.Set("updates", string("PARALL"));
    bpProps.Set("tol", props.tol);
    bpProps.Set("maxiter", props.maxiter);
    bpProps.Set("verbose", oneLess(props.verbose));
    BP bp(*this,bpProps);
    bp.init();
    return bp;
}

BP_dual CBP::getBP_dual() {
    PropertySet bpProps;
    bpProps.Set("updates", string("PARALL"));
    bpProps.Set("tol", props.tol);
    bpProps.Set("maxiter", props.maxiter);
    bpProps.Set("verbose", oneLess(props.verbose));
//     cerr << "In getBP_dual" << endl;
//     DAI_PV(bpProps);
    BP_dual bp_dual(*this,bpProps);
    return bp_dual;
}

vector<size_t> complement(vector<size_t>& xis, size_t n_states) {
    vector<size_t> cmp_xis(0);
    size_t j=0;
    for(size_t xi=0; xi<n_states; xi++) {
        while(j<xis.size() && xis[j]<xi) j++;
        if(j>=xis.size() || xis[j]>xi) cmp_xis.push_back(xi);
    }
    assert( xis.size()+cmp_xis.size() == n_states );
    return cmp_xis;
}

Real max(Real x,Real y) { return x>y?x:y; }

Real unSoftMax(Real lz, Real cmp_lz) {
    double m = max(lz, cmp_lz);
    lz -= m; cmp_lz -= m;
    double p = exp(lz)/(exp(lz)+exp(cmp_lz));
    return p;
}

Real logSumExp(Real lz, Real cmp_lz) {
    double m = max(lz, cmp_lz);
    lz -= m; cmp_lz -= m;
    return m+log(exp(lz)+exp(cmp_lz));
}

Real dist(const vector<Factor>& b1, const vector<Factor>& b2, size_t nv) {
    Real d=0.0;
    for(size_t k=0; k<nv; k++) {
        d += dist( b1[k], b2[k], Prob::DISTLINF );
    }
    return d;
}


void CBP::runRecurse(BP_dual &bp_dual,
                     double orig_logZ,
                     vector<size_t> clamped_vars_list,
                     set<size_t> clamped_vars,
                     size_t &num_leaves,
                     double &sum_level,
                     Real &lz_out, vector<Factor>& beliefs_out) {
    // choose a variable/states to clamp:
    size_t i;
    vector<size_t> xis;
    Real maxVar=0.0;
    bool found;
    bool clampVar = (Clamping()==Properties::ClampType::CLAMP_VAR);

    // XXX fix to just pass orig_logZ

    if(Recursion()==Properties::RecurseType::REC_LOGZ && recTol()>0 &&
       exp(bp_dual.logZ()-orig_logZ) < recTol()) {
        found = false;
    } else {
        found = chooseNextClampVar(bp_dual,
                                   clamped_vars_list,
                                   clamped_vars,
                                   i, xis, &maxVar);
    }

    if(!found) {
        num_leaves++;
        sum_level += clamped_vars_list.size();
        beliefs_out = bp_dual.beliefs();
        lz_out = bp_dual.logZ();
        return;
    }
    
    if(clampVar) {
        foreach(size_t xi, xis) { assert(/*0<=xi &&*/ xi<var(i).states()); }
    } else {
        foreach(size_t xI, xis) { assert(/*0<=xI &&*/ xI<factor(i).states()); }
    }
    // - otherwise, clamp and recurse, saving margin estimates for each
    // clamp setting. afterwards, combine estimates.

    // compute complement of 'xis'
    vector<size_t> cmp_xis=complement(xis, clampVar?var(i).states():factor(i).states());

    // XXX could do this more efficiently with a nesting version of
    // saveProbs/undoProbs
    Real lz; vector<Factor> b;
    BP_dual bp_dual_c(bp_dual);
    if(clampVar) {
        _clamp((FactorGraph&)bp_dual_c, var(i), xis);
        bp_dual_c.init(var(i));
    } else {
        _clampFactor((FactorGraph&)bp_dual_c, i, xis);
        bp_dual_c.init(factor(i).vars());
    }
    bp_dual_c.run();
    _iters += bp_dual_c.Iterations();

    lz = bp_dual_c.logZ();
    b = bp_dual_c.beliefs();

    Real cmp_lz; vector<Factor> cmp_b;
    BP_dual cmp_bp_dual_c(bp_dual);
    if(clampVar) {
        _clamp(cmp_bp_dual_c,var(i),cmp_xis);
        cmp_bp_dual_c.init(var(i));
    } else {
        _clampFactor(cmp_bp_dual_c,i,cmp_xis);
        cmp_bp_dual_c.init(factor(i).vars());
    }
    cmp_bp_dual_c.run();
    _iters += cmp_bp_dual_c.Iterations();

    cmp_lz = cmp_bp_dual_c.logZ();
    cmp_b = cmp_bp_dual_c.beliefs();

    double p = unSoftMax(lz, cmp_lz);
    Real bp__d=0.0;
    
    if(Recursion()==Properties::RecurseType::REC_BDIFF && recTol() > 0) {
        vector<Factor> combined_b(mixBeliefs(p,b,cmp_b));
        Real new_lz = logSumExp(lz,cmp_lz);
        bp__d = dist(bp_dual.beliefs(),combined_b,nrVars());
        if(exp(new_lz-orig_logZ)*bp__d < recTol()) {
            num_leaves++;
            sum_level += clamped_vars_list.size();
            beliefs_out = combined_b;
            lz_out = new_lz;
            return;
        }
    }

    // either we are not doing REC_BDIFF or the distance was large
    // enough to recurse:

    runRecurse(bp_dual_c, orig_logZ, 
               clamped_vars_list,
               clamped_vars,
               num_leaves, sum_level, lz, b);
    runRecurse(cmp_bp_dual_c, orig_logZ, 
               clamped_vars_list,
               clamped_vars,
               num_leaves, sum_level, cmp_lz, cmp_b);

    p = unSoftMax(lz,cmp_lz);

    beliefs_out = mixBeliefs(p, b, cmp_b);
    lz_out = logSumExp(lz,cmp_lz);

    if(props.verbose>=2) {
        Real d = dist( bp_dual.beliefs(), beliefs_out, nrVars() );
        cerr << "Distance (clamping " << i << "): " << d;
        if(Recursion()==Properties::RecurseType::REC_BDIFF)
            cerr << "; bp_dual predicted " << bp__d;
        cerr << "; max adjoint = " << maxVar << "; level = " << clamped_vars_list.size() << endl;
    }
}

// 'xis' must be sorted
bool CBP::chooseNextClampVar(BP_dual& bp,
                             vector<size_t> &clamped_vars_list,
                             set<size_t> &clamped_vars,
                             size_t &i, vector<size_t> &xis, Real *maxVarOut) {
    Real tiny=1.0e-14;
    if(props.verbose>=3) {
        cerr << "clamped_vars_list" << clamped_vars_list << endl;
    }
    if(clamped_vars_list.size() >= maxClampLevel()) {
        return false;
    }
    if(ChooseMethod()==Properties::ChooseMethodType::CHOOSE_RANDOM) {
        if(Clamping()==Properties::ClampType::CLAMP_VAR) {
            int t=0, t1=100;
            do {
                i = rnd_multi(nrVars());
                t++;
            } while(abs(bp.belief1(i).p().max()-1) < tiny &&
                    t < t1);
            if(t==t1) {
                return false;
                //             die("Too many levels requested in CBP");
            }
            // only pick probable values for variable
            size_t xi;
            do {
                xi = rnd_multi(var(i).states());
                t++;
            } while(bp.belief1(i).p()[xi] < tiny && t<t1);
            assert(t<t1);
            xis.resize(1, xi);
            //         assert(!_clamped_vars.count(i)); // not true for >2-ary variables
            DAI_IFVERB(2, endl<<"CHOOSE_RANDOM chose variable "<<i<<" state "<<xis[0]<<endl);
        } else {
            int t=0, t1=100;
            do {
                i = rnd_multi(nrFactors());
                t++;
            } while(abs(bp.belief2(i).p().max()-1) < tiny &&
                    t < t1);
            if(t==t1) {
                return false;
                //             die("Too many levels requested in CBP");
            }
            // only pick probable values for variable
            size_t xi;
            do {
                xi = rnd_multi(factor(i).states());
                t++;
            } while(bp.belief2(i).p()[xi] < tiny && t<t1);
            assert(t<t1);
            xis.resize(1, xi);
            //         assert(!_clamped_vars.count(i)); // not true for >2-ary variables
            DAI_IFVERB(2, endl<<"CHOOSE_RANDOM chose factor "<<i<<" state "<<xis[0]<<endl);
        }
    } else if(ChooseMethod()==Properties::ChooseMethodType::CHOOSE_BP_ALL) {
      // try clamping each variable manually
      assert(Clamping()==Properties::ClampType::CLAMP_VAR);
      Real max_diff=0.0;
      int win_k=-1, win_xk=-1;
      for(size_t k=0; k<nrVars(); k++) {
        for(size_t xk=0; xk<var(k).states(); xk++) {
          if(bp.belief1(k)[xk]<tiny) continue;
          BP_dual bp1(bp);
          bp1.clamp(var(k), xk);
          bp1.init(var(k));
          bp1.run();
          Real diff=0;
          for(size_t j=0; j<nrVars(); j++) {
            diff += dist(bp.belief1(j), bp1.belief1(j), Prob::DISTL1);
          }
          if(diff>max_diff) {
            max_diff=diff; win_k=k; win_xk=xk;
          }
        }
      }
      assert(win_k>=0); assert(win_xk>=0);
      i = win_k; xis.resize(1, win_xk);
    } else if(ChooseMethod()==Properties::ChooseMethodType::CHOOSE_BBP) {
        Real mvo; if(!maxVarOut) maxVarOut = &mvo;
        bool clampVar = (Clamping()==Properties::ClampType::CLAMP_VAR);
        pair<size_t, size_t> cv =
            bbpFindClampVar(bp,
                            clampVar,
                            clamped_vars_list.size(),
                            BBP_cost_function(),getProperties(),maxVarOut);

        // if slope isn't big enough then don't clamp
        if(*maxVarOut < minMaxAdj()) return false;

        size_t xi=cv.second;
        i = cv.first;
#define VAR_INFO (clampVar?"variable ":"factor ")                       \
            << i << " state " << xi                                     \
            << " (p=" << (clampVar?bp.belief1(i)[xi]:bp.belief2(i)[xi]) \
            << ", entropy = " << (clampVar?bp.belief1(i):bp.belief2(i)).entropy() \
            << ", maxVar = "<< *maxVarOut << ")" 
        Prob b = (clampVar?bp.belief1(i).p():bp.belief2(i).p());
        if(b[xi] < tiny) {
            cerr << "Warning, bbpFindClampVar found unlikely "
                 << VAR_INFO << endl;
            return false;
        }
        if(abs(b[xi]-1) < tiny) {
            cerr << "Warning, bbpFindClampVar found overly likely "
                 << VAR_INFO << endl;
            return false;
        }

        xis.resize(1,xi);
        if(clampVar) {
            assert(/*0<=xi &&*/ xi<var(i).states());
        } else {
            assert(/*0<=xi &&*/ xi<factor(i).states());
        }
        DAI_IFVERB(2, "CHOOSE_BBP (num clamped = " << clamped_vars_list.size()
               << ") chose " << i << " state " << xi << endl);
    } else {
        abort();
    }
    clamped_vars_list.push_back(i);
    clamped_vars.insert(i);
    return true;
}

// void CBP::clamp(const Var & n, size_t i) {
//   FactorGraph::clamp(n,i);
//   _clamped_vars.insert(findVar(n));
//   _clamped_vars_list.push_back(findVar(n));
// }

// void CBP::clamp(const Var & n, const vector<size_t> &is) {
//   FactorGraph::clamp(n,is);
//   _clamped_vars.insert(findVar(n));
//   _clamped_vars_list.push_back(findVar(n));
// }

void CBP::printDebugInfo() {
    DAI_PV(_beliefs1);
    DAI_PV(_beliefs2);
    DAI_PV(_logZ);
}

//----------------------------------------------------------------

bool doBBPTest=false;
bool doBBPGraph=false;
size_t bbpGraphLevel=3;

#define BPP_INIT_GIBBS 1

/// function which takes a factor graph as input, runs Gibbs and BP_dual,
/// creates and runs a BBP object, finds best variable, returns
/// (variable,state) pair for clamping
// pair<size_t, size_t> bbpFindClampVar(const CBP &fg, bbp_cfn_t cfn, const Properties &props, Real *maxVarOut) {
pair<size_t, size_t> bbpFindClampVar(BP_dual &in_bp_dual, bool clampVar,
    size_t numClamped, bbp_cfn_t cfn, const PropertySet &props, Real *maxVarOut) {
#if BPP_INIT_GIBBS
    vector<size_t> state = getGibbsState(in_bp_dual, 100);
    in_bp_dual.init(state);
    in_bp_dual.run();
#endif
  
    Real ourTol = doBBPTest ? 1.0e-11 : 1.0e-3;
    if(0) {
        PropertySet bp_Props;
        bp_Props.Set("updates", string("PARALL"));
        //   bp_Props.Set("tol", props.GetAs<double>("tol"));
        bp_Props.Set("tol", ourTol);
        bp_Props.Set("maxiter", props.GetAs<size_t>("maxiter"));
        bp_Props.Set("verbose", oneLess(props.GetAs<size_t>("verbose")));
        //   bp_Props.ConvertTo<BP_dual::UpdateType>("updates");
        //   DAI_PV(bp_Props.GetAs<BP_dual::UpdateType>("updates"));
        BP_dual bp_dual(in_bp_dual, bp_Props);
#if BPP_INIT_GIBBS
        bp_dual.init(state);
#endif
        bp_dual.run();
    }

    if(doBBPGraph && numClamped == bbpGraphLevel) {
        cerr << "Writing BBP graph data" << endl;
        makeBBPGraph(in_bp_dual,cfn);
        doBBPGraph=false; // only do it once
        cerr << "Done writing BBP graph data" << endl;
    }
    if(doBBPTest) {
        double err = numericBBPTest(in_bp_dual, cfn, /*bbp tol*/ ourTol, /*h*/ 1.0e-5);
        cerr << "Error from numericBBPTest: " << err << endl;
    }
    Real tic1=toc();
    BBP bbp(in_bp_dual);
    bbp.maxIter() = props.GetAs<size_t>("maxiter");
#if BPP_INIT_GIBBS
    gibbsInitBBPCostFnAdj(bbp, in_bp_dual, cfn, &state);
#else
    gibbsInitBBPCostFnAdj(bbp, in_bp_dual, cfn, NULL);
#endif
    Real tic2=toc();
    bbp.run(ourTol);
    if(props.GetAs<size_t>("verbose") >= 3) {
        cerr << "BBP took " << toc()-tic1 << " seconds (BBP.run = " << toc()-tic2 << " seconds), "
             << bbp.doneIters() << " iterations" << endl;
    }
  
    // find and return the (variable,state) with the largest adj_psi_1
    size_t argmax_var=0;
    size_t argmax_var_state=0;
    Real max_var=0;
    if(clampVar) {
        for(size_t i=0; i<in_bp_dual.nrVars(); i++) {
            Prob adj_psi_1 = bbp.adj_psi_1(i);
            if(0) {
                // helps to account for amount of movement possible in variable
                // i's beliefs? seems not..
                adj_psi_1 *= in_bp_dual.belief1(i).entropy();
            }
            // try to compensate for effect on same variable (doesn't work)
            //     adj_psi_1[gibbs.state()[i]] -= bp_dual.belief1(i)[gibbs.state()[i]]/10;
            pair<size_t,Real> argmax_state = adj_psi_1.argmax();

            if(i==0 || argmax_state.second>max_var) {
                argmax_var = i;
                max_var = argmax_state.second;
                argmax_var_state = argmax_state.first;
            }
        }
        assert(/*0 <= argmax_var_state &&*/
               argmax_var_state < in_bp_dual.var(argmax_var).states());
    } else {
        for(size_t I=0; I<in_bp_dual.nrFactors(); I++) {
            Prob adj_psi_2 = bbp.adj_psi_2(I);
            if(0) {
                // helps to account for amount of movement possible in variable
                // i's beliefs? seems not..
                adj_psi_2 *= in_bp_dual.belief2(I).entropy();
            }
            // try to compensate for effect on same variable (doesn't work)
            //     adj_psi_1[gibbs.state()[i]] -= bp_dual.belief1(i)[gibbs.state()[i]]/10;
            pair<size_t,Real> argmax_state = adj_psi_2.argmax();

            if(I==0 || argmax_state.second>max_var) {
                argmax_var = I;
                max_var = argmax_state.second;
                argmax_var_state = argmax_state.first;
            }
        }
        assert(/*0 <= argmax_var_state &&*/
               argmax_var_state < in_bp_dual.factor(argmax_var).states());
    }
    if(maxVarOut) *maxVarOut = max_var;
    return make_pair(argmax_var,argmax_var_state);
}

} // end of namespace dai
