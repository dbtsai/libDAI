#include <dai/bp.h>
#include <dai/bbp.h>
#include <dai/gibbs.h>

// for makeBBPGraph: {
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
// }

namespace dai {

using namespace std;

#define rnd_multi(x) rnd_int(0,(x)-1)

/// function to compute \sld\~w from w, Z_w, \sld w
Prob unnormAdjoint(const Prob &w, Real Z_w, const Prob &adj_w) {
    assert(w.size()==adj_w.size());
    Prob adj_w_unnorm(w.size(),0.0);
    Real s=0.0;
    for(size_t i=0; i<w.size(); i++) {
        s += w[i]*adj_w[i];
    }
    for(size_t i=0; i<w.size(); i++) {
        adj_w_unnorm[i] = (adj_w[i]-s)/Z_w;
    }
    return adj_w_unnorm;
}

/// Function to turn Gibbs state into b1_adj
/// calls bbp.init with calculated adjoints
void gibbsInitBBPCostFnAdj(BBP& bbp, const BP_dual &fg, bbp_cfn_t cfn_type, const vector<size_t>* stateP) {
    if(cfn_type==bbp_cfn_t::bbp_cfn_t::cfn_bethe_ent) {
        vector<Prob> b1_adj;
        vector<Prob> b2_adj;
        vector<Prob> psi1_adj;
        vector<Prob> psi2_adj;
        b1_adj.reserve(fg.nrVars());
        psi1_adj.reserve(fg.nrVars());
        b2_adj.reserve(fg.nrFactors());
        psi2_adj.reserve(fg.nrFactors());
        for(size_t i=0; i<fg.nrVars(); i++) {
            size_t dim = fg.var(i).states();
            int c = fg.nbV(i).size();
            Prob p(dim,0.0);
            for(size_t xi=0; xi<dim; xi++) {
                p[xi] = (1-c)*(1+log(fg.belief1(i)[xi]));
            }
            b1_adj.push_back(p);

            for(size_t xi=0; xi<dim; xi++) {
                p[xi] = -fg.belief1(i)[xi];
            }
            psi1_adj.push_back(p);
        }
        for(size_t I=0; I<fg.nrFactors(); I++) {
            size_t dim = fg.factor(I).states();
            Prob p(dim,0.0);
            for(size_t xI=0; xI<dim; xI++) {
                p[xI] = 1 + log(fg.belief2(I)[xI]/fg.factor(I).p()[xI]);
            }
            b2_adj.push_back(p);

            for(size_t xI=0; xI<dim; xI++) {
                p[xI] = -fg.belief2(I)[xI]/fg.factor(I).p()[xI];
            }
            psi2_adj.push_back(p);
        }
        bbp.init(b1_adj, b2_adj, psi1_adj, psi2_adj);
    } else if(cfn_type==bbp_cfn_t::cfn_factor_ent) {
        vector<Prob> b2_adj;
        b2_adj.reserve(fg.nrFactors());
        for(size_t I=0; I<fg.nrFactors(); I++) {
            size_t dim = fg.factor(I).states();
            Prob p(dim,0.0);
            for(size_t xI=0; xI<dim; xI++) {
                double bIxI = fg.belief2(I)[xI];
                if(bIxI<1.0e-15) {
                    p[xI] = -1.0e10;
                } else {
                    p[xI] = 1+log(bIxI);
                }
            }
            b2_adj.push_back(p);
        }
        bbp.init(get_zero_adj_1(fg), b2_adj);
    } else if(cfn_type==bbp_cfn_t::cfn_var_ent) {
        vector<Prob> b1_adj;
        b1_adj.reserve(fg.nrVars());
        for(size_t i=0; i<fg.nrVars(); i++) {
            size_t dim = fg.var(i).states();
            Prob p(dim,0.0);
            for(size_t xi=0; xi<fg.var(i).states(); xi++) {
                double bixi = fg.belief1(i)[xi];
                if(bixi<1.0e-15) {
                    p[xi] = -1.0e10;
                } else {
                    p[xi] = 1+log(bixi);
                }
            }
            b1_adj.push_back(p);
        }
        bbp.init(b1_adj);
    } else { // cost functions that use Gibbs sample: cfn_b, cfn_b2, cfn_exp
        vector<size_t> state;
        if(stateP==NULL) {
            state = getGibbsState(fg,2*fg.doneIters());
        } else {
            state = *stateP;
        }
        assert(state.size()==fg.nrVars());

        vector<Prob> b1_adj;
        b1_adj.reserve(fg.nrVars());
        for(size_t i=0; i<state.size(); i++) {
            size_t n = fg.var(i).states();
            Prob delta(n,0.0);
            assert(/*0<=state[i] &&*/ state[i] < n);
            double b = fg.belief1(i)[state[i]];
            if(cfn_type==bbp_cfn_t::cfn_b) {
                delta[state[i]] = 1.0;
            } else if(cfn_type==bbp_cfn_t::cfn_b2) {
                delta[state[i]] = b;
            } else if(cfn_type==bbp_cfn_t::cfn_exp) {
                delta[state[i]] = exp(b);
            } else { abort(); }
            b1_adj.push_back(delta);
        }
        bbp.init(b1_adj);
    }
}

/// This function returns the actual value of the cost function whose
/// gradient with respect to singleton beliefs is given by
/// gibbsToB1Adj on the same arguments
Real gibbsCostFn(const BP_dual &fg, bbp_cfn_t cfn_type, const vector<size_t> *stateP) {
    double cf=0.0;
    if(cfn_type==bbp_cfn_t::cfn_bethe_ent) { // ignores state
        cf = -fg.logZ();
    } else if(cfn_type==bbp_cfn_t::cfn_var_ent) { // ignores state
        for(size_t i=0; i<fg.nrVars(); i++) {
            cf += -fg.belief1(i).entropy();
        }
    } else {
        assert(stateP != NULL);
        vector<size_t> state = *stateP;
        assert(state.size()==fg.nrVars());
        for(size_t i=0; i<fg.nrVars(); i++) {
            double b = fg.belief1(i)[state[i]];
            if(cfn_type==bbp_cfn_t::cfn_b) {
                cf += b;
            } else if(cfn_type==bbp_cfn_t::cfn_b2) {
                cf += b*b/2;
            } else if(cfn_type==bbp_cfn_t::cfn_exp) {
                cf += exp(b);
            } else { abort(); }
        }
    }
    return cf;
}

vector<Prob> get_zero_adj_2(const BP_dual& bp_dual) {
    vector<Prob> adj_2;
    adj_2.reserve(bp_dual.nrFactors());
    for(size_t I=0; I<bp_dual.nrFactors(); I++) {
        adj_2.push_back(Prob(bp_dual.factor(I).states(),Real(0.0)));
    }
    return adj_2;
}

vector<Prob> get_zero_adj_1(const BP_dual& bp_dual) {
    vector<Prob> adj_1;
    adj_1.reserve(bp_dual.nrVars());
    for(size_t i=0; i<bp_dual.nrVars(); i++) {
        adj_1.push_back(Prob(bp_dual.var(i).states(),Real(0.0)));
    }
    return adj_1;
}

#define foreach_iI(_fg)                                                 \
    for(size_t i=0, I=0, iI=0; (iI<_fg->nr_edges()) && ((i=_fg->edge(iI).first) || 1) && ((I=_fg->edge(iI).second) || 1); iI++)
    
#define foreach_iIj(_fg)                                        \
    for(size_t i=0, I=0, j=0, iIj=0; (iIj<_VFV_ind.size()) &&   \
            ((i=get<0>(_VFV_ind[iIj])) || 1) &&                 \
            ((I=get<1>(_VFV_ind[iIj])) || 1) &&                 \
            ((j=get<2>(_VFV_ind[iIj])) || 1); iIj++)

#define foreach_IiJ(_fg)                                        \
    for(size_t I=0, i=0, J=0, IiJ=0; (IiJ<_FVF_ind.size()) &&   \
            ((I=get<0>(_FVF_ind[IiJ])) || 1) &&                 \
            ((i=get<1>(_FVF_ind[IiJ])) || 1) &&                 \
            ((J=get<2>(_FVF_ind[IiJ])) || 1); IiJ++)


#define LOOP_ij(body) {                                 \
        size_t i_states = _fg->var(i).states();         \
        size_t j_states = _fg->var(j).states();         \
        if(_fg->var(i) > _fg->var(j)) {                 \
            size_t xij=0;                               \
            for(size_t xi=0; xi<i_states; xi++) {       \
                for(size_t xj=0; xj<j_states; xj++) {   \
                    body;                               \
                    xij++;                              \
                }                                       \
            }                                           \
        } else {                                        \
            size_t xij=0;                               \
            for(size_t xj=0; xj<j_states; xj++) {       \
                for(size_t xi=0; xi<i_states; xi++) {   \
                    body;                               \
                    xij++;                              \
                }                                       \
            }                                           \
        }                                               \
    }

void BBP::RegenerateInds() {
    // initialise _VFV and _VFV_ind
    {
        size_t ind=0;
        _VFV.resize(_fg->nrVars());
        _VFV_ind.clear();
        for(size_t i=0; i<_fg->nrVars(); i++) {
            foreach(size_t I, _fg->nbV(i)) {
                vector<size_t>& v = _VFV[i][I];
                v.resize(_fg->nrVars());
                foreach(size_t j, _fg->nbF(I)) {
                    if(j!=i) {
                        v[j] = ind++;
                        _VFV_ind.push_back(tuple<size_t,size_t,size_t>(i,I,j));
                    }
                }
            }
        }
    }
    // initialise _FVF
    {
        size_t ind=0;
        _FVF.resize(_fg->nrFactors());
        _FVF_ind.clear();
        for(size_t I=0; I<_fg->nrFactors(); I++) {
            foreach(size_t i, _fg->nbF(I)) {
                vector<size_t>& v = _FVF[I][i];
                v.resize(_fg->nrFactors());
                foreach(size_t J, _fg->nbV(i)) {
                    if(J!=I) {
                        v[J] = ind++;
                        _FVF_ind.push_back(tuple<size_t,size_t,size_t>(I,i,J));
                    }
                }
            }
        }
    }
}

void BBP::RegenerateT() {
    _T.clear();
    foreach_iI(_fg) {
        Prob prod(_fg->var(i).states(),1.0);
        foreach(size_t J, _fg->nbV(i)) {
            if(J!=I) {
                prod *= _bp_dual->msgM(J,i);
            }
        }
        _T.push_back(prod);
    }
}

void BBP::RegenerateU() {
    _U.clear();
    foreach_iI(_fg) {
        //     Prob prod(_fg->var(i).states(), 1.0);
        Prob prod(_fg->factor(I).states(), 1.0);
        foreach(size_t j, _fg->nbF(I)) {
            if(i!=j) {
                Prob n_jI(_bp_dual->msgN(j,I));
                const BP_dual::_ind_t* ind = &(_bp_dual->index(j,I));
                // multiply prod with n_jI
                for(size_t x_I = 0; x_I < prod.size(); x_I++)
                    prod[x_I] *= n_jI[(*ind)[x_I]];
                //         prod *= _bp_dual->msgN(j,I);
            }
        }
        _U.push_back(prod);
    }
}

void BBP::RegenerateS() {
    _S.clear();
    foreach_iIj(_fg) {
        Factor prod(_fg->factor(I));
        foreach(size_t k, _fg->nbF(I)) {
            if(k!=i && k!=j) {
                if(0) {
                    prod *= Factor(_fg->var(k), _bp_dual->msgN(k,I));
                } else {
                    const BP_dual::_ind_t& ind = _bp_dual->index(k,I);
                    Prob p(_bp_dual->msgN(k,I));
                    for(size_t xI=0; xI<prod.states(); xI++) {
                        prod.p()[xI] *= p[ind[xI]];
                    }
                }
            }
        }
        // XXX improve efficiency? 
        // "Marginalize" onto i|j (unnormalized)
        Prob marg;
        marg = prod.marginal(VarSet(_fg->var(i)) | VarSet(_fg->var(j)), false).p();
        _S.push_back(marg);
    }
}

void BBP::RegenerateR() {
    _R.clear();
    foreach_IiJ(_fg) {
        Prob prod(_fg->var(i).states(), 1.0);
        foreach(size_t K, _fg->nbV(i)) {
            if(K!=I && K!=J) {
                prod *= _bp_dual->msgM(K,i);
            }
        }
        _R.push_back(prod);
    }
}

void BBP::RegenerateInputs() {
    _adj_b_1_unnorm.clear();
    for(size_t i=0; i<_fg->nrVars(); i++) {
        _adj_b_1_unnorm.push_back(unnormAdjoint(_bp_dual->belief1(i).p(), _bp_dual->belief1Z(i), _adj_b_1[i]));
    }
    _adj_b_2_unnorm.clear();
    for(size_t I=0; I<_fg->nrFactors(); I++) {
        _adj_b_2_unnorm.push_back(unnormAdjoint(_bp_dual->belief2(I).p(), _bp_dual->belief2Z(I), _adj_b_2[I]));
    }
}

/// initialise members for factor adjoints (call after RegenerateInputs)
void BBP::RegeneratePsiAdjoints() {
    _adj_psi_1.clear();
    _adj_psi_1.reserve(_fg->nrVars());
    for(size_t i=0; i<_fg->nrVars(); i++) {
        Prob p(_adj_b_1_unnorm[i]);
        assert(p.size()==_fg->var(i).states());
        foreach(size_t I, _fg->nbV(i)) {
            p *= _bp_dual->msgM(I,i);
        }
        p += _init_adj_psi_1[i];
        _adj_psi_1.push_back(p);
    }
    _adj_psi_2.clear();
    _adj_psi_2.reserve(_fg->nrFactors());
    for(size_t I=0; I<_fg->nrFactors(); I++) {
        Prob p(_adj_b_2_unnorm[I]);
        assert(p.size()==_fg->factor(I).states());
        foreach(size_t i, _fg->nbF(I)) {
            Prob n_iI(_bp_dual->msgN(i,I));
            const BP_dual::_ind_t* ind = &(_bp_dual->index(i,I));
            // multiply prod with n_jI
            for(size_t x_I = 0; x_I < p.size(); x_I++)
                p[x_I] *= n_iI[(*ind)[x_I]];
        }
        p += _init_adj_psi_2[I];
        _adj_psi_2.push_back(p);
    }
}

/// initialise members for messages adjoints (call after RegenerateInputs)
void BBP::RegenerateMessageAdjoints() {
    _adj_n.resize(_fg->nr_edges());
    _adj_m.resize(_fg->nr_edges());
    _adj_n_unnorm.resize(_fg->nr_edges());
    _adj_m_unnorm.resize(_fg->nr_edges());
    _new_adj_n.clear();
    _new_adj_n.reserve(_fg->nr_edges());
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        size_t i = _fg->edge(iI).first;
        size_t I = _fg->edge(iI).second;
        Prob prod(_fg->factor(I).p());
        prod *= _adj_b_2_unnorm[I];
        foreach(size_t j, _fg->nbF(I)) {
            if(j!=i) { // for all j in I\i
                Prob n_jI(_bp_dual->msgN(j,I));; 
                const BP_dual::_ind_t* ind = &(_bp_dual->index(j,I));
                // multiply prod with n_jI
                for(size_t x_I = 0; x_I < prod.size(); x_I++)
                    prod[x_I] *= n_jI[(*ind)[x_I]];
            }
        }
        Prob marg( _fg->var(i).states(), 0.0 );
        // ind is the precalculated Index(i,I) i.e. to x_I == k corresponds x_i == ind[k]
        const BP_dual::_ind_t* ind = &(_bp_dual->index(i,I));
        for( size_t r = 0; r < prod.size(); r++ )
            marg[(*ind)[r]] += prod[r];

        _new_adj_n.push_back(marg);
        upMsgN(iI); // propagate new _new_adj_n to _adj_n and _adj_n_unnorm
    }
    _new_adj_m.clear();
    _new_adj_m.reserve(_fg->nr_edges());
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        size_t i = _fg->edge(iI).first;
        size_t I = _fg->edge(iI).second;
        //     Prob prod(_fg->var(i).states(),1.0);
        Prob prod(_adj_b_1_unnorm[i]);
        assert(prod.size()==_fg->var(i).states());
        foreach(size_t J, _fg->nbV(i)) {
            if(J!=I) { // for all J in i\I
                prod *= _bp_dual->msgM(J,i);
            }
        }
        _new_adj_m.push_back(prod);
        upMsgM(iI); // propagate new _new_adj_m to _adj_m and _adj_m_unnorm
    }
}

void BBP::Regenerate() {
    RegenerateInds();
    RegenerateT();
    RegenerateU();
    RegenerateS();
    RegenerateR();
    RegenerateInputs();
    RegeneratePsiAdjoints();
    RegenerateMessageAdjoints();
    _iters = 0;
}
  
void BBP::calcNewN(size_t iI) {
    size_t i=_fg->edge(iI).first;
    size_t I=_fg->edge(iI).second;
    _adj_psi_1[i] += T(i,I)*_adj_n_unnorm[iI];
    _trip_end_t vfv_iI = _VFV[i][I];
    Prob &new_adj_n_iI = _new_adj_n[iI];
    new_adj_n_iI = Prob(_fg->var(i).states(),0.0);
    foreach(size_t j, _fg->nbF(I)) {
        if(j!=i) {
            size_t iIj = vfv_iI[j];
            Prob &p = _S[iIj];
            size_t jI = _fg->VV2E(j,I);
            LOOP_ij(
                new_adj_n_iI[xi] += p[xij]*_adj_m_unnorm[jI][xj];
            );
        }
    }
}

void BBP::calcNewM(size_t iI) {
    size_t i=_fg->edge(iI).first;
    size_t I=_fg->edge(iI).second;

    Prob p(U(I,i));
    Prob adj(_adj_m_unnorm[iI]);
    const BP_dual::_ind_t* ind = &(_bp_dual->index(i,I));
    for(size_t x_I = 0; x_I < p.size(); x_I++)
        p[x_I] *= adj[(*ind)[x_I]];
    _adj_psi_2[I] += p;

    _new_adj_m[iI] = Prob(_fg->var(i).states(),0.0);
    _trip_end_t fvf_Ii = _FVF[I][i];
    foreach(size_t J, _fg->nbV(i)) {
        if(J!=I) {
            size_t iJ = _fg->VV2E(i,J);
            _new_adj_m[iI] += _R[fvf_Ii[J]]*_adj_n_unnorm[iJ];
        }
    }
}

void BBP::upMsgM(size_t iI) {
    _adj_m[iI] = _new_adj_m[iI];
    _adj_m_unnorm[iI] = unnormAdjoint(_bp_dual->msgM(iI), _bp_dual->zM(iI), _adj_m[iI]);
}

void BBP::upMsgN(size_t iI) {
    _adj_n[iI] = _new_adj_n[iI];
    _adj_n_unnorm[iI] = unnormAdjoint(_bp_dual->msgN(iI), _bp_dual->zN(iI), _adj_n[iI]);
}

void BBP::doParUpdate() {
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        calcNewM(iI);
        calcNewN(iI);
    }
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        upMsgM(iI);
        upMsgN(iI);
    }
}

Real BBP::getUnMsgMag() {
    Real s=0.0;
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        s += _adj_m_unnorm[iI].sumAbs();
        s += _adj_n_unnorm[iI].sumAbs();
    }
    return s/_fg->nr_edges();
}

void BBP::getMsgMags(Real &s, Real &new_s) {
    s=0.0; new_s=0.0;
    for(size_t iI=0; iI<_fg->nr_edges(); iI++) {
        s += _adj_m[iI].sumAbs();
        s += _adj_n[iI].sumAbs();
        new_s += _new_adj_m[iI].sumAbs();
        new_s += _new_adj_n[iI].sumAbs();
    }
    Real t = _fg->nr_edges();
    s /= t; new_s /= t;
}

/// run until change is less than given tolerance
void BBP::run(Real tol) {
    size_t minIters = 2;
    do {
        _iters++;
        doParUpdate();
    } while((_iters < minIters || getUnMsgMag() > tol) && _iters < maxIter());
    if(_iters==maxIter()) {
        Real s, new_s;
        getMsgMags(s,new_s);
        cerr << "Warning: BBP didn't converge in " << _iters
             << " iterations (unnorm message magnitude = " << getUnMsgMag()
             << ", norm message mags = " << s << " -> " << new_s
             << ")" << endl;
    }
}

vector<size_t> getGibbsState(const BP_dual& fg, size_t iters) {
    PropertySet gibbsProps;
    gibbsProps.Set("iters", iters);
    gibbsProps.Set("verbose", oneLess(fg.Verbose()));
    Gibbs gibbs(fg, gibbsProps);
    gibbs.run();
    return gibbs.state();
}

void testUnnormAdjoint(int n);

/// given a state for each variable, use numerical derivatives
/// (multiplying a factor containing a variable by psi_1 adjustments)
/// to verify accuracy of _adj_psi_1.
/// 'h' controls size of perturbation.
/// 'bbpTol' controls tolerance of BBP run.
double numericBBPTest(const BP_dual& bp_dual, bbp_cfn_t cfn, double bbpTol, double h) {
    //   testUnnormAdjoint(4);

    vector<size_t> state = getGibbsState(bp_dual,2*bp_dual.doneIters());
    // calculate the value of the unperturbed cost function
    Real cf0 = gibbsCostFn(bp_dual, cfn, &state);

    // run BBP to estimate adjoints
    BBP bbp(bp_dual);
    gibbsInitBBPCostFnAdj(bbp, bp_dual, cfn, &state);
    bbp.run(bbpTol);

    Real d=0;

    if(1) {
        // verify bbp.adj_psi_1

        // for each variable i
        for(size_t i=0; i<bp_dual.nrVars(); i++) {
            vector<double> adj_est;
            // for each value xi
            for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
                // create a copy of bp_dual which multiplies this into some factor containing the variable
                BP_dual bp_dual_prb(bp_dual);
                assert(bp_dual_prb.nbV(i).size()>0);

                // create a perturbed psi_1[i] (there is no single-variable
                // psi in libDAI so we initialize to all 1 and then multiply
                // into nearest multivariate factor)
                size_t n = bp_dual.var(i).states();
                Prob psi_1_prb(n, 1.0);
                psi_1_prb[xi] += h;
                psi_1_prb.normalize();

                // update bp_dual_prb
                size_t I = bp_dual_prb.nbV(i)[0]; // use the first factor in list of neighbours of i
                bp_dual_prb.factor(I) *= Factor(bp_dual.var(i), psi_1_prb);
                bp_dual_prb.init(bp_dual.var(i)); // reset messages involving i
                //       bp_dual_prb.init();
    
                // run the copy to convergence
                bp_dual_prb.run();
    
                // calculate the new value of the cost function
                Real cf_prb = gibbsCostFn(bp_dual_prb, cfn, &state);

                // use this to estimate the adjoint for i
                // XXX why is it off by a factor of 2?
                adj_est.push_back((cf_prb-cf0)/h);
            }
            Prob p_adj_est(adj_est.begin(), adj_est.end());
            // compare this numerical estimate to the BBP estimate; sum the distances
            cerr << "i: " << i
                 << ", p_adj_est: " << p_adj_est
                 << ", bbp.adj_psi_1(i): " << bbp.adj_psi_1(i) << endl;
            d += dist(p_adj_est, bbp.adj_psi_1(i), Prob::DISTL1);
        }
    }
    if(1) {
        // verify bbp.adj_n and bbp.adj_m

        // We actually want to check the responsiveness of objective
        // function to changes in the final messages. But at the end of a
        // BBP run, the message adjoints are for the initial messages (and
        // they should be close to zero, see paper). So this resets the
        // BBP adjoints to the refer to the desired final messages
        bbp.RegenerateMessageAdjoints();

        // for each variable i
        for(size_t i=0; i<bp_dual.nrVars(); i++) {
            // for each factor I ~ i
            foreach(size_t I, bp_dual.nbV(i)) {
                vector<double> adj_n_est;
                // for each value xi
                for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
                    BP_dual bp_dual_prb(bp_dual);
                    // make h-sized change to newMsgN
                    bp_dual_prb.newMsgN(i,I)[xi] += h;
                    // recalculate beliefs
                    bp_dual_prb.CalcBeliefs();
                    // get cost function value
                    Real cf_prb = gibbsCostFn(bp_dual_prb, cfn, &state);
                    // add it to list of adjoints
                    adj_n_est.push_back((cf_prb-cf0)/h);
                }
        
                vector<double> adj_m_est;
                // for each value xi
                for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
                    BP_dual bp_dual_prb(bp_dual);
                    // make h-sized change to newMsgM
                    bp_dual_prb.newMsgM(I,i)[xi] += h;
                    // recalculate beliefs
                    bp_dual_prb.CalcBeliefs();
                    // get cost function value
                    Real cf_prb = gibbsCostFn(bp_dual_prb, cfn, &state);
                    // add it to list of adjoints
                    adj_m_est.push_back((cf_prb-cf0)/h);
                }

                Prob p_adj_n_est(adj_n_est.begin(), adj_n_est.end());
                // compare this numerical estimate to the BBP estimate; sum the distances
                cerr << "i: " << i << ", I: " << I
                     << ", adj_n_est: " << p_adj_n_est
                     << ", bbp.adj_n(i,I): " << bbp.adj_n(i,I) << endl;
                d += dist(p_adj_n_est, bbp.adj_n(i,I), Prob::DISTL1);

                Prob p_adj_m_est(adj_m_est.begin(), adj_m_est.end());
                // compare this numerical estimate to the BBP estimate; sum the distances
                cerr << "i: " << i << ", I: " << I
                     << ", adj_m_est: " << p_adj_m_est
                     << ", bbp.adj_m(I,i): " << bbp.adj_m(I,i) << endl;
                d += dist(p_adj_m_est, bbp.adj_m(I,i), Prob::DISTL1);
            }
        }
    }
    if(0) {
        // verify bbp.adj_b_1
        for(size_t i=0; i<bp_dual.nrVars(); i++) {
            vector<double> adj_b_1_est;
            // for each value xi
            for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
                BP_dual bp_dual_prb(bp_dual);

                // make h-sized change to b_1(i)[x_i]
                bp_dual_prb._beliefs.b1[i][xi] += h;

                // get cost function value
                Real cf_prb = gibbsCostFn(bp_dual_prb, cfn, &state);

                // add it to list of adjoints
                adj_b_1_est.push_back((cf_prb-cf0)/h);
            }
            Prob p_adj_b_1_est(adj_b_1_est.begin(), adj_b_1_est.end());
            // compare this numerical estimate to the BBP estimate; sum the distances
            cerr << "i: " << i
                 << ", adj_b_1_est: " << p_adj_b_1_est
                 << ", bbp.adj_b_1(i): " << bbp.adj_b_1(i) << endl;
            d += dist(p_adj_b_1_est, bbp.adj_b_1(i), Prob::DISTL1);
        }
    }

    // return total of distances
    return d;
}

void testUnnormAdjoint(int n) {
    double h = 1.0e-5;
    Prob q_unnorm(n);
    // create a random q_unnorm
    for(int i=0; i<n; i++) {
        q_unnorm[i] = exp(4*rnd_stdnormal());
    }
    // normalize it to get q and Z_q
    Prob q(q_unnorm);
    double Z_q = q.normalize();
    // cost function V just selects one (random) element of q
    int vk = rnd_multi(n);
    double V = q[vk];
    // calculate adj_q for this V
    Prob adj_q(n, 0.0);
    adj_q[vk] = 1;
    // use unnormAdjoint to get adj_q_unnorm
    Prob adj_q_unnorm = unnormAdjoint(q, Z_q, adj_q);
    // with perturbations of q_unnorm, test that adj_q_unnorm is correct
    Prob adj_q_unnorm_numeric(n, 0.0);
    for(int i=0; i<n; i++) {
        Prob q_unnorm_prb = q_unnorm;
        q_unnorm_prb[i] += h;
        Prob q_prb(q_unnorm_prb);
        q_prb.normalize();
        DAI_PV(q_unnorm_prb);
        DAI_PV(q_prb);
        double V_prb = q_prb[vk];
        adj_q_unnorm_numeric[i] = (V_prb-V)/h;
    }
    DAI_PV(q_unnorm);
    DAI_PV(q);
    DAI_PV(Z_q);
    DAI_PV(vk);
    DAI_PV(V);
    DAI_PV(adj_q_unnorm);
    DAI_PV(adj_q_unnorm_numeric);
}

void makeBBPGraph(const BP_dual& bp_dual, bbp_cfn_t cfn) {
    double tiny=1.0e-7;
    vector<size_t> state = getGibbsState(bp_dual,2*bp_dual.doneIters());
    const char *graphs_dir = "./bbp-data";
    mkdir(graphs_dir, -1);
    size_t num_samples=500;
    for(size_t i=0; i<bp_dual.nrVars(); i++) {
        for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
            if(bp_dual.belief1(i)[xi]<tiny || abs(bp_dual.belief1(i)[xi]-1)<tiny)
                continue;
            char *fn;
            asprintf(&fn,"%s/i:%d-xi:%d.out",graphs_dir,(int)i,(int)xi);
            ofstream os(fn, ios_base::out|ios_base::trunc);
      
            for(size_t n=0; n<num_samples; n++) {
                double c=((double)n)/(num_samples-1);
                Prob psi_1_prb(bp_dual.var(i).states(),1.0-c);
                psi_1_prb[xi] = 1.0;

                BP_dual bp_dual_prb(bp_dual);
                assert(bp_dual_prb.nbV(i).size()>0);
                size_t I = bp_dual_prb.nbV(i)[0]; // use the first factor in list of neighbours of i
                bp_dual_prb.factor(I) *= Factor(bp_dual.var(i), psi_1_prb);
                bp_dual_prb.init(bp_dual.var(i)); // reset messages involving i
    
                // run the copy to convergence
                bp_dual_prb.run();
        
                Real cf = gibbsCostFn(bp_dual_prb, cfn, &state);

                os << n << "\t" << cf << endl;
            }
      
            os.close();
            free(fn);
        }
    }
}

} // end of namespace dai
