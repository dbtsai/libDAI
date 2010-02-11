/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2007       Bastian Wemmenhove
 *  Copyright (C) 2007-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2007       Radboud University Nijmegen, The Netherlands
 */


#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <dai/mr.h>
#include <dai/bp.h>
#include <dai/jtree.h>
#include <dai/util.h>


namespace dai {


using namespace std;


const char *MR::Name = "MR";


void MR::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );
    DAI_ASSERT( opts.hasKey("verbose") );
    DAI_ASSERT( opts.hasKey("updates") );
    DAI_ASSERT( opts.hasKey("inits") );

    props.tol = opts.getStringAs<Real>("tol");
    props.verbose = opts.getStringAs<size_t>("verbose");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");
    props.inits = opts.getStringAs<Properties::InitType>("inits");
}


PropertySet MR::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "verbose", props.verbose );
    opts.Set( "updates", props.updates );
    opts.Set( "inits", props.inits );
    return opts;
}


string MR::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "verbose=" << props.verbose << ",";
    s << "updates=" << props.updates << ",";
    s << "inits=" << props.inits << "]";
    return s.str();
}


void MR::init(size_t Nin, Real *_w, Real *_th) {
    size_t i,j;

    N = Nin;

    con.resize(N);
    nb.resize(N);
    tJ.resize(N);
    for(i=0; i<N; i++ ) {
        nb[i].resize(kmax);
        tJ[i].resize(kmax);
        con[i]=0;
        for(j=0; j<N; j++ )
            if( _w[i*N+j] != 0.0 ) {
                nb[i][con[i]] = j;
                tJ[i][con[i]] = tanh(_w[i*N+j]);
                con[i]++;
            }
    }

    theta.resize(N);
    for(i=0; i<N; i++)
      theta[i] = _th[i];
}


Real MR::init_cor_resp() {
    size_t j,k,l, runx,i2;
    Real variab1, variab2;
    Real md, maxdev;
    Real thbJsite[kmax];
    Real xinter;
    Real rinter;
    Real res[kmax];
    size_t s2;
    size_t flag;
    size_t concav;
    size_t runs = 3000;
    Real eps = 0.2;
    size_t cavity;

    vector<vector<Real> > tJ_org;
    vector<vector<size_t> > nb_org;
    vector<size_t> con_org;
    vector<Real> theta_org;

    vector<Real> xfield(N*kmax,0.0);
    vector<Real> rfield(N*kmax,0.0);
    vector<Real> Hfield(N,0.0);
    vector<Real> devs(N*kmax,0.0);
    vector<Real> devs2(N*kmax,0.0);
    vector<Real> dev(N,0.0);
    vector<Real> avmag(N,0.0);

    // save original tJ, nb
    nb_org = nb;
    tJ_org = tJ;
    con_org = con;
    theta_org = theta;

    maxdev = 0.0;
    for(cavity=0; cavity<N; cavity++){    // for each spin to be removed
        con = con_org;
        concav=con[cavity];

        nb = nb_org;
        tJ = tJ_org;

        //  Adapt the graph variables nb[], tJ[] and con[]
        for(size_t i=0; i<con[cavity]; i++) {
            size_t ij = nb[cavity][i];
            flag=0;
            j=0;
            do{
                if(nb[ij][j]==cavity){
                    while(j<(con[ij]-1)){
                        nb[ij][j]=nb[ij][j+1];
                        tJ[ij][j] = tJ[ij][j+1];
                        j++;
                    }
                flag=1;
                }
                j++;
            } while(flag==0);
        }
        for(size_t i=0; i<con[cavity]; i++)
            con[nb[cavity][i]]--;
        con[cavity] = 0;


        // Do everything starting from the new graph********

        makekindex();
        theta = theta_org;

        for(size_t i=0; i<kmax*N; i++)
            xfield[i] = 3.0*(2*rnd_uniform()-1.);

        for(i2=0; i2<concav; i2++){ // Subsequently apply a field to each cavity spin ****************

            s2 = nb[cavity][i2];    // identify the index of the cavity spin
            for(size_t i=0; i<con[s2]; i++)
                rfield[kmax*s2+i] = 1.;

            runx=0;
            do {      // From here start the response and belief propagation
                runx++;
                md=0.0;
                for(k=0; k<N; k++){
                    if(k!=cavity) {
                        for(size_t i=0; i<con[k]; i++)
                            thbJsite[i] = tJ[k][i];
                        for(l=0; l<con[k]; l++){
                            xinter = 1.;
                            rinter = 0.;
                            if(k==s2) rinter += 1.;
                            for(j=0; j<con[k]; j++)
                                if(j!=l){
                                    variab2 = tanh(xfield[kmax*nb[k][j]+kindex[k][j]]);
                                    variab1 = thbJsite[j]*variab2;
                                    xinter *= (1.+variab1)/(1.-variab1);

                                    rinter += thbJsite[j]*rfield[kmax*nb[k][j]+kindex[k][j]]*(1-variab2*variab2)/(1-variab1*variab1);
                                }

                            variab1 = 0.5*log(xinter);
                            xinter = variab1 + theta[k];
                            devs[kmax*k+l] = xinter-xfield[kmax*k+l];
                            xfield[kmax*k+l] = xfield[kmax*k+l]+devs[kmax*k+l]*eps;
                            if( fabs(devs[kmax*k+l]) > md )
                                md = fabs(devs[kmax*k+l]);

                            devs2[kmax*k+l] = rinter-rfield[kmax*k+l];
                            rfield[kmax*k+l] = rfield[kmax*k+l]+devs2[kmax*k+l]*eps;
                            if( fabs(devs2[kmax*k+l]) > md )
                                md = fabs(devs2[kmax*k+l]);
                        }
                    }
                }
            } while((md > props.tol)&&(runx<runs)); // Precision condition reached -> BP and RP finished
            if(runx==runs)
                if( props.verbose >= 2 )
                    cerr << "init_cor_resp: Convergence not reached (md=" << md << ")..." << endl;
            if(md > maxdev)
                maxdev = md;

            // compute the observables (i.e. magnetizations and responses)******

            for(size_t i=0; i<concav; i++){
                rinter = 0.;
                xinter = 1.;
                if(i!=i2)
                    for(j=0; j<con[nb[cavity][i]]; j++){
                        variab2 = tanh(xfield[kmax*nb[nb[cavity][i]][j]+kindex[nb[cavity][i]][j]]);
                        variab1 = tJ[nb[cavity][i]][j]*variab2;
                        rinter +=  tJ[nb[cavity][i]][j]*rfield[kmax*nb[nb[cavity][i]][j]+kindex[nb[cavity][i]][j]]*(1-variab2*variab2)/(1-variab1*variab1);
                        xinter *= (1.+variab1)/(1.-variab1);
                    }
                xinter = tanh(0.5*log(xinter)+theta[nb[cavity][i]]);
                res[i] = rinter*(1-xinter*xinter);
            }

            // *******************

            for(size_t i=0; i<concav; i++)
                if(nb[cavity][i]!=s2)
            //      if(i!=i2)
                    cors[cavity][i2][i] = res[i];
                else
                    cors[cavity][i2][i] = 0;
        } // close for i2 = 0...concav
    }

    // restore nb, tJ, con
    tJ = tJ_org;
    nb = nb_org;
    con = con_org;
    theta = theta_org;

    return maxdev;
}


Real MR::T(size_t i, sub_nb A) {
    sub_nb _nbi_min_A(con[i]);
    _nbi_min_A.set();
    _nbi_min_A &= ~A;

    Real res = theta[i];
    for( size_t _j = 0; _j < _nbi_min_A.size(); _j++ )
        if( _nbi_min_A.test(_j) )
            res += atanh(tJ[i][_j] * M[i][_j]);
    return tanh(res);
}


Real MR::T(size_t i, size_t _j) {
    sub_nb j(con[i]);
    j.set(_j);
    return T(i,j);
}


Real MR::Omega(size_t i, size_t _j, size_t _l) {
    sub_nb jl(con[i]);
    jl.set(_j);
    jl.set(_l);
    Real Tijl = T(i,jl);
    return Tijl / (1.0 + tJ[i][_l] * M[i][_l] * Tijl);
}


Real MR::Gamma(size_t i, size_t _j, size_t _l1, size_t _l2) {
    sub_nb jll(con[i]);
    jll.set(_j);
    Real Tij = T(i,jll);
    jll.set(_l1);
    jll.set(_l2);
    Real Tijll = T(i,jll);

    return (Tijll - Tij) / (1.0 + tJ[i][_l1] * tJ[i][_l2] * M[i][_l1] * M[i][_l2] + tJ[i][_l1] * M[i][_l1] * Tijll + tJ[i][_l2] * M[i][_l2] * Tijll);
}


Real MR::Gamma(size_t i, size_t _l1, size_t _l2) {
    sub_nb ll(con[i]);
    Real Ti = T(i,ll);
    ll.set(_l1);
    ll.set(_l2);
    Real Till = T(i,ll);

    return (Till - Ti) / (1.0 + tJ[i][_l1] * tJ[i][_l2] * M[i][_l1] * M[i][_l2] + tJ[i][_l1] * M[i][_l1] * Till + tJ[i][_l2] * M[i][_l2] * Till);
}


Real MR::_tJ(size_t i, sub_nb A) {
    sub_nb::size_type _j = A.find_first();
    if( _j == sub_nb::npos )
        return 1.0;
    else
        return tJ[i][_j] * _tJ(i, A.reset(_j));
}


Real MR::appM(size_t i, sub_nb A) {
    sub_nb::size_type _j = A.find_first();
    if( _j == sub_nb::npos )
        return 1.0;
    else {
        sub_nb A_j(A); A_j.reset(_j);

        Real result = M[i][_j] * appM(i, A_j);
        for( size_t _k = 0; _k < A_j.size(); _k++ )
            if( A_j.test(_k) ) {
                sub_nb A_jk(A_j); A_jk.reset(_k);
                result += cors[i][_j][_k] * appM(i,A_jk);
            }

        return result;
    }
}


void MR::sum_subs(size_t j, sub_nb A, Real *sum_even, Real *sum_odd) {
    *sum_even = 0.0;
    *sum_odd = 0.0;

    sub_nb B(A.size());
    do {
        if( B.count() % 2 )
            *sum_odd += _tJ(j,B) * appM(j,B);
        else
            *sum_even += _tJ(j,B) * appM(j,B);

        // calc next subset B
        size_t bit = 0;
        for( ; bit < A.size(); bit++ )
            if( A.test(bit) ) {
                if( B.test(bit) )
                    B.reset(bit);
                else {
                    B.set(bit);
                    break;
                }
            }
    } while (!B.none());
}


void MR::solvemcav() {
    Real sum_even, sum_odd;
    Real maxdev;
    size_t maxruns = 1000;

    makekindex();
    for(size_t i=0; i<N; i++)
        for(size_t _j=0; _j<con[i]; _j++)
            M[i][_j]=0.1;

    size_t run=0;
    do {
        maxdev=0.0;
        run++;
        for(size_t i=0; i<N; i++){ // for all i
            for(size_t _j=0; _j<con[i]; _j++){ // for all j in N_i
                size_t _i = kindex[i][_j];
                size_t j = nb[i][_j];
                DAI_ASSERT( nb[j][_i] == i );

                Real newM = 0.0;
                if( props.updates == Properties::UpdateType::FULL ) {
                    // find indices in nb[j] that do not correspond with i
                    sub_nb _nbj_min_i(con[j]);
                    _nbj_min_i.set();
                    _nbj_min_i.reset(kindex[i][_j]);

                    // find indices in nb[i] that do not correspond with j
                    sub_nb _nbi_min_j(con[i]);
                    _nbi_min_j.set();
                    _nbi_min_j.reset(_j);

                    sum_subs(j, _nbj_min_i, &sum_even, &sum_odd);
                    newM = (tanh(theta[j]) * sum_even + sum_odd) / (sum_even + tanh(theta[j]) * sum_odd);

                    sum_subs(i, _nbi_min_j, &sum_even, &sum_odd);
                    Real denom = sum_even + tanh(theta[i]) * sum_odd;
                    Real numer = 0.0;
                    for(size_t _k=0; _k<con[i]; _k++) if(_k != _j) {
                        sub_nb _nbi_min_jk(_nbi_min_j);
                        _nbi_min_jk.reset(_k);
                        sum_subs(i, _nbi_min_jk, &sum_even, &sum_odd);
                        numer += tJ[i][_k] * cors[i][_j][_k] * (tanh(theta[i]) * sum_even + sum_odd);
                    }
                    newM -= numer / denom;
                } else if( props.updates == Properties::UpdateType::LINEAR ) {
                    newM = T(j,_i);
                    for(size_t _l=0; _l<con[i]; _l++) if( _l != _j )
                        newM -= Omega(i,_j,_l) * tJ[i][_l] * cors[i][_j][_l];
                    for(size_t _l1=0; _l1<con[j]; _l1++) if( _l1 != _i )
                        for( size_t _l2=_l1+1; _l2<con[j]; _l2++) if( _l2 != _i)
                            newM += Gamma(j,_i,_l1,_l2) * tJ[j][_l1] * tJ[j][_l2] * cors[j][_l1][_l2];
                }

                Real dev = newM - M[i][_j];
//              dev *= 0.02;
                if( fabs(dev) >= maxdev )
                    maxdev = fabs(dev);

                newM = M[i][_j] + dev;
                if( fabs(newM) > 1.0 )
                    newM = sign(newM);
                M[i][_j] = newM;
            }
        }
    } while((maxdev>props.tol)&&(run<maxruns));

    _iters = run;
    if( maxdev > _maxdiff )
        _maxdiff = maxdev;

    if(run==maxruns){
        if( props.verbose >= 1 )
            cerr << "solve_mcav: Convergence not reached (maxdev=" << maxdev << ")..." << endl;
    }
}


void MR::solveM() {
    for(size_t i=0; i<N; i++) {
        if( props.updates == Properties::UpdateType::FULL ) {
            // find indices in nb[i]
            sub_nb _nbi(con[i]);
            _nbi.set();

            // calc numerator1 and denominator1
            Real sum_even, sum_odd;
            sum_subs(i, _nbi, &sum_even, &sum_odd);

            Mag[i] = (tanh(theta[i]) * sum_even + sum_odd) / (sum_even + tanh(theta[i]) * sum_odd);

        } else if( props.updates == Properties::UpdateType::LINEAR ) {
            sub_nb empty(con[i]);
            Mag[i] = T(i,empty);

            for(size_t _l1=0; _l1<con[i]; _l1++)
                for( size_t _l2=_l1+1; _l2<con[i]; _l2++)
                    Mag[i] += Gamma(i,_l1,_l2) * tJ[i][_l1] * tJ[i][_l2] * cors[i][_l1][_l2];
        }
        if(fabs(Mag[i])>1.)
            Mag[i] = sign(Mag[i]);
    }
}


Real MR::init_cor() {
    Real md = 0.0;
    for( size_t i = 0; i < nrVars(); i++ ) {
        vector<Factor> pairq;
        if( props.inits == Properties::InitType::CLAMPING ) {
            BP bpcav(*this, PropertySet()("updates", string("SEQMAX"))("tol", (Real)1.0e-9)("maxiter", (size_t)10000)("verbose", (size_t)0)("logdomain", false));
            bpcav.makeCavity( i );
            pairq = calcPairBeliefs( bpcav, delta(i), false, true );
            md = std::max( md, bpcav.maxDiff() );
        } else if( props.inits == Properties::InitType::EXACT ) {
            JTree jtcav(*this, PropertySet()("updates", string("HUGIN"))("verbose", (size_t)0) );
            jtcav.makeCavity( i );
            pairq = calcPairBeliefs( jtcav, delta(i), false, true );
        }
        for( size_t jk = 0; jk < pairq.size(); jk++ ) {
            VarSet::const_iterator kit = pairq[jk].vars().begin();
            size_t j = findVar( *(kit) );
            size_t k = findVar( *(++kit) );
            pairq[jk].normalize();
            Real cor = (pairq[jk][3] - pairq[jk][2] - pairq[jk][1] + pairq[jk][0]) - (pairq[jk][3] + pairq[jk][2] - pairq[jk][1] - pairq[jk][0]) * (pairq[jk][3] - pairq[jk][2] + pairq[jk][1] - pairq[jk][0]);
            for( size_t _j = 0; _j < con[i]; _j++ ) if( nb[i][_j] == j )
                for( size_t _k = 0; _k < con[i]; _k++ ) if( nb[i][_k] == k ) {
                    cors[i][_j][_k] = cor;
                    cors[i][_k][_j] = cor;
                }
        }
    }
    return md;
}


string MR::identify() const {
    return string(Name) + printProperties();
}


Real MR::run() {
    if( supported ) {
        if( props.verbose >= 1 )
            cerr << "Starting " << identify() << "...";

        double tic = toc();

        M.resize(N);
        for(size_t i=0; i<N; i++)
          M[i].resize(kmax);

        cors.resize(N);
        for(size_t i=0; i<N; i++)
          cors[i].resize(kmax);
        for(size_t i=0; i<N; i++)
          for(size_t j=0; j<kmax; j++)
            cors[i][j].resize(kmax);

        kindex.resize(N);
        for(size_t i=0; i<N; i++)
          kindex[i].resize(kmax);

        if( props.inits == Properties::InitType::RESPPROP ) {
            Real md = init_cor_resp();
            if( md > _maxdiff )
                _maxdiff = md;
        } else if( props.inits == Properties::InitType::EXACT ) {
            Real md = init_cor();
            if( md > _maxdiff )
                _maxdiff = md;
        }
        else if( props.inits == Properties::InitType::CLAMPING ) {
            Real md = init_cor();
            if( md > _maxdiff )
                _maxdiff = md;
        }

        solvemcav();

        Mag.resize(N);
        solveM();

        if( props.verbose >= 1 )
            cerr << Name << " needed " << toc() - tic << " seconds." << endl;

        return _maxdiff;
    } else
        return 1.0;
}


void MR::makekindex() {
    for(size_t i=0; i<N; i++)
        for(size_t j=0; j<con[i]; j++) {
            size_t ij = nb[i][j];       // ij is the j'th neighbour of spin i
            size_t k=0;
            while( nb[ij][k] != i )
                k++;
            kindex[i][j] = k;   // the j'th neighbour of spin i has spin i as its k'th neighbour
        }
}


Factor MR::beliefV( size_t i ) const {
    if( supported ) {
        Prob x(2);
        x[0] = 0.5 - Mag[i] / 2.0;
        x[1] = 0.5 + Mag[i] / 2.0;

        return Factor( var(i), x );
    } else
        return Factor();
}


vector<Factor> MR::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); i++ )
        result.push_back( belief( var(i) ) );
    return result;
}



MR::MR( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), supported(true), _maxdiff(0.0), _iters(0) {
    setProperties( opts );

    // check whether all vars in fg are binary
    // check whether connectivity is <= kmax
    for( size_t i = 0; i < fg.nrVars(); i++ )
        if( (fg.var(i).states() > 2) || (fg.delta(i).size() > kmax) ) {
            supported = false;
            break;
        }

    if( !supported )
        DAI_THROWE(NOT_IMPLEMENTED,"MR only supports binary variables with low connectivity");

    // check whether all interactions are pairwise or single
    for( size_t I = 0; I < fg.nrFactors(); I++ )
        if( fg.factor(I).vars().size() > 2 ) {
            supported = false;
            break;
        }

    if( !supported )
        DAI_THROWE(NOT_IMPLEMENTED,"MR does not support higher order interactions (only single and pairwise are supported)");

    // create w and th
    size_t Nin = fg.nrVars();

    Real *w = new Real[Nin*Nin];
    Real *th = new Real[Nin];

    for( size_t i = 0; i < Nin; i++ ) {
        th[i] = 0.0;
        for( size_t j = 0; j < Nin; j++ )
            w[i*Nin+j] = 0.0;
    }

    for( size_t I = 0; I < fg.nrFactors(); I++ ) {
        const Factor &psi = fg.factor(I);
        if( psi.vars().size() == 1 ) {
            size_t i = fg.findVar( *(psi.vars().begin()) );
            th[i] += 0.5 * log(psi[1] / psi[0]);
        } else if( psi.vars().size() == 2 ) {
            size_t i = fg.findVar( *(psi.vars().begin()) );
            VarSet::const_iterator jit = psi.vars().begin();
            size_t j = fg.findVar( *(++jit) );

            w[i*Nin+j] += 0.25 * log(psi[3] * psi[0] / (psi[2] * psi[1]));
            w[j*Nin+i] += 0.25 * log(psi[3] * psi[0] / (psi[2] * psi[1]));

            th[i] += 0.25 * log(psi[3] / psi[2] * psi[1] / psi[0]);
            th[j] += 0.25 * log(psi[3] / psi[1] * psi[2] / psi[0]);
        }
    }

    init(Nin, w, th);

    delete th;
    delete w;
}


} // end of namespace dai
