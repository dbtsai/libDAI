/*  Copyright (C) 2009  Frederik Eaton [frederik at ofb dot net]

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


#include <dai/bp.h>
#include <dai/bbp.h>
#include <dai/gibbs.h>
#include <dai/util.h>
#include <dai/bipgraph.h>


namespace dai {


using namespace std;


typedef BipartiteGraph::Neighbor Neighbor;


Prob unnormAdjoint( const Prob &w, Real Z_w, const Prob &adj_w ) {
    assert( w.size() == adj_w.size() );
    Prob adj_w_unnorm( w.size(), 0.0 );
    Real s = 0.0;
    for( size_t i = 0; i < w.size(); i++ )
        s += w[i] * adj_w[i];
    for( size_t i = 0; i < w.size(); i++ )
        adj_w_unnorm[i] = (adj_w[i] - s) / Z_w;
    return adj_w_unnorm;
//  THIS WOULD BE ABOUT 50% SLOWER:  return (adj_w - (w * adj_w).sum()) / Z_w;
}


std::vector<size_t> getGibbsState( const InfAlg &ia, size_t iters ) {
    PropertySet gibbsProps;
    gibbsProps.Set("iters", iters);
    gibbsProps.Set("verbose", size_t(0));
    Gibbs gibbs( ia.fg(), gibbsProps );
    gibbs.run();
    return gibbs.state();
}


size_t getFactorEntryForState( const FactorGraph &fg, size_t I, const vector<size_t> &state ) {
    size_t f_entry = 0;
    for( int _j = fg.nbF(I).size() - 1; _j >= 0; _j-- ) {
        // note that iterating over nbF(I) yields the same ordering
        // of variables as iterating over factor(I).vars()
        size_t j = fg.nbF(I)[_j];
        f_entry *= fg.var(j).states();
        f_entry += state[j];
    }
    return f_entry;
}


#define LOOP_ij(body) {                             \
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
    // initialise _indices
    //     typedef std::vector<size_t>        _ind_t;
    //     std::vector<std::vector<_ind_t> >  _indices; 
    _indices.resize( _fg->nrVars() );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        _indices[i].clear();
        _indices[i].reserve( _fg->nbV(i).size() );
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            _ind_t index;
            index.reserve( _fg->factor(I).states() );
            for( IndexFor k(_fg->var(i), _fg->factor(I).vars()); k >= 0; ++k )
                index.push_back( k );
            _indices[i].push_back( index );
        }
    }
}


void BBP::RegenerateT() {
    // _T[i][_I]
    _T.resize( _fg->nrVars() );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        _T[i].resize( _fg->nbV(i).size() );
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            Prob prod( _fg->var(i).states(), 1.0 );
            foreach( const Neighbor &J, _fg->nbV(i) )
                if( J.node != I.node )
                    prod *= _bp_dual.msgM( i, J.iter );
            _T[i][I.iter] = prod;
        }
    }
}


void BBP::RegenerateU() {
    // _U[I][_i]
    _U.resize( _fg->nrFactors() );
    for( size_t I = 0; I < _fg->nrFactors(); I++ ) {
        _U[I].resize( _fg->nbF(I).size() );
        foreach( const Neighbor &i, _fg->nbF(I) ) {
            Prob prod( _fg->factor(I).states(), 1.0 );
            foreach( const Neighbor &j, _fg->nbF(I) )
                if( i.node != j.node ) {
                    Prob n_jI( _bp_dual.msgN( j, j.dual ) );
                    const _ind_t &ind = _index( j, j.dual );
                    // multiply prod by n_jI
                    for( size_t x_I = 0; x_I < prod.size(); x_I++ )
                        prod[x_I] *= n_jI[ind[x_I]];
                }
            _U[I][i.iter] = prod;
        }
    }
}


void BBP::RegenerateS() {
    // _S[i][_I][_j]
    _S.resize( _fg->nrVars() );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        _S[i].resize( _fg->nbV(i).size() );
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            _S[i][I.iter].resize( _fg->nbF(I).size() );
            foreach( const Neighbor &j, _fg->nbF(I) )
                if( i != j ) {
                    Factor prod( _fg->factor(I) );
                    foreach( const Neighbor &k, _fg->nbF(I) ) {
                        if( k != i && k.node != j.node ) {
                            const _ind_t &ind = _index( k, k.dual );
                            Prob p( _bp_dual.msgN( k, k.dual ) );
                            for( size_t x_I = 0; x_I < prod.states(); x_I++ )
                                prod.p()[x_I] *= p[ind[x_I]];
                        }
                    }
                    // "Marginalize" onto i|j (unnormalized)
                    Prob marg;
                    marg = prod.marginal( VarSet(_fg->var(i), _fg->var(j)), false ).p();
                    _S[i][I.iter][j.iter] = marg;
                }
        }
    }
}


void BBP::RegenerateR() {
    // _R[I][_i][_J]
    _R.resize( _fg->nrFactors() );
    for( size_t I = 0; I < _fg->nrFactors(); I++ ) {
        _R[I].resize( _fg->nbF(I).size() );
        foreach( const Neighbor &i, _fg->nbF(I) ) {
            _R[I][i.iter].resize( _fg->nbV(i).size() );
            foreach( const Neighbor &J, _fg->nbV(i) ) {
                if( I != J ) {
                    Prob prod( _fg->var(i).states(), 1.0 );
                    foreach( const Neighbor &K, _fg->nbV(i) )
                        if( K.node != I && K.node != J.node ) 
                            prod *= _bp_dual.msgM( i, K.iter );
                    _R[I][i.iter][J.iter] = prod;
                }
            }
        }
    }
}


void BBP::RegenerateInputs() {
    _adj_b_V_unnorm.clear();
    _adj_b_V_unnorm.reserve( _fg->nrVars() );
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        _adj_b_V_unnorm.push_back( unnormAdjoint( _bp_dual.beliefV(i).p(), _bp_dual.beliefVZ(i), _adj_b_V[i] ) );
    _adj_b_F_unnorm.clear();
    _adj_b_F_unnorm.reserve( _fg->nrFactors() );
    for( size_t I = 0; I < _fg->nrFactors(); I++ )
        _adj_b_F_unnorm.push_back( unnormAdjoint( _bp_dual.beliefF(I).p(), _bp_dual.beliefFZ(I), _adj_b_F[I] ) );
}


void BBP::RegeneratePsiAdjoints() {
    _adj_psi_V.clear();
    _adj_psi_V.reserve( _fg->nrVars() );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        Prob p( _adj_b_V_unnorm[i] );
        assert( p.size() == _fg->var(i).states() );
        foreach( const Neighbor &I, _fg->nbV(i) )
            p *= _bp_dual.msgM( i, I.iter );
        p += _init_adj_psi_V[i];
        _adj_psi_V.push_back( p );
    }
    _adj_psi_F.clear();
    _adj_psi_F.reserve( _fg->nrFactors() );
    for( size_t I = 0; I < _fg->nrFactors(); I++ ) {
        Prob p( _adj_b_F_unnorm[I] );
        assert( p.size() == _fg->factor(I).states() );
        foreach( const Neighbor &i, _fg->nbF(I) ) {
            Prob n_iI( _bp_dual.msgN( i, i.dual ) );
            const _ind_t& ind = _index( i, i.dual );
            // multiply prod with n_jI
            for( size_t x_I = 0; x_I < p.size(); x_I++ )
                p[x_I] *= n_iI[ind[x_I]];
        }
        p += _init_adj_psi_F[I];
        _adj_psi_F.push_back( p );
    }
}


void BBP::RegenerateParMessageAdjoints() {
    size_t nv = _fg->nrVars();
    _adj_n.resize( nv );
    _adj_m.resize( nv );
    _adj_n_unnorm.resize( nv );
    _adj_m_unnorm.resize( nv );
    _new_adj_n.resize( nv );
    _new_adj_m.resize( nv );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        size_t n_i = _fg->nbV(i).size();
        _adj_n[i].resize( n_i );
        _adj_m[i].resize( n_i );
        _adj_n_unnorm[i].resize( n_i );
        _adj_m_unnorm[i].resize( n_i );
        _new_adj_n[i].resize( n_i );
        _new_adj_m[i].resize( n_i );
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            { // calculate adj_n
                Prob prod( _fg->factor(I).p() );
                prod *= _adj_b_F_unnorm[I];
                foreach( const Neighbor &j, _fg->nbF(I) )
                    if( i != j ) {
                        Prob n_jI( _bp_dual.msgN( j, j.dual ) );
                        const _ind_t &ind = _index( j, j.dual );
                        // multiply prod with n_jI
                        for( size_t x_I = 0; x_I < prod.size(); x_I++ )
                            prod[x_I] *= n_jI[ind[x_I]];
                    }
                Prob marg( _fg->var(i).states(), 0.0 );
                const _ind_t &ind = _index( i, I.iter );
                for( size_t r = 0; r < prod.size(); r++ )
                    marg[ind[r]] += prod[r];
                _new_adj_n[i][I.iter] = marg;
                upMsgN( i, I.iter );
            }

            { // calculate adj_m
                Prob prod( _adj_b_V_unnorm[i] );
                assert( prod.size() == _fg->var(i).states() );
                foreach( const Neighbor &J, _fg->nbV(i) )
                    if( J.node != I.node )
                        prod *= _bp_dual.msgM(i,J.iter);
                _new_adj_m[i][I.iter] = prod;
                upMsgM( i, I.iter );
            }
        }
    }
}


void BBP::RegenerateSeqMessageAdjoints() {
    size_t nv = _fg->nrVars();
    _adj_m.resize( nv );
    _adj_m_unnorm.resize( nv );
    _new_adj_m.resize( nv );
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        size_t n_i = _fg->nbV(i).size();
        _adj_m[i].resize( n_i );
        _adj_m_unnorm[i].resize( n_i );
        _new_adj_m[i].resize( n_i );
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            // calculate adj_m
            Prob prod( _adj_b_V_unnorm[i] );
            assert( prod.size() == _fg->var(i).states() );
            foreach( const Neighbor &J, _fg->nbV(i) )
                if( J.node != I.node )
                    prod *= _bp_dual.msgM( i, J.iter );
            _adj_m[i][I.iter] = prod;
            calcUnnormMsgM( i, I.iter );
            _new_adj_m[i][I.iter] = Prob( _fg->var(i).states(), 0.0 );
        }
    }
    for( size_t i = 0; i < _fg->nrVars(); i++ ) {
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            // calculate adj_n
            Prob prod( _fg->factor(I).p() );
            prod *= _adj_b_F_unnorm[I];
            foreach( const Neighbor &j, _fg->nbF(I) )
                if( i != j ) {
                    Prob n_jI( _bp_dual.msgN( j, j.dual) );
                    const _ind_t& ind = _index( j, j.dual );
                    // multiply prod with n_jI
                    for( size_t x_I = 0; x_I < prod.size(); x_I++ )
                        prod[x_I] *= n_jI[ind[x_I]];
                }
            Prob marg( _fg->var(i).states(), 0.0 );
            const _ind_t &ind = _index( i, I.iter );
            for( size_t r = 0; r < prod.size(); r++ )
                marg[ind[r]] += prod[r];
            sendSeqMsgN( i, I.iter,marg );
        }
    }
}


void BBP::calcNewN( size_t i, size_t _I ) {
    _adj_psi_V[i] += T(i,_I) * _adj_n_unnorm[i][_I];
    Prob &new_adj_n_iI = _new_adj_n[i][_I];
    new_adj_n_iI = Prob( _fg->var(i).states(), 0.0 );
    size_t I = _fg->nbV(i)[_I];
    foreach( const Neighbor &j, _fg->nbF(I) )
        if( j != i ) {
            const Prob &p = _S[i][_I][j.iter];
            const Prob &_adj_m_unnorm_jI = _adj_m_unnorm[j][j.dual];
            LOOP_ij(
                new_adj_n_iI[xi] += p[xij] * _adj_m_unnorm_jI[xj];
            );
            /* THE FOLLOWING WOULD BE ABOUT TWICE AS SLOW:
            Var vi = _fg->var(i);
            Var vj = _fg->var(j);
            new_adj_n_iI = (Factor(VarSet(vi, vj), p) * Factor(vj,_adj_m_unnorm_jI)).marginal(vi,false).p();
            */
        }
}


void BBP::calcNewM( size_t i, size_t _I ) {
    const Neighbor &I = _fg->nbV(i)[_I];
    Prob p( U(I, I.dual) );
    const Prob &adj = _adj_m_unnorm[i][_I];
    const _ind_t &ind = _index(i,_I);
    for( size_t x_I = 0; x_I < p.size(); x_I++ )
        p[x_I] *= adj[ind[x_I]];
    _adj_psi_F[I] += p;
    /* THE FOLLOWING WOULD BE SLIGHTLY SLOWER:
    _adj_psi_F[I] += (Factor( _fg->factor(I).vars(), U(I, I.dual) ) * Factor( _fg->var(i), _adj_m_unnorm[i][_I] )).p();
    */

    _new_adj_m[i][_I] = Prob( _fg->var(i).states(), 0.0 );
    foreach( const Neighbor &J, _fg->nbV(i) )
        if( J != I )
            _new_adj_m[i][_I] += _R[I][I.dual][J.iter] * _adj_n_unnorm[i][J.iter];
}


void BBP::calcUnnormMsgN( size_t i, size_t _I ) {
    _adj_n_unnorm[i][_I] = unnormAdjoint( _bp_dual.msgN(i,_I), _bp_dual.zN(i,_I), _adj_n[i][_I] );
}


void BBP::calcUnnormMsgM( size_t i, size_t _I ) {
    _adj_m_unnorm[i][_I] = unnormAdjoint( _bp_dual.msgM(i,_I), _bp_dual.zM(i,_I), _adj_m[i][_I] );
}


void BBP::upMsgN( size_t i, size_t _I ) {
    _adj_n[i][_I] = _new_adj_n[i][_I];
    calcUnnormMsgN( i, _I );
}


void BBP::upMsgM( size_t i, size_t _I ) {
    _adj_m[i][_I] = _new_adj_m[i][_I];
    calcUnnormMsgM( i, _I );
}


void BBP::doParUpdate() {
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            calcNewM( i, I.iter );
            calcNewN( i, I.iter );
        }
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            upMsgM( i, I.iter );
            upMsgN( i, I.iter );
        }
}


void BBP::incrSeqMsgM( size_t i, size_t _I, const Prob &p ) {
/*    if( props.clean_updates )
        _new_adj_m[i][_I] += p;
    else {*/
        _adj_m[i][_I] += p;
        calcUnnormMsgM(i, _I);
//    }
}


#if 0
Real pv_thresh=1000;
#endif


/*
void BBP::updateSeqMsgM( size_t i, size_t _I ) {
    if( props.clean_updates ) {
#if 0
        if(_new_adj_m[i][_I].sumAbs() > pv_thresh ||
           _adj_m[i][_I].sumAbs() > pv_thresh) {

            DAI_DMSG("in updateSeqMsgM");
            DAI_PV(i);
            DAI_PV(_I);
            DAI_PV(_adj_m[i][_I]);
            DAI_PV(_new_adj_m[i][_I]);
        }
#endif
        _adj_m[i][_I] += _new_adj_m[i][_I];
        calcUnnormMsgM( i, _I );
        _new_adj_m[i][_I].fill( 0.0 );
    }
}
*/

void BBP::setSeqMsgM( size_t i, size_t _I, const Prob &p ) {
    _adj_m[i][_I] = p;
    calcUnnormMsgM( i, _I );
}


void BBP::sendSeqMsgN( size_t i, size_t _I, const Prob &f ) {
    Prob f_unnorm = unnormAdjoint( _bp_dual.msgN(i,_I), _bp_dual.zN(i,_I), f );
    const Neighbor &I = _fg->nbV(i)[_I];
    assert( I.iter == _I );
    _adj_psi_V[i] += f_unnorm * T( i, _I );
#if 0
    if(f_unnorm.sumAbs() > pv_thresh) {
        DAI_DMSG("in sendSeqMsgN");
        DAI_PV(i);
        DAI_PV(I);
        DAI_PV(_I);
        DAI_PV(_bp_dual.msgN(i,_I));
        DAI_PV(_bp_dual.zN(i,_I));
        DAI_PV(_bp_dual.msgM(i,_I));
        DAI_PV(_bp_dual.zM(i,_I));
        DAI_PV(_fg->factor(I).p());
    }
#endif
    foreach( const Neighbor &J, _fg->nbV(i) ) {
        if( J.node != I.node ) {
#if 0
            if(f_unnorm.sumAbs() > pv_thresh) {
                DAI_DMSG("in sendSeqMsgN loop");
                DAI_PV(J);
                DAI_PV(f_unnorm);
                DAI_PV(_R[J][J.dual][_I]);
                DAI_PV(f_unnorm * _R[J][J.dual][_I]);
            }
#endif
            incrSeqMsgM( i, J.iter, f_unnorm * R(J, J.dual, _I) );
        }
    }
}


void BBP::sendSeqMsgM( size_t j, size_t _I ) {
    const Neighbor &I = _fg->nbV(j)[_I];

//     DAI_PV(j);
//     DAI_PV(I);
//     DAI_PV(_adj_m_unnorm_jI);
//     DAI_PV(_adj_m[j][_I]);
//     DAI_PV(_bp_dual.zM(j,_I));

    size_t _j = I.dual;
    const Prob &_adj_m_unnorm_jI = _adj_m_unnorm[j][_I];
    Prob um( U(I, _j) );
    const _ind_t &ind = _index(j, _I);
    for( size_t x_I = 0; x_I < um.size(); x_I++ )
        um[x_I] *= _adj_m_unnorm_jI[ind[x_I]];
    um *= 1 - props.damping;
    _adj_psi_F[I] += um;
    
    /* THE FOLLOWING WOULD BE SLIGHTLY SLOWER:
    _adj_psi_F[I] += (Factor( _fg->factor(I).vars(), U(I, _j) ) * Factor( _fg->var(j), _adj_m_unnorm[j][_I] )).p() * (1.0 - props.damping);
    */

//     DAI_DMSG("in sendSeqMsgM");
//     DAI_PV(j);
//     DAI_PV(I);
//     DAI_PV(_I);
//     DAI_PV(_fg->nbF(I).size());
    foreach( const Neighbor &i, _fg->nbF(I) ) {
        if( i.node != j ) {
            const Prob &S = _S[i][i.dual][_j];
            Prob msg( _fg->var(i).states(), 0.0 );
            LOOP_ij(
                msg[xi] += S[xij] * _adj_m_unnorm_jI[xj];
            );
            msg *= 1.0 - props.damping;
            /* THE FOLLOWING WOULD BE ABOUT TWICE AS SLOW:
            Var vi = _fg->var(i);
            Var vj = _fg->var(j);
            msg = (Factor(VarSet(vi,vj), S) * Factor(vj,_adj_m_unnorm_jI)).marginal(vi,false).p() * (1.0 - props.damping);
            */
#if 0
            if(msg.sumAbs() > pv_thresh) {
                DAI_DMSG("in sendSeqMsgM loop");

                DAI_PV(j);
                DAI_PV(I);
                DAI_PV(_I);
                DAI_PV(_fg->nbF(I).size());
                DAI_PV(_fg->factor(I).p());
                DAI_PV(_S[i][i.dual][_j]);

                DAI_PV(i);
                DAI_PV(i.dual);
                DAI_PV(msg);
                DAI_PV(_fg->nbV(i).size());
            }
#endif
            assert( _fg->nbV(i)[i.dual].node == I );
            sendSeqMsgN( i, i.dual, msg );
        }
    }
    setSeqMsgM( j, _I, _adj_m[j][_I] * props.damping );
}


Real BBP::getUnMsgMag() {
    Real s = 0.0;
    size_t e = 0;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            s += _adj_m_unnorm[i][I.iter].sumAbs();
            s += _adj_n_unnorm[i][I.iter].sumAbs();
            e++;
        }
    return s / e;
}


void BBP::getMsgMags( Real &s, Real &new_s ) {
    s = 0.0; 
    new_s = 0.0;
    size_t e = 0;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            s += _adj_m[i][I.iter].sumAbs();
            s += _adj_n[i][I.iter].sumAbs();
            new_s += _new_adj_m[i][I.iter].sumAbs();
            new_s += _new_adj_n[i][I.iter].sumAbs();
            e++;
        }
    s /= e; 
    new_s /= e;
}

// tuple<size_t,size_t,Real> BBP::getArgMaxPsi1Adj() {
//     size_t argmax_var=0;
//     size_t argmax_var_state=0;
//     Real max_var=0;
//     for( size_t i = 0; i < _fg->nrVars(); i++ ) {
//         pair<size_t,Real> argmax_state = adj_psi_V(i).argmax();
//         if(i==0 || argmax_state.second>max_var) {
//             argmax_var = i;
//             max_var = argmax_state.second;
//             argmax_var_state = argmax_state.first;
//         }
//     }
//     assert(/*0 <= argmax_var_state &&*/
//            argmax_var_state < _fg->var(argmax_var).states());
//     return tuple<size_t,size_t,Real>(argmax_var,argmax_var_state,max_var);
// }


void BBP::getArgmaxMsgM( size_t &out_i, size_t &out__I, Real &mag ) {
    bool found = false;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) ) {
            Real thisMag = _adj_m[i][I.iter].sumAbs();
            if( !found || mag < thisMag ) {
                found = true;
                mag = thisMag;
                out_i = i;
                out__I = I.iter;
            }
        }
    assert( found );
}


Real BBP::getMaxMsgM() {
    size_t dummy;
    Real mag;
    getArgmaxMsgM( dummy, dummy, mag );
    return mag;
}


Real BBP::getTotalMsgM() {
    Real mag = 0.0;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) )
            mag += _adj_m[i][I.iter].sumAbs();
    return mag;
}


Real BBP::getTotalNewMsgM() {
    Real mag = 0.0;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) )
            mag += _new_adj_m[i][I.iter].sumAbs();
    return mag;
}


Real BBP::getTotalMsgN() {
    Real mag = 0.0;
    for( size_t i = 0; i < _fg->nrVars(); i++ )
        foreach( const Neighbor &I, _fg->nbV(i) )
            mag += _adj_n[i][I.iter].sumAbs();
    return mag;
}


void BBP::Regenerate() {
    RegenerateInds();
    RegenerateT();
    RegenerateU();
    RegenerateS();
    RegenerateR();
    RegenerateInputs();
    RegeneratePsiAdjoints();
    if( props.updates == Properties::UpdateType::PAR )
        RegenerateParMessageAdjoints();
    else
        RegenerateSeqMessageAdjoints();
    _iters = 0;
}


std::vector<Prob> BBP::getZeroAdjF( const FactorGraph &fg ) {
    vector<Prob> adj_2;
    adj_2.reserve( fg.nrFactors() );
    for( size_t I = 0; I < fg.nrFactors(); I++ )
        adj_2.push_back( Prob( fg.factor(I).states(), 0.0 ) );
    return adj_2;
}


std::vector<Prob> BBP::getZeroAdjV( const FactorGraph &fg ) {
    vector<Prob> adj_1;
    adj_1.reserve( fg.nrVars() );
    for( size_t i = 0; i < fg.nrVars(); i++ )
        adj_1.push_back( Prob( fg.var(i).states(), 0.0 ) );
    return adj_1;
}


void BBP::run() {
    typedef BBP::Properties::UpdateType UT;
    Real &tol = props.tol;
    UT &updates = props.updates;

    Real tic = toc();
    switch( (size_t)updates ) {
        case UT::SEQ_MAX: {
            size_t i, _I;
            Real mag;
            do {
                _iters++;
                getArgmaxMsgM( i, _I, mag );
                sendSeqMsgM( i, _I );
            } while( mag > tol && _iters < props.maxiter );

            if( _iters >= props.maxiter )
                if( props.verbose >= 1 )
                    cerr << "Warning: BBP didn't converge in " << _iters << " iterations (greatest message magnitude = " << mag << ")" << endl;
            break;
        } case UT::SEQ_FIX: {
            Real mag;
            do {
                _iters++;
                mag = getTotalMsgM();
                if( mag < tol ) 
                    break;

                for( size_t i = 0; i < _fg->nrVars(); i++ )
                    foreach( const Neighbor &I, _fg->nbV(i) )
                        sendSeqMsgM( i, I.iter );
/*                for( size_t i = 0; i < _fg->nrVars(); i++ )
                    foreach( const Neighbor &I, _fg->nbV(i) )
                        updateSeqMsgM( i, I.iter );*/
            } while( mag > tol && _iters < props.maxiter );

            if( _iters >= props.maxiter )
                if( props.verbose >= 1 )
                    cerr << "Warning: BBP didn't converge in " << _iters << " iterations (greatest message magnitude = " << mag << ")" << endl;
            break;
        } case UT::SEQ_BP_REV:
          case UT::SEQ_BP_FWD: {
            const BP *bp = static_cast<const BP*>(_ia);
            vector<pair<size_t, size_t> > sentMessages = bp->getSentMessages();
            size_t totalMessages = sentMessages.size();
            if( totalMessages == 0 ) {
                cerr << "Asked for updates = " << updates << " but no BP messages; did you forget to set recordSentMessages?" << endl;
                DAI_THROW(INTERNAL_ERROR);
            }
            if( updates==UT::SEQ_BP_FWD )
                reverse( sentMessages.begin(), sentMessages.end() );
//          DAI_PV(sentMessages.size());
//          DAI_PV(_iters);
//          DAI_PV(props.maxiter);
            while( sentMessages.size() > 0 && _iters < props.maxiter ) {
//              DAI_PV(sentMessages.size());
//              DAI_PV(_iters);
                _iters++;
                pair<size_t, size_t> e = sentMessages.back();
                sentMessages.pop_back();
                size_t i = e.first, _I = e.second;
                sendSeqMsgM( i, _I );
            }
            if( _iters >= props.maxiter )
                if( props.verbose >= 1 )
                    cerr << "Warning: BBP updates limited to " << props.maxiter << " iterations, but using UpdateType " << updates << " with " << totalMessages << " messages" << endl;
            break;
        } case UT::PAR: {
            do {
                _iters++;
                doParUpdate();
            } while( (_iters < 2 || getUnMsgMag() > tol) && _iters < props.maxiter );
            if( _iters == props.maxiter ) {
                Real s, new_s;
                getMsgMags( s, new_s );
                if( props.verbose >= 1 )
                    cerr << "Warning: BBP didn't converge in " << _iters << " iterations (unnorm message magnitude = " << getUnMsgMag() << ", norm message mags = " << s << " -> " << new_s << ")" << endl;
            }
            break;
        }
    }
    if( props.verbose >= 3 )
        cerr << "BBP::run() took " << toc()-tic << " seconds " << doneIters() << " iterations" << endl;
}


double numericBBPTest( const InfAlg &bp, const vector<size_t> *state, const PropertySet &bbp_props, bbp_cfn_t cfn, double h ) {
    // calculate the value of the unperturbed cost function
    Real cf0 = getCostFn( bp, cfn, state );

    // run BBP to estimate adjoints
    BBP bbp( &bp, bbp_props );
    initBBPCostFnAdj( bbp, bp, cfn, state );
    bbp.run();

    Real d = 0;
    const FactorGraph& fg = bp.fg();

    if( 1 ) {
        // verify bbp.adj_psi_V

        // for each variable i
        for( size_t i = 0; i < fg.nrVars(); i++ ) {
            vector<double> adj_est;
            // for each value xi
            for( size_t xi = 0; xi < fg.var(i).states(); xi++ ) {
                // Clone 'bp' (which may be any InfAlg)
                InfAlg *bp_prb = bp.clone();

                // perturb it
                size_t n = bp_prb->fg().var(i).states();
                Prob psi_1_prb( n, 1.0 );
                psi_1_prb[xi] += h;
//                 psi_1_prb.normalize();
                size_t I = bp_prb->fg().nbV(i)[0]; // use first factor in list of neighbors of i
                bp_prb->fg().factor(I) *= Factor( bp_prb->fg().var(i), psi_1_prb );
                
                // call 'init' on the perturbed variables
                bp_prb->init( bp_prb->fg().var(i) );
                
                // run copy to convergence
                bp_prb->run();

                // calculate new value of cost function
                Real cf_prb = getCostFn( *bp_prb, cfn, state );

                // use to estimate adjoint for i
                adj_est.push_back( (cf_prb - cf0) / h );
                
                // free cloned InfAlg
                delete bp_prb;
            }
            Prob p_adj_est( adj_est.begin(), adj_est.end() );
            // compare this numerical estimate to the BBP estimate; sum the distances
            cout << "i: " << i
                 << ", p_adj_est: " << p_adj_est
                 << ", bbp.adj_psi_V(i): " << bbp.adj_psi_V(i) << endl;
            d += dist( p_adj_est, bbp.adj_psi_V(i), Prob::DISTL1 );
        }
    }
    /*    if(1) {
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
                    Real cf_prb = getCostFn(bp_dual_prb, cfn, &state);
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
                    Real cf_prb = getCostFn(bp_dual_prb, cfn, &state);
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
    */
    /*    if(0) {
        // verify bbp.adj_b_V
        for(size_t i=0; i<bp_dual.nrVars(); i++) {
            vector<double> adj_b_V_est;
            // for each value xi
            for(size_t xi=0; xi<bp_dual.var(i).states(); xi++) {
                BP_dual bp_dual_prb(bp_dual);

                // make h-sized change to b_1(i)[x_i]
                bp_dual_prb._beliefs.b1[i][xi] += h;

                // get cost function value
                Real cf_prb = getCostFn(bp_dual_prb, cfn, &state);

                // add it to list of adjoints
                adj_b_V_est.push_back((cf_prb-cf0)/h);
            }
            Prob p_adj_b_V_est(adj_b_V_est.begin(), adj_b_V_est.end());
            // compare this numerical estimate to the BBP estimate; sum the distances
            cerr << "i: " << i
                 << ", adj_b_V_est: " << p_adj_b_V_est
                 << ", bbp.adj_b_V(i): " << bbp.adj_b_V(i) << endl;
            d += dist(p_adj_b_V_est, bbp.adj_b_V(i), Prob::DISTL1);
        }
    }
    */

    // return total of distances
    return d;
}


bool needGibbsState( bbp_cfn_t cfn ) {
    switch( (size_t)cfn ) {
        case bbp_cfn_t::CFN_GIBBS_B:
        case bbp_cfn_t::CFN_GIBBS_B2:
        case bbp_cfn_t::CFN_GIBBS_EXP:
        case bbp_cfn_t::CFN_GIBBS_B_FACTOR:
        case bbp_cfn_t::CFN_GIBBS_B2_FACTOR:
        case bbp_cfn_t::CFN_GIBBS_EXP_FACTOR:
            return true;
        default:
            return false;
    }
}


void initBBPCostFnAdj( BBP &bbp, const InfAlg &ia, bbp_cfn_t cfn_type, const vector<size_t> *stateP ) {
    const FactorGraph &fg = ia.fg();
    
    switch( (size_t)cfn_type ) {
        case bbp_cfn_t::CFN_BETHE_ENT: {
            vector<Prob> b1_adj;
            vector<Prob> b2_adj;
            vector<Prob> psi1_adj;
            vector<Prob> psi2_adj;
            b1_adj.reserve( fg.nrVars() );
            psi1_adj.reserve( fg.nrVars() );
            b2_adj.reserve( fg.nrFactors() );
            psi2_adj.reserve( fg.nrFactors() );
            for( size_t i = 0; i < fg.nrVars(); i++ ) {
                size_t dim = fg.var(i).states();
                int c = fg.nbV(i).size();
                Prob p(dim,0.0);
                for( size_t xi = 0; xi < dim; xi++ )
                    p[xi] = (1 - c) * (1 + log( ia.beliefV(i)[xi] ));
                b1_adj.push_back( p );

                for( size_t xi = 0; xi < dim; xi++ )
                    p[xi] = -ia.beliefV(i)[xi];
                psi1_adj.push_back( p );
            }
            for( size_t I = 0; I < fg.nrFactors(); I++ ) {
                size_t dim = fg.factor(I).states();
                Prob p( dim, 0.0 );
                for( size_t xI = 0; xI < dim; xI++ )
                    p[xI] = 1 + log( ia.beliefF(I)[xI] / fg.factor(I).p()[xI] );
                b2_adj.push_back( p );

                for( size_t xI = 0; xI < dim; xI++ )
                    p[xI] = -ia.beliefF(I)[xI] / fg.factor(I).p()[xI];
                psi2_adj.push_back( p );
            }
            bbp.init( b1_adj, b2_adj, psi1_adj, psi2_adj );
            break;
        } case bbp_cfn_t::CFN_FACTOR_ENT: {
            vector<Prob> b2_adj;
            b2_adj.reserve( fg.nrFactors() );
            for( size_t I = 0; I < fg.nrFactors(); I++ ) {
                size_t dim = fg.factor(I).states();
                Prob p( dim, 0.0 );
                for( size_t xI = 0; xI < dim; xI++ ) {
                    double bIxI = ia.beliefF(I)[xI];
                    if( bIxI < 1.0e-15 )
                        p[xI] = -1.0e10;
                    else
                        p[xI] = 1 + log( bIxI );
                }
                b2_adj.push_back(p);
            }
            bbp.init( bbp.getZeroAdjV(fg), b2_adj );
            break;
        } case bbp_cfn_t::CFN_VAR_ENT: {
            vector<Prob> b1_adj;
            b1_adj.reserve( fg.nrVars() );
            for( size_t i = 0; i < fg.nrVars(); i++ ) {
                size_t dim = fg.var(i).states();
                Prob p( dim, 0.0 );
                for( size_t xi = 0; xi < fg.var(i).states(); xi++ ) {
                    double bixi = ia.beliefV(i)[xi];
                    if( bixi < 1.0e-15 )
                        p[xi] = -1.0e10;
                    else
                        p[xi] = 1 + log( bixi );
                }
                b1_adj.push_back( p );
            }
            bbp.init( b1_adj );
            break;
        } case bbp_cfn_t::CFN_GIBBS_B:
          case bbp_cfn_t::CFN_GIBBS_B2:
          case bbp_cfn_t::CFN_GIBBS_EXP: {
            // cost functions that use Gibbs sample, summing over variable marginals
            vector<size_t> state;
            if( stateP == NULL )
                state = getGibbsState( ia, 2*ia.Iterations() );
            else
                state = *stateP;
            assert( state.size() == fg.nrVars() );

            vector<Prob> b1_adj;
            b1_adj.reserve(fg.nrVars());
            for( size_t i = 0; i < state.size(); i++ ) {
                size_t n = fg.var(i).states();
                Prob delta( n, 0.0 );
                assert(/*0<=state[i] &&*/ state[i] < n);
                double b = ia.beliefV(i)[state[i]];
                switch( (size_t)cfn_type ) {
                    case bbp_cfn_t::CFN_GIBBS_B:
                        delta[state[i]] = 1.0; 
                        break;
                    case bbp_cfn_t::CFN_GIBBS_B2:
                        delta[state[i]] = b; 
                        break;
                    case bbp_cfn_t::CFN_GIBBS_EXP:
                        delta[state[i]] = exp(b); 
                        break;
                    default: 
                        DAI_THROW(INTERNAL_ERROR);
                }
                b1_adj.push_back( delta );
            }
            bbp.init( b1_adj );
            break;
        } case bbp_cfn_t::CFN_GIBBS_B_FACTOR:
          case bbp_cfn_t::CFN_GIBBS_B2_FACTOR:
          case bbp_cfn_t::CFN_GIBBS_EXP_FACTOR: {
            // cost functions that use Gibbs sample, summing over factor marginals
            vector<size_t> state;
            if( stateP == NULL )
                state = getGibbsState( ia, 2*ia.Iterations() );
            else
                state = *stateP;
            assert( state.size() == fg.nrVars() );

            vector<Prob> b2_adj;
            b2_adj.reserve( fg.nrVars() );
            for( size_t I = 0; I <  fg.nrFactors(); I++ ) {
                size_t n = fg.factor(I).states();
                Prob delta( n, 0.0 );

                size_t x_I = getFactorEntryForState( fg, I, state );
                assert(/*0<=x_I &&*/ x_I < n);

                double b = ia.beliefF(I)[x_I];
                switch( (size_t)cfn_type ) {
                    case bbp_cfn_t::CFN_GIBBS_B_FACTOR:
                        delta[x_I] = 1.0; 
                        break;
                    case bbp_cfn_t::CFN_GIBBS_B2_FACTOR:
                        delta[x_I] = b; 
                        break;
                    case bbp_cfn_t::CFN_GIBBS_EXP_FACTOR:
                        delta[x_I] = exp( b );
                        break;
                    default: 
                        DAI_THROW(INTERNAL_ERROR);
                }
                b2_adj.push_back( delta );
            }
            bbp.init( bbp.getZeroAdjV(fg), b2_adj );
            break;
        } default: 
            DAI_THROW(UNKNOWN_ENUM_VALUE);
    }
}


Real getCostFn( const InfAlg &ia, bbp_cfn_t cfn_type, const vector<size_t> *stateP ) {
    double cf = 0.0;
    const FactorGraph &fg = ia.fg();

    switch( (size_t)cfn_type ) {
        case bbp_cfn_t::CFN_BETHE_ENT: // ignores state
            cf = -ia.logZ();
            break;
        case bbp_cfn_t::CFN_VAR_ENT: // ignores state
            for( size_t i = 0; i < fg.nrVars(); i++ )
                cf += -ia.beliefV(i).entropy();
            break;
        case bbp_cfn_t::CFN_FACTOR_ENT: // ignores state
            for( size_t I = 0; I < fg.nrFactors(); I++ )
                cf += -ia.beliefF(I).entropy();
            break;
        case bbp_cfn_t::CFN_GIBBS_B:
        case bbp_cfn_t::CFN_GIBBS_B2:
        case bbp_cfn_t::CFN_GIBBS_EXP: {
            assert( stateP != NULL );
            vector<size_t> state = *stateP;
            assert( state.size() == fg.nrVars() );
            for( size_t i = 0; i < fg.nrVars(); i++ ) {
                double b = ia.beliefV(i)[state[i]];
                switch( (size_t)cfn_type ) {
                    case bbp_cfn_t::CFN_GIBBS_B:
                        cf += b;
                        break;
                    case bbp_cfn_t::CFN_GIBBS_B2:
                        cf += b * b / 2.0;
                        break;
                    case bbp_cfn_t::CFN_GIBBS_EXP:
                        cf += exp( b );
                        break;
                    default:
                        DAI_THROW(INTERNAL_ERROR);
                }
            }
            break;
        } case bbp_cfn_t::CFN_GIBBS_B_FACTOR:
          case bbp_cfn_t::CFN_GIBBS_B2_FACTOR:
          case bbp_cfn_t::CFN_GIBBS_EXP_FACTOR: {
            assert( stateP != NULL );
            vector<size_t> state = *stateP;
            assert( state.size() == fg.nrVars() );
            for( size_t I = 0; I < fg.nrFactors(); I++ ) {
                size_t x_I = getFactorEntryForState( fg, I, state );
                double b = ia.beliefF(I)[x_I];
                switch( (size_t)cfn_type ) {
                    case bbp_cfn_t::CFN_GIBBS_B_FACTOR:
                        cf += b;
                        break;
                    case bbp_cfn_t::CFN_GIBBS_B2_FACTOR:
                        cf += b * b / 2.0;
                        break;
                    case bbp_cfn_t::CFN_GIBBS_EXP_FACTOR:
                        cf += exp( b );
                        break;
                    default:
                        DAI_THROW(INTERNAL_ERROR);
                }
            }
            break;
        } default: 
            cerr << "Unknown cost function " << cfn_type << endl;
            DAI_THROW(UNKNOWN_ENUM_VALUE);
    }
    return cf;
}


} // end of namespace dai


/* {{{ GENERATED CODE: DO NOT EDIT. Created by 
    ./scripts/regenerate-properties include/dai/bbp.h src/bbp.cpp 
*/
namespace dai {

void BBP::Properties::set(const PropertySet &opts)
{
    const std::set<PropertyKey> &keys = opts.allKeys();
    std::set<PropertyKey>::const_iterator i;
    bool die=false;
    for(i=keys.begin(); i!=keys.end(); i++) {
        if(*i == "verbose") continue;
        if(*i == "maxiter") continue;
        if(*i == "tol") continue;
        if(*i == "damping") continue;
        if(*i == "updates") continue;
        cerr << "BBP: Unknown property " << *i << endl;
        die=true;
    }
    if(die) {
        DAI_THROW(UNKNOWN_PROPERTY_TYPE);
    }
    if(!opts.hasKey("verbose")) {
        cerr << "BBP: Missing property \"verbose\" for method \"BBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("maxiter")) {
        cerr << "BBP: Missing property \"maxiter\" for method \"BBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("tol")) {
        cerr << "BBP: Missing property \"tol\" for method \"BBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("damping")) {
        cerr << "BBP: Missing property \"damping\" for method \"BBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("updates")) {
        cerr << "BBP: Missing property \"updates\" for method \"BBP\"" << endl;
        die=true;
    }
    if(die) {
        DAI_THROW(NOT_ALL_PROPERTIES_SPECIFIED);
    }
    verbose = opts.getStringAs<size_t>("verbose");
    maxiter = opts.getStringAs<size_t>("maxiter");
    tol = opts.getStringAs<double>("tol");
    damping = opts.getStringAs<double>("damping");
    updates = opts.getStringAs<UpdateType>("updates");
}
PropertySet BBP::Properties::get() const {
    PropertySet opts;
    opts.Set("verbose", verbose);
    opts.Set("maxiter", maxiter);
    opts.Set("tol", tol);
    opts.Set("damping", damping);
    opts.Set("updates", updates);
    return opts;
}
string BBP::Properties::toString() const {
    stringstream s(stringstream::out);
    s << "[";
    s << "verbose=" << verbose << ",";
    s << "maxiter=" << maxiter << ",";
    s << "tol=" << tol << ",";
    s << "damping=" << damping << ",";
    s << "updates=" << updates << ",";
    s << "]";
    return s.str();
}
} // end of namespace dai
/* }}} END OF GENERATED CODE */
