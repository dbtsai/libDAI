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


#include <iostream>
#include <sstream>
#include <algorithm>

#include <dai/bp_dual.h>
#include <dai/util.h>
#include <dai/bipgraph.h>


namespace dai {


using namespace std;


typedef BipartiteGraph::Neighbor Neighbor;


void BP_dual::regenerateMessages() {
    size_t nv = fg().nrVars();
    _msgs.Zn.resize(nv);
    _msgs.Zm.resize(nv);
    _msgs.m.resize(nv);
    _msgs.n.resize(nv);
    for( size_t i = 0; i < nv; i++ ) {
        size_t nvf = fg().nbV(i).size();
        _msgs.Zn[i].resize(nvf, 1.0);
        _msgs.Zm[i].resize(nvf, 1.0);
        size_t states = fg().var(i).states();
        _msgs.n[i].resize(nvf, Prob(states));
        _msgs.m[i].resize(nvf, Prob(states));
    }
}


void BP_dual::regenerateBeliefs() {
    _beliefs.b1.clear();
    _beliefs.b1.reserve(fg().nrVars());
    _beliefs.Zb1.resize(fg().nrVars(), 1.0);
    _beliefs.b2.clear();
    _beliefs.b2.reserve(fg().nrFactors());
    _beliefs.Zb2.resize(fg().nrFactors(), 1.0);

    for( size_t i = 0; i < fg().nrVars(); i++ )
        _beliefs.b1.push_back( Prob( fg().var(i).states() ) );
    for( size_t I = 0; I < fg().nrFactors(); I++ )
        _beliefs.b2.push_back( Prob( fg().factor(I).states() ) );
}


void BP_dual::init() {
    regenerateMessages();
    regenerateBeliefs();
    calcMessages();
    calcBeliefs();
}


void BP_dual::calcMessages() {
    // calculate 'n' messages from "factor marginal / factor"
    vector<Factor> bs;
    size_t nf = fg().nrFactors();
    for( size_t I = 0; I < nf; I++ )
        bs.push_back(_ia->beliefF(I));
    assert(nf == bs.size());
    for( size_t I = 0; I < nf; I++ ) {
        Factor f = bs[I];
        f /= fg().factor(I);
        foreach(const Neighbor &i, fg().nbF(I))
            msgN(i, i.dual) = f.marginal(fg().var(i)).p();
    }
    // calculate 'm' messages and normalizers from 'n' messages
    for( size_t i = 0; i < fg().nrVars(); i++ )
        foreach(const Neighbor &I, fg().nbV(i))
            calcNewM(i, I.iter);
    // recalculate 'n' messages and normalizers from 'm' messages
    for( size_t i = 0; i < fg().nrVars(); i++ ) {
        foreach(const Neighbor &I, fg().nbV(i)) {
            Prob oldN = msgN(i,I.iter);
            calcNewN(i, I.iter);
            Prob newN = msgN(i,I.iter);
#if 0
            // check that new 'n' messages match old ones
            if((oldN-newN).maxAbs() > 1.0e-5) {
                cerr << "New 'n' messages don't match old: " <<
                    "(i,I) = (" << i << ", " << I << 
                    ") old = " << oldN << ", new = " << newN << endl;
                DAI_THROW(INTERNAL_ERROR);
            }
#endif
        }
    }
}


void BP_dual::calcBeliefV(size_t i) {
    Prob prod( fg().var(i).states(), 1.0 );
    foreach(const Neighbor &I, fg().nbV(i))
        prod *= msgM(i,I.iter);
    _beliefs.Zb1[i] = prod.normalize();
    _beliefs.b1[i] = prod;
}


void BP_dual::calcBeliefF(size_t I) {
    Prob prod( fg().factor(I).p() );
    foreach(const Neighbor &j, fg().nbF(I)) {
        IndexFor ind (fg().var(j), fg().factor(I).vars() );
        Prob n(msgN(j,j.dual));
        for(size_t x=0; ind >= 0; x++, ++ind)
            prod[x] *= n[ind];
    }
    _beliefs.Zb2[I] = prod.normalize();
    _beliefs.b2[I] = prod;
}


// called after run()
void BP_dual::calcBeliefs() {
    for( size_t i = 0; i < fg().nrVars(); i++ )
        calcBeliefV(i);  // calculate b_i
    for( size_t I = 0; I < fg().nrFactors(); I++ )
        calcBeliefF(I);  // calculate b_I
}


void BP_dual::calcNewM(size_t i, size_t _I) {
    // calculate updated message I->i
    const Neighbor &I = fg().nbV(i)[_I];
    Prob prod( fg().factor(I).p() );
    foreach(const Neighbor &j, fg().nbF(I)) {
        if( j != i ) {     // for all j in I \ i
            Prob n(msgN(j,j.dual));
            IndexFor ind(fg().var(j), fg().factor(I).vars());
            for(size_t x=0; ind >= 0; x++, ++ind)
                prod[x] *= n[ind];
        }
    }
    // Marginalize onto i
    Prob marg( fg().var(i).states(), 0.0 );
    // ind is the precalculated Index(i,I) i.e. to x_I == k corresponds x_i == ind[k]
    IndexFor ind(fg().var(i), fg().factor(I).vars());
    for(size_t x=0; ind >= 0; x++, ++ind)
        marg[ind] += prod[x];
    
    _msgs.Zm[i][_I] = marg.normalize();
    _msgs.m[i][_I] = marg;
}


void BP_dual::calcNewN(size_t i, size_t _I) {
    // calculate updated message i->I
    const Neighbor &I = fg().nbV(i)[_I];
    Prob prod(fg().var(i).states(), 1.0);
    foreach(const Neighbor &J, fg().nbV(i)) {
        if(J.node != I.node) // for all J in i \ I
            prod *= msgM(i,J.iter);
    }
    _msgs.Zn[i][_I] = prod.normalize();
    _msgs.n[i][_I] = prod;
}


} // end of namespace dai
