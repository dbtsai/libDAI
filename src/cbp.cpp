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
#include <map>
#include <set>
#include <algorithm>

#include <dai/util.h>
#include <dai/properties.h>

#include <dai/bp.h>
#include <dai/cbp.h>
#include <dai/bbp.h>


namespace dai {


using namespace std;
using boost::shared_ptr;


const char *CBP::Name = "CBP";


void CBP::setBeliefs( const std::vector<Factor> &bs, double logZ ) {
    size_t i = 0;
    _beliefsV.clear(); 
    _beliefsV.reserve( nrVars() );
    _beliefsF.clear(); 
    _beliefsF.reserve( nrFactors() );
    for( i = 0; i < nrVars(); i++ ) 
        _beliefsV.push_back( bs[i] );
    for( ; i < nrVars() + nrFactors(); i++ ) 
        _beliefsF.push_back( bs[i] );
    _logZ = logZ;
}


void CBP::construct() {
    _beliefsV.clear();
    _beliefsV.reserve(nrVars());
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefsV.push_back( Factor(var(i)).normalized() );

    _beliefsF.clear();
    _beliefsF.reserve(nrFactors());
    for( size_t I = 0; I < nrFactors(); I++ ) {
        Factor f = factor(I);
        f.fill(1); f.normalize();
        _beliefsF.push_back( f );
    }

    // to compute average level
    _sum_level = 0;
    _num_leaves = 0;

    _maxdiff = 0;
    _iters = 0;

    if( props.clamp_outfile.length() > 0 ) {
        _clamp_ofstream = shared_ptr<ofstream>(new ofstream( props.clamp_outfile.c_str(), ios_base::out|ios_base::trunc ));
        *_clamp_ofstream << "# COUNT LEVEL VAR STATE" << endl;
    }
}


/// Calculates a vector of mixtures p * b + (1-p) * c
static vector<Factor> mixBeliefs( Real p, const vector<Factor> &b, const vector<Factor> &c ) {
    vector<Factor> out;
    assert( b.size() == c.size() );
    out.reserve( b.size() );
    Real pc = 1 - p;
    for( size_t i = 0; i < b.size(); i++ )
        // probably already normalized, but do it again just in case
        out.push_back( b[i].normalized() * p + c[i].normalized() * pc );
    return out;
}


double CBP::run() {
    size_t seed = props.rand_seed;
    if( seed > 0) 
        rnd_seed( seed );

    InfAlg *bp = getInfAlg();
    bp->init();
    bp->run();
    _iters += bp->Iterations();

    vector<Factor> beliefs_out;
    Real lz_out;
    size_t choose_count=0;
    runRecurse( bp, bp->logZ(), vector<size_t>(0), _num_leaves, choose_count, _sum_level, lz_out, beliefs_out );
    if( props.verbose >= 1 )
        cerr << "CBP average levels = " << (_sum_level / _num_leaves) << ", leaves = " << _num_leaves << endl;
    setBeliefs( beliefs_out, lz_out );
    return 0.0;
}


/// \todo Eventually this allow other inference algorithms that BP to be selected, via a property
InfAlg* CBP::getInfAlg() {
    PropertySet bpProps;
    bpProps.Set("updates", props.updates);
    bpProps.Set("tol", props.tol);
    bpProps.Set("maxiter", props.maxiter);
    bpProps.Set("verbose", props.verbose);
    bpProps.Set("logdomain", false);
    bpProps.Set("damping", 0.0);
    BP *bp = new BP( *this, bpProps );
    bp->recordSentMessages = true;
    bp->init();
    return bp;
}


void CBP::runRecurse( InfAlg *bp, double orig_logZ, vector<size_t> clamped_vars_list, size_t &num_leaves,
                      size_t &choose_count, double &sum_level, Real &lz_out, vector<Factor>& beliefs_out) {
    // choose a variable/states to clamp:
    size_t i;
    vector<size_t> xis;
    Real maxVar = 0.0;
    bool found;
    bool clampingVar = (Clamping() == Properties::ClampType::CLAMP_VAR);

    if( Recursion() == Properties::RecurseType::REC_LOGZ && recTol() > 0 && exp( bp->logZ() - orig_logZ ) < recTol() )
        found = false;
    else
        found = chooseNextClampVar( bp, clamped_vars_list, i, xis, &maxVar );

    if( !found ) {
        num_leaves++;
        sum_level += clamped_vars_list.size();
        beliefs_out = bp->beliefs();
        lz_out = bp->logZ();
        return;
    }

    choose_count++;
    if( props.clamp_outfile.length() > 0 )
        *_clamp_ofstream << choose_count << "\t" << clamped_vars_list.size() << "\t" << i << "\t" << xis[0] << endl;
    
    if( clampingVar )
        foreach( size_t xi, xis ) 
            assert(/*0<=xi &&*/ xi < var(i).states() );
    else
        foreach( size_t xI, xis ) 
            assert(/*0<=xI &&*/ xI < factor(i).states() );
    // - otherwise, clamp and recurse, saving margin estimates for each
    // clamp setting. afterwards, combine estimates.

    // compute complement of 'xis'
    vector<size_t> cmp_xis = complement( xis, clampingVar ? var(i).states() : factor(i).states() );

    /// \todo could do this more efficiently with a nesting version of saveProbs/undoProbs
    Real lz; 
    vector<Factor> b;
    InfAlg *bp_c = bp->clone();
    if( clampingVar ) {
        bp_c->fg().clampVar( i, xis );
        bp_c->init( var(i) );
    } else {
        bp_c->fg().clampFactor( i, xis );
        bp_c->init( factor(i).vars() );
    }
    bp_c->run();
    _iters += bp_c->Iterations();

    lz = bp_c->logZ();
    b = bp_c->beliefs();

    Real cmp_lz; 
    vector<Factor> cmp_b;
    InfAlg *cmp_bp_c = bp->clone();
    if( clampingVar ) {
        cmp_bp_c->fg().clampVar( i, cmp_xis );
        cmp_bp_c->init(var(i));
    } else {
        cmp_bp_c->fg().clampFactor( i, cmp_xis );
        cmp_bp_c->init( factor(i).vars() );
    }
    cmp_bp_c->run();
    _iters += cmp_bp_c->Iterations();

    cmp_lz = cmp_bp_c->logZ();
    cmp_b = cmp_bp_c->beliefs();

    double p = unSoftMax( lz, cmp_lz );
    Real bp__d = 0.0;
    
    if( Recursion() == Properties::RecurseType::REC_BDIFF && recTol() > 0 ) {
        vector<Factor> combined_b( mixBeliefs( p, b, cmp_b ) );
        Real new_lz = logSumExp( lz,cmp_lz );
        bp__d = dist( bp->beliefs(), combined_b, nrVars() );
        if( exp( new_lz - orig_logZ) * bp__d < recTol() ) {
            num_leaves++;
            sum_level += clamped_vars_list.size();
            beliefs_out = combined_b;
            lz_out = new_lz;
            return;
        }
    }

    // either we are not doing REC_BDIFF or the distance was large
    // enough to recurse:
    runRecurse( bp_c, orig_logZ, clamped_vars_list, num_leaves, choose_count, sum_level, lz, b );
    runRecurse( cmp_bp_c, orig_logZ, clamped_vars_list, num_leaves, choose_count, sum_level, cmp_lz, cmp_b );

    p = unSoftMax( lz, cmp_lz );

    beliefs_out = mixBeliefs( p, b, cmp_b );
    lz_out = logSumExp( lz, cmp_lz );

    if( props.verbose >= 2 ) {
        Real d = dist( bp->beliefs(), beliefs_out, nrVars() );
        cerr << "Distance (clamping " << i << "): " << d;
        if( Recursion() == Properties::RecurseType::REC_BDIFF )
            cerr << "; bp_dual predicted " << bp__d;
        cerr << "; max_adjoint = " << maxVar << "; logZ = " << lz_out << " (in " << bp->logZ() << ") (orig " << orig_logZ << "); p = " << p << "; level = " << clamped_vars_list.size() << endl;
    }

    delete bp_c;
    delete cmp_bp_c;
}


// 'xis' must be sorted
bool CBP::chooseNextClampVar( InfAlg *bp, vector<size_t> &clamped_vars_list, size_t &i, vector<size_t> &xis, Real *maxVarOut ) {
    Real tiny = 1.0e-14;
    if( props.verbose >= 3 )
        cerr << "clamped_vars_list" << clamped_vars_list << endl;
    if( clamped_vars_list.size() >= maxClampLevel() )
        return false;
    if( ChooseMethod() == Properties::ChooseMethodType::CHOOSE_RANDOM ) {
        if( Clamping() == Properties::ClampType::CLAMP_VAR ) {
            int t = 0, t1 = 100;
            do {
                i = rnd( nrVars() );
                t++;
            } while( abs( bp->beliefV(i).p().max() - 1) < tiny && t < t1 );
            if( t == t1 ) {
                return false;
                // die("Too many levels requested in CBP");
            }
            // only pick probable values for variable
            size_t xi;
            do {
                xi = rnd( var(i).states() );
                t++;
            } while( bp->beliefV(i).p()[xi] < tiny && t < t1 );
            assert( t < t1 );
            xis.resize( 1, xi );
            //         assert(!_clamped_vars.count(i)); // not true for >2-ary variables
            DAI_IFVERB(2, "CHOOSE_RANDOM at level " << clamped_vars_list.size() << " chose variable " << i << " state " << xis[0] << endl);
        } else {
            int t = 0, t1 = 100;
            do {
                i = rnd( nrFactors() );
                t++;
            } while( abs( bp->beliefF(i).p().max() - 1) < tiny && t < t1 );
            if( t == t1 )
                return false;
                // die("Too many levels requested in CBP");
            // only pick probable values for variable
            size_t xi;
            do {
                xi = rnd( factor(i).states() );
                t++;
            } while( bp->beliefF(i).p()[xi] < tiny && t < t1 );
            assert( t < t1 );
            xis.resize( 1, xi );
            //         assert(!_clamped_vars.count(i)); // not true for >2-ary variables
            DAI_IFVERB(2, endl<<"CHOOSE_RANDOM chose factor "<<i<<" state "<<xis[0]<<endl);
        }
    } else if( ChooseMethod() == Properties::ChooseMethodType::CHOOSE_MAXENT ) {
        if( Clamping() == Properties::ClampType::CLAMP_VAR ) {
            Real max_ent = -1.0;
            int win_k = -1, win_xk = -1;
            for( size_t k = 0; k < nrVars(); k++ ) {
                Real ent=bp->beliefV(k).entropy();
                if( max_ent < ent ) {
                    max_ent = ent;
                    win_k = k;
                    win_xk = bp->beliefV(k).p().argmax().first;
                }
            }
            assert( win_k >= 0 ); 
            assert( win_xk >= 0 );
            i = win_k; 
            xis.resize( 1, win_xk );
            DAI_IFVERB(2, endl<<"CHOOSE_MAXENT chose variable "<<i<<" state "<<xis[0]<<endl);
            if( bp->beliefV(i).p()[xis[0]] < tiny ) {
                DAI_IFVERB(2, "Warning: CHOOSE_MAXENT found unlikely state, not recursing");
                return false;
            }
        } else {
            Real max_ent = -1.0;
            int win_k = -1, win_xk = -1;
            for( size_t k = 0; k < nrFactors(); k++ ) {
                Real ent = bp->beliefF(k).entropy();
                if( max_ent < ent ) {
                    max_ent = ent;
                    win_k = k;
                    win_xk = bp->beliefF(k).p().argmax().first;
                }
            }
            assert( win_k >= 0 ); 
            assert( win_xk >= 0 );
            i = win_k; 
            xis.resize( 1, win_xk );
            DAI_IFVERB(2, endl<<"CHOOSE_MAXENT chose factor "<<i<<" state "<<xis[0]<<endl);
            if( bp->beliefF(i).p()[xis[0]] < tiny ) {
                DAI_IFVERB(2, "Warning: CHOOSE_MAXENT found unlikely state, not recursing");
                return false;
            }
        }
    } else if( ChooseMethod()==Properties::ChooseMethodType::CHOOSE_BP_L1 ||
               ChooseMethod()==Properties::ChooseMethodType::CHOOSE_BP_CFN ) {
        bool doL1 = (ChooseMethod() == Properties::ChooseMethodType::CHOOSE_BP_L1);
        vector<size_t> state;
        if( !doL1 && needGibbsState( BBP_cost_function() ) )
            state = getGibbsState( *bp, 2*bp->Iterations() );
        // try clamping each variable manually
        assert( Clamping() == Properties::ClampType::CLAMP_VAR );
        Real max_cost = 0.0;
        int win_k = -1, win_xk = -1;
        for( size_t k = 0; k < nrVars(); k++ ) {
            for( size_t xk = 0; xk < var(k).states(); xk++ ) {
                if( bp->beliefV(k)[xk] < tiny ) 
                    continue;
                InfAlg *bp1 = bp->clone();
                bp1->clamp( k, xk );
                bp1->init( var(k) );
                bp1->run();
                Real cost = 0;
                if( doL1 )
                    for( size_t j = 0; j < nrVars(); j++ )
                        cost += dist( bp->beliefV(j), bp1->beliefV(j), Prob::DISTL1 );
                else
                    cost = getCostFn( *bp1, BBP_cost_function(), &state );
                if( cost > max_cost || win_k == -1 ) {
                    max_cost = cost;
                    win_k = k;
                    win_xk = xk;
                }
                delete bp1;
            }
        }
        assert( win_k >= 0 ); 
        assert( win_xk >= 0 );
        i = win_k; 
        xis.resize( 1, win_xk );
    } else if( ChooseMethod() == Properties::ChooseMethodType::CHOOSE_BBP ) {
        Real mvo; 
        if( !maxVarOut ) 
            maxVarOut = &mvo;
        bool clampingVar = (Clamping() == Properties::ClampType::CLAMP_VAR);
        pair<size_t, size_t> cv = bbpFindClampVar( *bp, clampingVar, props.bbp_props, BBP_cost_function(), &mvo );

        // if slope isn't big enough then don't clamp
        if( mvo < minMaxAdj() )
            return false;

        size_t xi = cv.second;
        i = cv.first;
#define VAR_INFO (clampingVar?"variable ":"factor ")                       \
            << i << " state " << xi                                     \
            << " (p=" << (clampingVar?bp->beliefV(i)[xi]:bp->beliefF(i)[xi]) \
            << ", entropy = " << (clampingVar?bp->beliefV(i):bp->beliefF(i)).entropy() \
            << ", maxVar = "<< mvo << ")" 
        Prob b = ( clampingVar ? bp->beliefV(i).p() : bp->beliefF(i).p());
        if( b[xi] < tiny ) {
            cerr << "Warning, at level " << clamped_vars_list.size() << ", bbpFindClampVar found unlikely " << VAR_INFO << endl;
            return false;
        }
        if( abs(b[xi] - 1) < tiny ) {
            cerr << "Warning at level " << clamped_vars_list.size() << ", bbpFindClampVar found overly likely " << VAR_INFO << endl;
            return false;
        }

        xis.resize( 1, xi );
        if( clampingVar )
            assert(/*0<=xi &&*/ xi < var(i).states() );
        else
            assert(/*0<=xi &&*/ xi < factor(i).states() );
        DAI_IFVERB(2, "CHOOSE_BBP (num clamped = " << clamped_vars_list.size() << ") chose " << i << " state " << xi << endl);
    } else
        DAI_THROW(INTERNAL_ERROR);
    clamped_vars_list.push_back( i );
    return true;
}


void CBP::printDebugInfo() {
    DAI_PV(_beliefsV);
    DAI_PV(_beliefsF);
    DAI_PV(_logZ);
}


//----------------------------------------------------------------


/** Function which takes a factor graph as input, runs Gibbs and BP_dual,
 *  creates and runs a BBP object, finds best variable, returns
 *  (variable,state) pair for clamping
 */
pair<size_t, size_t> bbpFindClampVar( const InfAlg &in_bp, bool clampingVar, const PropertySet &bbp_props, bbp_cfn_t cfn, Real *maxVarOut ) {
    BBP bbp( &in_bp, bbp_props );
    initBBPCostFnAdj( bbp, in_bp, cfn, NULL );
    bbp.run();
  
    // find and return the (variable,state) with the largest adj_psi_V
    size_t argmax_var = 0;
    size_t argmax_var_state = 0;
    Real max_var = 0;
    if( clampingVar ) {
        for( size_t i = 0; i < in_bp.fg().nrVars(); i++ ) {
            Prob adj_psi_V = bbp.adj_psi_V(i);
            if(0) {
                // helps to account for amount of movement possible in variable
                // i's beliefs? seems not..
                adj_psi_V *= in_bp.beliefV(i).entropy();
            }
            if(0){
//                 adj_psi_V *= Prob(in_bp.fg().var(i).states(),1.0)-in_bp.beliefV(i).p();
                adj_psi_V *= in_bp.beliefV(i).p();
            }
            // try to compensate for effect on same variable (doesn't work)
            //     adj_psi_V[gibbs.state()[i]] -= bp_dual.beliefV(i)[gibbs.state()[i]]/10;
            pair<size_t,Real> argmax_state = adj_psi_V.argmax();

            if( i == 0 || argmax_state.second > max_var ) {
                argmax_var = i;
                max_var = argmax_state.second;
                argmax_var_state = argmax_state.first;
            }
        }
        assert(/*0 <= argmax_var_state &&*/
               argmax_var_state < in_bp.fg().var(argmax_var).states() );
    } else {
        for( size_t I = 0; I < in_bp.fg().nrFactors(); I++ ) {
            Prob adj_psi_F = bbp.adj_psi_F(I);
            if(0) {
                // helps to account for amount of movement possible in variable
                // i's beliefs? seems not..
                adj_psi_F *= in_bp.beliefF(I).entropy();
            }
            // try to compensate for effect on same variable (doesn't work)
            //     adj_psi_V[gibbs.state()[i]] -= bp_dual.beliefV(i)[gibbs.state()[i]]/10;
            pair<size_t,Real> argmax_state = adj_psi_F.argmax();

            if( I == 0 || argmax_state.second > max_var ) {
                argmax_var = I;
                max_var = argmax_state.second;
                argmax_var_state = argmax_state.first;
            }
        }
        assert(/*0 <= argmax_var_state &&*/
               argmax_var_state < in_bp.fg().factor(argmax_var).states() );
    }
    if( maxVarOut ) 
        *maxVarOut = max_var;
    return make_pair( argmax_var, argmax_var_state );
}


vector<size_t> complement( vector<size_t> &xis, size_t n_states ) {
    vector<size_t> cmp_xis( 0 );
    size_t j = 0;
    for( size_t xi = 0; xi < n_states; xi++ ) {
        while( j < xis.size() && xis[j] < xi )
            j++;
        if( j >= xis.size() || xis[j] > xi )
            cmp_xis.push_back(xi);
    }
    assert( xis.size()+cmp_xis.size() == n_states );
    return cmp_xis;
}


Real unSoftMax( Real a, Real b ) {
    if( a > b )
        return 1.0 / (1.0 + exp(b-a));
    else
        return exp(a-b) / (exp(a-b) + 1.0);
}


Real logSumExp( Real a, Real b ) {
    if( a > b )
        return a + log1p( exp( b-a ) );
    else
        return b + log1p( exp( a-b ) );
}


Real dist( const vector<Factor> &b1, const vector<Factor> &b2, size_t nv ) {
    Real d = 0.0;
    for( size_t k = 0; k < nv; k++ )
        d += dist( b1[k], b2[k], Prob::DISTLINF );
    return d;
}


} // end of namespace dai


/* {{{ GENERATED CODE: DO NOT EDIT. Created by 
    ./scripts/regenerate-properties include/dai/cbp.h src/cbp.cpp 
*/
namespace dai {

void CBP::Properties::set(const PropertySet &opts)
{
    const std::set<PropertyKey> &keys = opts.allKeys();
    std::set<PropertyKey>::const_iterator i;
    bool die=false;
    for(i=keys.begin(); i!=keys.end(); i++) {
        if(*i == "verbose") continue;
        if(*i == "tol") continue;
        if(*i == "updates") continue;
        if(*i == "maxiter") continue;
        if(*i == "rec_tol") continue;
        if(*i == "max_levels") continue;
        if(*i == "min_max_adj") continue;
        if(*i == "choose") continue;
        if(*i == "recursion") continue;
        if(*i == "clamp") continue;
        if(*i == "bbp_props") continue;
        if(*i == "bbp_cfn") continue;
        if(*i == "rand_seed") continue;
        if(*i == "clamp_outfile") continue;
        cerr << "CBP: Unknown property " << *i << endl;
        die=true;
    }
    if(die) {
        DAI_THROW(UNKNOWN_PROPERTY_TYPE);
    }
    if(!opts.hasKey("tol")) {
        cerr << "CBP: Missing property \"tol\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("updates")) {
        cerr << "CBP: Missing property \"updates\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("maxiter")) {
        cerr << "CBP: Missing property \"maxiter\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("rec_tol")) {
        cerr << "CBP: Missing property \"rec_tol\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("min_max_adj")) {
        cerr << "CBP: Missing property \"min_max_adj\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("choose")) {
        cerr << "CBP: Missing property \"choose\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("recursion")) {
        cerr << "CBP: Missing property \"recursion\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("clamp")) {
        cerr << "CBP: Missing property \"clamp\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("bbp_props")) {
        cerr << "CBP: Missing property \"bbp_props\" for method \"CBP\"" << endl;
        die=true;
    }
    if(!opts.hasKey("bbp_cfn")) {
        cerr << "CBP: Missing property \"bbp_cfn\" for method \"CBP\"" << endl;
        die=true;
    }
    if(die) {
        DAI_THROW(NOT_ALL_PROPERTIES_SPECIFIED);
    }
    if(opts.hasKey("verbose")) {
        verbose = opts.getStringAs<size_t>("verbose");
    } else {
        verbose = 0;
    }
    tol = opts.getStringAs<double>("tol");
    updates = opts.getStringAs<UpdateType>("updates");
    maxiter = opts.getStringAs<size_t>("maxiter");
    rec_tol = opts.getStringAs<double>("rec_tol");
    if(opts.hasKey("max_levels")) {
        max_levels = opts.getStringAs<size_t>("max_levels");
    } else {
        max_levels = 10;
    }
    min_max_adj = opts.getStringAs<double>("min_max_adj");
    choose = opts.getStringAs<ChooseMethodType>("choose");
    recursion = opts.getStringAs<RecurseType>("recursion");
    clamp = opts.getStringAs<ClampType>("clamp");
    bbp_props = opts.getStringAs<PropertySet>("bbp_props");
    bbp_cfn = opts.getStringAs<bbp_cfn_t>("bbp_cfn");
    if(opts.hasKey("rand_seed")) {
        rand_seed = opts.getStringAs<size_t>("rand_seed");
    } else {
        rand_seed = 0;
    }
    if(opts.hasKey("clamp_outfile")) {
        clamp_outfile = opts.getStringAs<std::string>("clamp_outfile");
    } else {
        clamp_outfile = "";
    }
}
PropertySet CBP::Properties::get() const {
    PropertySet opts;
    opts.Set("verbose", verbose);
    opts.Set("tol", tol);
    opts.Set("updates", updates);
    opts.Set("maxiter", maxiter);
    opts.Set("rec_tol", rec_tol);
    opts.Set("max_levels", max_levels);
    opts.Set("min_max_adj", min_max_adj);
    opts.Set("choose", choose);
    opts.Set("recursion", recursion);
    opts.Set("clamp", clamp);
    opts.Set("bbp_props", bbp_props);
    opts.Set("bbp_cfn", bbp_cfn);
    opts.Set("rand_seed", rand_seed);
    opts.Set("clamp_outfile", clamp_outfile);
    return opts;
}
string CBP::Properties::toString() const {
    stringstream s(stringstream::out);
    s << "[";
    s << "verbose=" << verbose << ",";
    s << "tol=" << tol << ",";
    s << "updates=" << updates << ",";
    s << "maxiter=" << maxiter << ",";
    s << "rec_tol=" << rec_tol << ",";
    s << "max_levels=" << max_levels << ",";
    s << "min_max_adj=" << min_max_adj << ",";
    s << "choose=" << choose << ",";
    s << "recursion=" << recursion << ",";
    s << "clamp=" << clamp << ",";
    s << "bbp_props=" << bbp_props << ",";
    s << "bbp_cfn=" << bbp_cfn << ",";
    s << "rand_seed=" << rand_seed << ",";
    s << "clamp_outfile=" << clamp_outfile;
    s << "]";
    return s.str();
}
} // end of namespace dai
/* }}} END OF GENERATED CODE */
