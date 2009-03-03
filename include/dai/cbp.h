#ifndef ____defined_libdai_cbp_h__
#define ____defined_libdai_cbp_h__

#include <dai/daialg.h>

#include <dai/cbp.h>
#include <dai/bbp.h>

using namespace std;

/* Class for "clamped belief propagation".

'chooseNextClamp' can be overridden (e.g. in BPC); this version
randomly chooses variables to clamp.
*/

namespace dai {

template<class T>
vector<T> concat(const vector<T>& u, const vector<T>& v) {
    vector<T> w;
    w.reserve(u.size()+v.size());
    for(size_t i=0; i<u.size(); i++)
        w.push_back(u[i]);
    for(size_t i=0; i<v.size(); i++)
        w.push_back(v[i]);
    return w;
}

class CBP : public DAIAlgFG {
    vector<Factor> _beliefs1;
    vector<Factor> _beliefs2;
    double _logZ;
    double _est_logZ;
    double _old_est_logZ;

    double _sum_level;
    size_t _num_leaves;

    size_t maxClampLevel() { return props.max_levels; }
    double minMaxAdj() { return props.min_max_adj; }
    double recTol() { return props.rec_tol; }


    bbp_cfn_t BBP_cost_function() {
      return props.bbp_cfn;
    }

    void printDebugInfo();

    void runRecurse(BP_dual &bp_dual,
                    double orig_logZ,
                    vector<size_t> clamped_vars_list,
                    set<size_t> clamped_vars,
                    size_t &num_leaves,
                    double &sum_level,
                    Real &lz_out, vector<Factor>& beliefs_out);

    BP getBP();
    BP_dual getBP_dual();

  public:
    size_t _iters;
    double _maxdiff;

    struct Properties {
      typedef BP::Properties::UpdateType UpdateType;
      DAI_ENUM(RecurseType,REC_FIXED,REC_LOGZ,REC_BDIFF);
      DAI_ENUM(ChooseMethodType,CHOOSE_RANDOM,CHOOSE_BBP,CHOOSE_BP_ALL);
      DAI_ENUM(ClampType,CLAMP_VAR,CLAMP_FACTOR);

      double tol;
      double rec_tol;
      size_t maxiter;
      size_t verbose;
      UpdateType updates;
      size_t max_levels;
      double min_max_adj;

      ChooseMethodType choose;
      RecurseType recursion;
      ClampType clamp;
      bbp_cfn_t bbp_cfn;
      size_t rand_seed;
    } props;

    Properties::ChooseMethodType ChooseMethod() { return props.choose; }
    Properties::RecurseType Recursion() { return props.recursion; }
    Properties::ClampType Clamping() { return props.clamp; }

    //----------------------------------------------------------------

    // construct CBP object from FactorGraph
    CBP(const FactorGraph &fg, const PropertySet &opts) : DAIAlgFG(fg) {
        setProperties(opts);
        construct();
    }

    static const char *Name;
    /// List of property names
    static const char *PropertyList[];

    /// @name General InfAlg interface
    //@{
    virtual CBP* clone() const { return new CBP(*this); }
//     virtual CBP* create() const { return new CBP(); }
    virtual CBP* create() const { throw "Unimplemented"; }
    virtual std::string identify() const { return string(Name) + printProperties(); }
    virtual Factor belief (const Var &n) const { return _beliefs1[findVar(n)]; }
    virtual Factor belief (const VarSet &) const { assert(0); }
    virtual vector<Factor> beliefs() const { return concat(_beliefs1, _beliefs2); }
    virtual Real logZ() const { return _logZ; }
    virtual void init() {};
    virtual void init( const VarSet & ) {};
    virtual double run();
    virtual double maxDiff() const { return _maxdiff; }
    virtual size_t Iterations() const { return _iters; }
    //@}


    //----------------------------------------------------------------

    virtual bool chooseNextClampVar(BP_dual& bp,
                                    vector<size_t> &clamped_vars_list,
                                    set<size_t> &clamped_vars,
                                    size_t &i, vector<size_t> &xis, Real *maxVarOut);

    //----------------------------------------------------------------

    void setBeliefsFrom(const BP& b);
    void setBeliefs(const vector<Factor>& bs, double logZ) {
        size_t i=0;
        _beliefs1.clear(); _beliefs1.reserve(nrVars());
        _beliefs2.clear(); _beliefs2.reserve(nrFactors());
        for(i=0; i<nrVars(); i++) _beliefs1.push_back(bs[i]);
        for(; i<nrVars()+nrFactors(); i++) _beliefs2.push_back(bs[i]);
        _logZ = logZ;
    }
    Factor belief1 (size_t i) const { return _beliefs1[i]; }
    Factor belief2 (size_t I) const { return _beliefs2[I]; }

    //----------------------------------------------------------------

    void construct();

    void setProperties( const PropertySet &opts );
    PropertySet getProperties() const;
    std::string printProperties() const;
};

pair<size_t, size_t> bbpFindClampVar(BP_dual &in_bp_dual, bool clampVar, size_t numClamped, bbp_cfn_t cfn, const PropertySet &props, Real *maxVarOut);

}

#endif
