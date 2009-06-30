#include <dai/util.h>

#include <dai/emalg.h>

namespace dai{

std::map< std::string, ParameterEstimation::ParamEstFactory>* 
ParameterEstimation::_registry = NULL;

void ParameterEstimation::loadDefaultRegistry() {
  _registry = new std::map< std::string, ParamEstFactory>();
  (*_registry)["ConditionalProbEstimation"] = CondProbEstimation::factory;
}

ParameterEstimation* ParameterEstimation::construct(const std::string& method,
						    const PropertySet& p) {
  if (_registry == NULL) {
    loadDefaultRegistry();
  }
  std::map< std::string, ParamEstFactory>::iterator i = _registry->find(method);
  if (i == _registry->end()) {
    DAI_THROW(UNKNOWN_PARAMETER_ESTIMATION_METHOD);
  }
  ParamEstFactory factory = i->second;
  return factory(p);
}

ParameterEstimation* CondProbEstimation::factory(const PropertySet& p) {
  size_t target_dimension =  p.getStringAs<size_t>("target_dim");
  size_t total_dimension = p.getStringAs<size_t>("total_dim");
  Real pseudo_count = 1;
  if (p.hasKey("pseudo_count")) {
    pseudo_count = p.getStringAs<Real>("pseudo_count");
  }
  Prob counts_vec(total_dimension, pseudo_count);
  return new CondProbEstimation(target_dimension, counts_vec);
}

CondProbEstimation::CondProbEstimation(size_t target_dimension, 
				       Prob pseudocounts) 
  : _target_dim(target_dimension),
    _stats(pseudocounts),
    _initial_stats(pseudocounts) {
  if (_stats.size() % _target_dim) {
    DAI_THROW(MALFORMED_PROPERTY);
  }
}

void CondProbEstimation::addSufficientStatistics(Prob& p) {
  _stats += p;
}

Prob CondProbEstimation::estimate() {
  for (size_t parent = 0; parent < _stats.size(); parent += _target_dim) {
    Real norm = 0;
    size_t top = parent + _target_dim;
    for (size_t i = parent; i < top; ++i) {
      norm += _stats[i];
    }
    if (norm != 0) {
      norm = 1 / norm;
    }
    for (size_t i = parent; i < top; ++i) {
      _stats[i] *= norm;
    }
  }
  Prob result = _stats;
  _stats = _initial_stats;
  return result;
}

Permute
SharedParameters::calculatePermutation(const std::vector< Var >& varorder,
				       const std::vector< size_t >& dims,
				       VarSet& outVS) {
  std::vector<long> labels(dims.size());
  
  // Check that the variable set is compatible
  if (varorder.size() != dims.size()) {
    DAI_THROW(INVALID_SHARED_PARAMETERS_ORDER);
  }
  
  // Collect all labels, and order them in vs
  for (size_t di = 0; di < dims.size(); ++di) {
    if (dims[di] != varorder[di].states()) {
      DAI_THROW(INVALID_SHARED_PARAMETERS_ORDER);
    }
    outVS |= varorder[di];
    labels[di] = varorder[di].label();
  }
  
  // Construct the sigma array for the permutation object
  std::vector<size_t> sigma(dims.size(), 0);
  VarSet::iterator set_iterator = outVS.begin();
  for  (size_t vs_i = 0; vs_i < dims.size(); ++vs_i, ++set_iterator) {
    std::vector< long >::iterator location = find(labels.begin(), labels.end(),
						  set_iterator->label());
    sigma[vs_i] = location - labels.begin();
  }
  
  return Permute(dims, sigma);
}

void SharedParameters::setPermsAndVarSetsFromVarOrders() {
  if (_varorders.size() == 0) {
    return;
  }
  FactorOrientations::const_iterator foi = _varorders.begin();
  std::vector< size_t > dims(foi->second.size());
  size_t total_dim = 1;
  for (size_t i = 0; i < dims.size(); ++i) {
    dims[i] = foi->second[i].states();
    total_dim *= dims[i];
  }
  
  // Construct the permutation objects
  for ( ; foi != _varorders.end(); ++foi) {
    VarSet vs;
    _perms[foi->first] = calculatePermutation(foi->second, dims, vs);
    _varsets[foi->first] = vs;
  }
  
  if  (_estimation == NULL || _estimation->probSize() != total_dim) {
    DAI_THROW(INVALID_SHARED_PARAMETERS_ORDER);
  }
}

SharedParameters::SharedParameters(std::istream& is,
				   const FactorGraph& fg_varlookup)
  : _varsets(),
    _perms(),
    _varorders(),
    _estimation(NULL),
    _deleteEstimation(1) 
{
  std::string est_method;
  PropertySet props;
  is >> est_method;
  is >> props;

  _estimation = ParameterEstimation::construct(est_method, props);

  size_t num_factors;
  is >> num_factors;
  for (size_t sp_i = 0; sp_i < num_factors; ++sp_i) {
    std::string line;
    std::vector< std::string > fields;
    size_t factor;
    std::vector< Var > var_order;
    std::istringstream iss;

    while(line.size() == 0 && getline(is, line));
    tokenizeString(line, fields, " \t");

    // Lookup the factor in the factorgraph
    if (fields.size() < 1) { 
      DAI_THROW(INVALID_SHARED_PARAMETERS_INPUT_LINE);
    }
    iss.str(fields[0]);
    iss >> factor;
    const VarSet& vs = fg_varlookup.factor(factor).vars();
    if (fields.size() != vs.size() + 1) {
      DAI_THROW(INVALID_SHARED_PARAMETERS_INPUT_LINE);
    }

    // Construct the vector of Vars
    for (size_t fi = 1; fi < fields.size(); ++fi) {
      // Lookup a single variable by label
      long label;
      std::istringstream labelparse(fields[fi]);
      labelparse >> label;
      VarSet::const_iterator vsi = vs.begin();
      for ( ; vsi != vs.end(); ++vsi) {
	if (vsi->label() == label) break;
      }
      if (vsi == vs.end()) {
	DAI_THROW(INVALID_SHARED_PARAMETERS_INPUT_LINE);
      }
      var_order.push_back(*vsi);
    }
    _varorders[factor] = var_order;
  }
  setPermsAndVarSetsFromVarOrders();
}

SharedParameters::SharedParameters(const SharedParameters& sp) 
  : _varsets(sp._varsets),
    _perms(sp._perms),
    _varorders(sp._varorders),
    _estimation(sp._estimation),
    _deleteEstimation(sp._deleteEstimation)
{
  if (_deleteEstimation) {
    _estimation = _estimation->clone();
  }
}

SharedParameters::SharedParameters(const FactorOrientations& varorders,
				   ParameterEstimation* estimation) 
  : _varsets(),
    _perms(),
    _varorders(varorders),
    _estimation(estimation),
    _deleteEstimation(0) 
{
  setPermsAndVarSetsFromVarOrders();
}

void SharedParameters::collectSufficientStatistics(InfAlg& alg) {
  std::map< FactorIndex, Permute >::iterator i = _perms.begin();
  for ( ; i != _perms.end(); ++i) {
    Permute& perm = i->second;
    VarSet& vs = _varsets[i->first];
    
    Factor b = alg.belief(vs);
    Prob p(b.states(), 0.0);
    for (size_t entry = 0; entry < b.states(); ++entry) {
      p[entry] = b[perm.convert_linear_index(entry)];
    }
    _estimation->addSufficientStatistics(p);
  }
}

void SharedParameters::setParameters(FactorGraph& fg) {
  Prob p = _estimation->estimate();
  std::map< FactorIndex, Permute >::iterator i = _perms.begin();
  for ( ; i != _perms.end(); ++i) {
    Permute& perm = i->second;
    VarSet& vs = _varsets[i->first];
    
    Factor f(vs, 0.0);
    for (size_t entry = 0; entry < f.states(); ++entry) {
      f[perm.convert_linear_index(entry)] = p[entry];
    }

    fg.setFactor(i->first, f);
  }
}

MaximizationStep::MaximizationStep (std::istream& is,
				    const FactorGraph& fg_varlookup ) 
  : _params()
{
  size_t num_params = -1;
  is >> num_params;
  _params.reserve(num_params);
  for (size_t i = 0; i < num_params; ++i) {
    SharedParameters p(is, fg_varlookup);
    _params.push_back(p);
  }
}


void MaximizationStep::addExpectations(InfAlg& alg) {
  for (size_t i = 0; i < _params.size(); ++i) {
    _params[i].collectSufficientStatistics(alg);
  }
}

void MaximizationStep::maximize(FactorGraph& fg) {
  for (size_t i = 0; i < _params.size(); ++i) {
    _params[i].setParameters(fg);
  }
}

const std::string EMAlg::MAX_ITERS_KEY("max_iters");
const std::string EMAlg::LOG_Z_TOL_KEY("log_z_tol");
const size_t EMAlg::MAX_ITERS_DEFAULT = 30;
const Real EMAlg::LOG_Z_TOL_DEFAULT = 0.01;

EMAlg::EMAlg(const Evidence& evidence, InfAlg& estep, std::istream& msteps_file)
  : _evidence(evidence),
    _estep(estep),
    _msteps(),
    _iters(0),
    _lastLogZ(),
    _max_iters(MAX_ITERS_DEFAULT),
    _log_z_tol(LOG_Z_TOL_DEFAULT)
{
  msteps_file.exceptions( std::istream::eofbit | std::istream::failbit 
			  | std::istream::badbit );
  size_t num_msteps = -1;
  msteps_file >> num_msteps;
  _msteps.reserve(num_msteps);
  for (size_t i = 0; i < num_msteps; ++i) {
    MaximizationStep m(msteps_file, estep.fg());
    _msteps.push_back(m);
  }
}	

void EMAlg::setTermConditions(const PropertySet* p) {
  if (NULL == p) {
    return;
  }
  if (p->hasKey(MAX_ITERS_KEY)) {
    _max_iters = p->getStringAs<size_t>(MAX_ITERS_KEY);
  }
  if (p->hasKey(LOG_Z_TOL_KEY)) {
    _log_z_tol = p->getStringAs<Real>(LOG_Z_TOL_KEY);
  }
}

bool EMAlg::hasSatisfiedTermConditions() const {
  if (_iters >= _max_iters) {
    return 1;
  } else if (_lastLogZ.size() < 3) { 
    // need at least 2 to calculate ratio
    // Also, throw away first iteration, as the parameters may not
    // have been normalized according to the estimation method
    return 0;
  } else {
    Real current = _lastLogZ[_lastLogZ.size() - 1];
    Real previous = _lastLogZ[_lastLogZ.size() - 2];
    if (previous == 0) return 0;
    Real diff = current - previous;
    if (diff < 0) {
      std::cerr << "Error: in EM log-likehood decreased from " << previous 
		<< " to " << current << std::endl;
      return 1;
    }
    return diff / abs(previous) <= _log_z_tol;
  }
}

Real EMAlg::iterate(MaximizationStep& mstep) {
  Evidence::const_iterator e = _evidence.begin();
  Real logZ = 0;

  // Expectation calculation
  for ( ; e != _evidence.end(); ++e) {
    InfAlg* clamped = _estep.clone();
    e->second.applyEvidence(*clamped);
    clamped->run();
    
    logZ += clamped->logZ();

    mstep.addExpectations(*clamped);

    delete clamped;
  }
  
  // Maximization of parameters
  mstep.maximize(_estep.fg());

  return logZ;
}

Real EMAlg::iterate() {
  Real likelihood;
  for (size_t i = 0; i < _msteps.size(); ++i) {
    likelihood = iterate(_msteps[i]);
  }
  _lastLogZ.push_back(likelihood);
  ++_iters;
  return likelihood;
}

void EMAlg::run() {
  while(!hasSatisfiedTermConditions()) {
    iterate();
  }
}

}
