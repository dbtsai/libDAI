/*
  Copyright 2009 Charles Vaske <cvaske@soe.ucsc.edu>
  University of California Santa Cruz

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __defined_libdai_emalg_h
#define __defined_libdai_emalg_h

#include<vector>
#include<map>

#include <dai/factor.h>
#include <dai/daialg.h>
#include <dai/evidence.h>
#include <dai/index.h>
#include <dai/properties.h>


/// \file 
/** \brief Defines classes related to Expectation Maximization:
 *  EMAlg, ParameterEstimate, and FactorOrientations
 */

namespace dai {

///Interface for a parameter estimation method. 
/** This parameter estimation interface is based on sufficient statistics. 
 *  Implementations are responsible for collecting data from a probability 
 *  vector passed to it from a SharedParameters container object.
 *
 *  Implementations of this interface should register a factory function
 *  via the static ParameterEstimation::registerMethod function.
 */
class ParameterEstimation {
public:  
  /// A pointer to a factory function.
  typedef ParameterEstimation* (*ParamEstFactory)(const PropertySet&);

  /// General factory method for construction of ParameterEstimation subclasses.
  static ParameterEstimation* construct(const std::string& method, 
					const PropertySet& p);
  /// Register a subclass with ParameterEstimation::construct.
  static void registerMethod(const std::string method, 
			     const ParamEstFactory f) {
    if (_registry == NULL) {
      loadDefaultRegistry();
    }
    (*_registry)[method] = f;
  }
  /// Virtual destructor for deleting pointers to derived classes.
  virtual ~ParameterEstimation() {}
  /// Estimate the factor using the accumulated sufficient statistics and reset.
  virtual Prob estimate() = 0;
  /// Accumulate the sufficient statistics for p.
  virtual void addSufficientStatistics(Prob& p) = 0;
  /// Returns the size of the Prob that is passed to addSufficientStatistics.
  virtual size_t probSize() const = 0;
  /// A virtual copy constructor.
  virtual ParameterEstimation* clone() const= 0;
private:
  static std::map< std::string, ParamEstFactory >* _registry;
  static void loadDefaultRegistry();
};

/// Estimates the parameters of a conditional probability, using pseudocounts.
class CondProbEstimation : private ParameterEstimation {
private:
  size_t _target_dim;
  Prob _stats;
  Prob _initial_stats;
public:
  /** For a conditional probability \f$ Pr( X | Y ) \f$, 
   *  \param target_dimension should equal \f$ | X | \f$
   *  \param pseudocounts has length \f$ |X| \cdot |Y| \f$
   */
  CondProbEstimation(size_t target_dimension, Prob pseudocounts);

  /// Virtual constructor, using a PropertySet.
  /** Some keys in the PropertySet are required:
   *     - target_dimension, which should be equal to \f$ | X | \f$
   *     - total_dimension, which sholud be equal to \f$ |X| \cdot |Y| \f$
   *  
   *  An optional key is:
   *     - pseudo_count which specifies the initial counts (defaults to 1)
   */
  static ParameterEstimation* factory(const PropertySet& p);
  /// Virtual destructor
  virtual ~CondProbEstimation() {}
  /// Returns an estimate of the conditional probability distribution.
  /** The format of the resulting Prob keeps all the values for 
   *  \f$ P(X | Y=a) \f$ sequential in teh array.
   */
  virtual Prob estimate();
  /// Accumulate sufficient statistics from the expectations in p.
  virtual void addSufficientStatistics(Prob& p);
  /// Returns the required size for arguments to addSufficientStatistics
  virtual size_t probSize() const { return _stats.size(); }
  /// Virtual copy constructor.
  virtual ParameterEstimation* clone() const {
    return new CondProbEstimation(_target_dim, _initial_stats);
  }
};

/** A single factor or set of factors whose parameters should be
 *  estimated.  Each factor's values are reordered to match a
 *  canonical variable ordering.  This canonical variable ordering
 *  will likely not be the order of variables required to make two
 *  factors parameters isomorphic.  Therefore, this ordering of the
 *  variables must be specified for ever factor to ensure that
 *  parameters can be shared between different factors during EM.
 */
class SharedParameters {
public:
  /// Convenience label for an index into a FactorGraph to a factor.
  typedef size_t FactorIndex;
  /// Convenience label for a grouping of factor orientations.
  typedef std::map< FactorIndex, std::vector< Var > > FactorOrientations;
private:
  std::map< FactorIndex, VarSet > _varsets;
  std::map< FactorIndex, Permute > _perms;
  FactorOrientations _varorders;
  ParameterEstimation* _estimation;
  bool _deleteEstimation;

  static Permute calculatePermutation(const std::vector< Var >& varorder,
				      const std::vector< size_t >& dims,
				      VarSet& outVS);
  void setPermsAndVarSetsFromVarOrders();
public:
  /// Copy constructor
  SharedParameters(const SharedParameters& sp);
  /// Constructor useful in programmatic settings 
  /** \param varorders  all the factor orientations for this parameter
      \param estimation a pointer to the parameter estimation method
   */ 
  SharedParameters(const FactorOrientations& varorders,
		   ParameterEstimation* estimation);

  /// Constructor for making an object from a stream
  SharedParameters(std::istream& is, const FactorGraph& fg_varlookup);

  /// Destructor
  ~SharedParameters() { if (_deleteEstimation) delete _estimation; }

  /// Collect the necessary statistics from expected values
  void collectSufficientStatistics(InfAlg& alg);

  /// Estimate and set the shared parameters
  void setParameters(FactorGraph& fg);
};

/** A maximization step groups together several parameter estimation
 * tasks into a single unit.
 */
class MaximizationStep { 
private:
  std::vector< SharedParameters > _params;
public:
  MaximizationStep() : _params() {}

  /// Construct an step object taht contains all these estimation probelms
  MaximizationStep(std::vector< SharedParameters >& maximizations) : 
    _params(maximizations) {}  

  /// Construct a step from an input stream
  MaximizationStep(std::istream& is, const FactorGraph& fg_varlookup);
  
  /** Collect the beliefs from this InfAlg as expectations for
   *  the next Maximization step.
   */
  void addExpectations(InfAlg& alg);

  /** Using all of the currently added expectations, make new factors 
   *  with maximized parameters and set them in the FactorGraph.
   */
  void maximize(FactorGraph& fg);
};

/// EMAlg performs Expectation Maximization to learn factor parameters.
/** This requires specifying:
 *     - Evidence (instances of observations from the graphical model),
 *     - InfAlg for performing the E-step, which includes the factor graph
 *     - a vector of MaximizationSteps steps to be performed
 *
 *  This implementation can peform incremental EM by using multiple 
 *  MaximizationSteps.  An expectation step is performed between execution
 *  of each MaximizationStep.  A call to iterate() will cycle through all
 *  MaximizationSteps.
 */  
class EMAlg {
private:
  /// All the data samples used during learning
  const Evidence& _evidence;
  
  /// How to do the expectation step
  InfAlg& _estep;

  /// The maximization steps to take
  std::vector<MaximizationStep> _msteps;

  size_t _iters;
  std::vector<Real> _lastLogZ;

public:
  /// Construct an EMAlg from all these objects
  EMAlg(const Evidence& evidence, InfAlg& estep, 
	std::vector<MaximizationStep>& msteps) 
    : _evidence(evidence),
      _estep(estep),
      _msteps(msteps),
      _iters(0),
      _lastLogZ() 
  {}
  
  /// Construct an EMAlg from an input stream
  EMAlg(const Evidence& evidence, InfAlg& estep, std::istream& mstep_file);

  /// Perform an iteration over all maximization steps
  Real iterate();
  /// Performs an iteration over a single MaximizationStep
  Real iterate(const MaximizationStep& mstep);

};

} // namespace dai

#endif
