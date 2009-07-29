/*  Copyright (C) 2009  Charles Vaske  [cvaske at soe dot ucsc dot edu]
    University of California Santa Cruz

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


#ifndef __defined_libdai_emalg_h
#define __defined_libdai_emalg_h


#include <vector>
#include <map>

#include <dai/factor.h>
#include <dai/daialg.h>
#include <dai/evidence.h>
#include <dai/index.h>
#include <dai/properties.h>


/// \file
/// \brief Defines classes related to Expectation Maximization: EMAlg, ParameterEstimation, CondProbEstimation and SharedParameters
/// \todo Describe EM file format


namespace dai {


/// Interface for a parameter estimation method. 
/** This parameter estimation interface is based on sufficient statistics. 
 *  Implementations are responsible for collecting data from a probability 
 *  vector passed to it from a SharedParameters container object.
 *
 *  Implementations of this interface should register a factory function
 *  via the static ParameterEstimation::registerMethod function.
 *  The default registry only contains CondProbEstimation, named
 *  "ConditionalProbEstimation".
 */
class ParameterEstimation {
    public:  
        /// A pointer to a factory function.
        typedef ParameterEstimation* (*ParamEstFactory)( const PropertySet& );

        /// Virtual destructor for deleting pointers to derived classes.
        virtual ~ParameterEstimation() {}
        /// Virtual copy constructor.
        virtual ParameterEstimation* clone() const = 0;

        /// General factory method for construction of ParameterEstimation subclasses.
        static ParameterEstimation* construct( const std::string &method, const PropertySet &p );

        /// Register a subclass so that it can be used with ParameterEstimation::construct.
        static void registerMethod( const std::string &method, const ParamEstFactory &f ) {
            if( _registry == NULL )
                loadDefaultRegistry();
            (*_registry)[method] = f;
        }

        /// Estimate the factor using the accumulated sufficient statistics and reset.
        virtual Prob estimate() = 0;

        /// Accumulate the sufficient statistics for p.
        virtual void addSufficientStatistics( const Prob &p ) = 0;

        /// Returns the size of the Prob that should be passed to addSufficientStatistics.
        virtual size_t probSize() const = 0;

    private:
        /// A static registry containing all methods registered so far.
        static std::map<std::string, ParamEstFactory> *_registry;

        /// Registers default ParameterEstimation subclasses (currently, only CondProbEstimation).
        static void loadDefaultRegistry();
};


/// Estimates the parameters of a conditional probability table, using pseudocounts.
class CondProbEstimation : private ParameterEstimation {
    private:
        /// Number of states of the variable of interest
        size_t _target_dim;
        /// Current pseudocounts
        Prob _stats;
        /// Initial pseudocounts
        Prob _initial_stats;

    public:
        /// Constructor
        /** For a conditional probability \f$ P( X | Y ) \f$, 
         *  \param target_dimension should equal \f$ | X | \f$
         *  \param pseudocounts has length \f$ |X| \cdot |Y| \f$
         */
        CondProbEstimation( size_t target_dimension, const Prob &pseudocounts );

        /// Virtual constructor, using a PropertySet.
        /** Some keys in the PropertySet are required.
         *  For a conditional probability \f$ P( X | Y ) \f$, 
         *     - target_dimension should be equal to \f$ | X | \f$
         *     - total_dimension should be equal to \f$ |X| \cdot |Y| \f$
         *  
         *  An optional key is:
         *     - pseudo_count, which specifies the initial counts (defaults to 1)
         */
        static ParameterEstimation* factory( const PropertySet &p );
        
        /// Virtual copy constructor
        virtual ParameterEstimation* clone() const { return new CondProbEstimation( _target_dim, _initial_stats ); }

        /// Virtual destructor
        virtual ~CondProbEstimation() {}
        
        /// Returns an estimate of the conditional probability distribution.
        /** The format of the resulting Prob keeps all the values for 
         *  \f$ P(X | Y=y) \f$ in sequential order in the array.
         */
        virtual Prob estimate();
        
        /// Accumulate sufficient statistics from the expectations in p.
        virtual void addSufficientStatistics( const Prob &p );
        
        /// Returns the required size for arguments to addSufficientStatistics.
        virtual size_t probSize() const { return _stats.size(); }
};


/// A single factor or set of factors whose parameters should be estimated.
/** To ensure that parameters can be shared between different factors during
 *  EM learning, each factor's values are reordered to match a desired variable 
 *  ordering. The ordering of the variables in a factor may therefore differ 
 *  from the canonical ordering used in libDAI. The SharedParameters
 *  class couples one or more factors (together with the specified orderings
 *  of the variables) with a ParameterEstimation object, taking care of the
 *  necessary permutations of the factor entries / parameters.
 */
class SharedParameters {
    public:
        /// Convenience label for an index into a factor in a FactorGraph.
        typedef size_t FactorIndex;
        /// Convenience label for a grouping of factor orientations.
        typedef std::map<FactorIndex, std::vector<Var> > FactorOrientations;

    private:
        /// Maps factor indices to the corresponding VarSets
        std::map<FactorIndex, VarSet> _varsets;
        /// Maps factor indices to the corresponding Permute objects that permute the desired ordering into the canonical ordering
        std::map<FactorIndex, Permute> _perms;
        /// Maps factor indices to the corresponding desired variable orderings
        FactorOrientations _varorders;
        /// Parameter estimation method to be used
        ParameterEstimation *_estimation;
        /// Indicates whether the object pointed to by _estimation should be deleted upon destruction
        bool _deleteEstimation;

        /// Calculates the permutation that permutes the variables in varorder into the canonical ordering
        static Permute calculatePermutation( const std::vector<Var> &varorder, VarSet &outVS );

        /// Initializes _varsets and _perms from _varorders
        void setPermsAndVarSetsFromVarOrders();

    public:
        /// Copy constructor
        SharedParameters( const SharedParameters &sp );

        /// Constructor 
        /** \param varorders  all the factor orientations for this parameter
         *  \param estimation a pointer to the parameter estimation method
         */ 
        SharedParameters( const FactorOrientations &varorders, ParameterEstimation *estimation );

        /// Constructor for making an object from a stream and a factor graph
        SharedParameters( std::istream &is, const FactorGraph &fg_varlookup );

        /// Destructor
        ~SharedParameters() { 
            if( _deleteEstimation )
                delete _estimation;
        }

        /// Collect the necessary statistics from expected values
        void collectSufficientStatistics( InfAlg &alg );

        /// Estimate and set the shared parameters
        void setParameters( FactorGraph &fg );
};


/// A MaximizationStep groups together several parameter estimation tasks into a single unit.
class MaximizationStep { 
    private:
        std::vector<SharedParameters> _params;

    public:
        /// Default constructor
        MaximizationStep() : _params() {}

        /// Constructor from a vector of SharedParameters objects
        MaximizationStep( std::vector<SharedParameters> &maximizations ) : _params(maximizations) {}  

        /// Constructor from an input stream and a corresponding factor graph
        MaximizationStep( std::istream &is, const FactorGraph &fg_varlookup );

        /// Collect the beliefs from this InfAlg as expectations for the next Maximization step.
        void addExpectations( InfAlg &alg );

        /// Using all of the currently added expectations, make new factors with maximized parameters and set them in the FactorGraph.
        void maximize( FactorGraph &fg );
};


/// EMAlg performs Expectation Maximization to learn factor parameters.
/** This requires specifying:
 *     - Evidence (instances of observations from the graphical model),
 *     - InfAlg for performing the E-step, which includes the factor graph,
 *     - a vector of MaximizationSteps steps to be performed.
 *
 *  This implementation can perform incremental EM by using multiple 
 *  MaximizationSteps.  An expectation step is performed between execution
 *  of each MaximizationStep.  A call to iterate() will cycle through all
 *  MaximizationSteps.
 *
 *  Having multiple and separate maximization steps allows for maximizing some
 *  parameters, performing another E step, and then maximizing separate
 *  parameters, which may result in faster convergence in some cases.
 */  
class EMAlg {
    private:
        /// All the data samples used during learning
        const Evidence &_evidence;

        /// How to do the expectation step
        InfAlg &_estep;

        /// The maximization steps to take
        std::vector<MaximizationStep> _msteps;

        /// Number of iterations done
        size_t _iters;

        /// History of likelihoods
        std::vector<Real> _lastLogZ;

        /// Maximum number of iterations
        size_t _max_iters;

        /// Convergence tolerance
        Real _log_z_tol;

    public:
        /// Key for setting maximum iterations @see setTermConditions
        static const std::string MAX_ITERS_KEY;
        /// Default maximum iterations @see setTermConditions
        static const size_t MAX_ITERS_DEFAULT;
        /// Key for setting likelihood termination condition @see setTermConditions
        static const std::string LOG_Z_TOL_KEY;
        /// Default likelihood tolerance @see setTermConditions
        static const Real LOG_Z_TOL_DEFAULT;

        /// Construct an EMAlg from all these objects
        EMAlg( const Evidence &evidence, InfAlg &estep, std::vector<MaximizationStep> &msteps, const PropertySet &termconditions ) 
          : _evidence(evidence), _estep(estep), _msteps(msteps), _iters(0), _lastLogZ(), _max_iters(MAX_ITERS_DEFAULT), _log_z_tol(LOG_Z_TOL_DEFAULT)
        { 
              setTermConditions( termconditions );
        }
  
        /// Construct an EMAlg from an Evidence object, an InfAlg object, and an input stream
        EMAlg( const Evidence &evidence, InfAlg &estep, std::istream &mstep_file );

        /// Change the conditions for termination
        /** There are two possible parameters in the PropertySet
         *    - max_iters maximum number of iterations
         *    - log_z_tol proportion of increase in logZ
         *
         *  \see hasSatisifiedTermConditions()
         */
        void setTermConditions( const PropertySet &p );

        /// Determine if the termination conditions have been met.
        /** There are two sufficient termination conditions:
         *    -# the maximum number of iterations has been performed
         *    -# the ratio of logZ increase over previous logZ is less than the 
         *       tolerance, i.e.,
         *       \f$ \frac{\log(Z_t) - \log(Z_{t-1})}{| \log(Z_{t-1}) | } < \mathrm{tol} \f$.
         */
        bool hasSatisfiedTermConditions() const;

        /// Returns number of iterations done so far
        size_t getCurrentIters() const { return _iters; }

        /// Perform an iteration over all maximization steps
        Real iterate();

        /// Perform an iteration over a single MaximizationStep
        Real iterate( MaximizationStep &mstep );

        /// Iterate until termination conditions are satisfied
        void run();
};


} // end of namespace dai


#endif
