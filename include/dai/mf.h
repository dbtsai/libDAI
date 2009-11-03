/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Defines class MF which implements the Mean Field algorithm


#ifndef __defined_libdai_mf_h
#define __defined_libdai_mf_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>


namespace dai {


/// Approximate inference algorithm "Mean Field"
/** The Mean Field algorithm iteratively calculates approximations of
 *  single variable marginals (beliefs). The update equation for 
 *  a single belief \f$b_i\f$ is given by:
 *    \f[ b_i^{\mathrm{new}}(x_i) \propto \prod_{I\in N_i} \exp \left( \sum_{x_{N_I \setminus \{i\}}} \log f_I(x_I) \prod_{j \in N_I \setminus \{i\}} b_j(x_j) \right) \f]
 *  These update equations are performed for all variables until convergence.
 */
class MF : public DAIAlgFG {
    private:
        /// Current approximations of single variable marginals
        std::vector<Factor>  _beliefs;
        /// Maximum difference encountered so far
        Real _maxdiff;
        /// Number of iterations needed
        size_t _iters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            Real tol;

            /// Damping constant
            Real damping;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
    /// \name Constructors/destructors
    //@{
        /// Default constructor
        MF() : DAIAlgFG(), _beliefs(), _maxdiff(0.0), _iters(0U), props() {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        MF( const FactorGraph &fg, const PropertySet &opts ) : DAIAlgFG(fg), _beliefs(), _maxdiff(0.0), _iters(0U), props() {
            setProperties( opts );
            construct();
        }
    //@}

    /// \name General InfAlg interface
    //@{
        virtual MF* clone() const { return new MF(*this); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    //@}

    /// \name Managing parameters (which are stored in MF::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding MF::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        /// Helper function for constructors
        void construct();

        /// Calculates an updated belief of variable \a i
        Factor calcNewBelief( size_t i );
};


} // end of namespace dai


#endif
