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
/// \brief Defines class BP
/// \todo Improve documentation


#ifndef __defined_libdai_bp_h
#define __defined_libdai_bp_h


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>


namespace dai {


/// Approximate inference algorithm "(Loopy) Belief Propagation"
class BP : public DAIAlgFG {
    private:
        typedef std::vector<size_t> ind_t;
        typedef std::multimap<Real, std::pair<std::size_t, std::size_t> > LutType;
        struct EdgeProp {
            ind_t  index;
            Prob   message;
            Prob   newMessage;
            Real   residual;
        };
        std::vector<std::vector<EdgeProp> > _edges;
        std::vector<std::vector<LutType::iterator> > _edge2lut;
        LutType _lut;
        /// Maximum difference encountered so far
        Real _maxdiff;
        /// Number of iterations needed
        size_t _iters;
        /// The history of message updates (only recorded if recordSentMessages is true)
        std::vector<std::pair<std::size_t, std::size_t> > _sentMessages;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible update schedules
            DAI_ENUM(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL);

            /// Enumeration of inference variants
            DAI_ENUM(InfType,SUMPROD,MAXPROD);

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            Real tol;

            /// Do updates in logarithmic domain?
            bool logdomain;

            /// Damping constant
            Real damping;

            /// Update schedule
            UpdateType updates;

            /// Type of inference: sum-product or max-product?
            InfType inference;
        } props;

        /// Name of this inference algorithm
        static const char *Name;

        /// Specifies whether the history of message updates should be recorded
        bool recordSentMessages;

    public:
        /// Default constructor
        BP() : DAIAlgFG(), _edges(), _edge2lut(), _lut(), _maxdiff(0.0), _iters(0U), _sentMessages(), props(), recordSentMessages(false) {}

        /// Copy constructor
        BP( const BP &x ) : DAIAlgFG(x), _edges(x._edges), _edge2lut(x._edge2lut),
            _lut(x._lut), _maxdiff(x._maxdiff), _iters(x._iters), _sentMessages(x._sentMessages),
            props(x.props), recordSentMessages(x.recordSentMessages)
        {
            for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
                _edge2lut[l->second.first][l->second.second] = l;
        }

        /// Assignment operator
        BP& operator=( const BP &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _edges = x._edges;
                _lut = x._lut;
                for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
                    _edge2lut[l->second.first][l->second.second] = l;
                _maxdiff = x._maxdiff;
                _iters = x._iters;
                _sentMessages = x._sentMessages;
                props = x.props;
                recordSentMessages = x.recordSentMessages;
            }
            return *this;
        }

        /// Construct from FactorGraph fg and PropertySet opts
        BP( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg), _edges(), _maxdiff(0.0), _iters(0U), _sentMessages(), props(), recordSentMessages(false) {
            setProperties( opts );
            construct();
        }

    /// \name General InfAlg interface
    //@{
        virtual BP* clone() const { return new BP(*this); }
        virtual std::string identify() const;
        virtual Factor belief( const Var &n ) const;
        virtual Factor belief( const VarSet &ns ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual Factor beliefF( size_t I ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    //@}

    /// \name Additional interface specific for BP
    //@{
        /// Calculates the joint state of all variables that has maximum probability
        /** Assumes that run() has been called and that props.inference == MAXPROD
         */
        std::vector<std::size_t> findMaximum() const;

        /// Returns history of sent messages
        const std::vector<std::pair<std::size_t, std::size_t> >& getSentMessages() const {
            return _sentMessages;
        }

        /// Clears history of sent messages
        void clearSentMessages() {
            _sentMessages.clear();
        }
    //@}

    /// \name Managing parameters (which are stored in BP::props)
    //@{
        /// Set parameters of this inference algorithm.
        /** The parameters are set according to \a opts. 
         *  The values can be stored either as std::string or as the type of the corresponding BP::props member.
         */
        void setProperties( const PropertySet &opts );
        /// Returns parameters of this inference algorithm converted into a PropertySet.
        PropertySet getProperties() const;
        /// Returns parameters of this inference algorithm formatted as a string in the format "[key1=val1,key2=val2,...,keyn=valn]".
        std::string printProperties() const;
    //@}

    private:
        const Prob & message(size_t i, size_t _I) const { return _edges[i][_I].message; }
        Prob & message(size_t i, size_t _I) { return _edges[i][_I].message; }
        Prob & newMessage(size_t i, size_t _I) { return _edges[i][_I].newMessage; }
        const Prob & newMessage(size_t i, size_t _I) const { return _edges[i][_I].newMessage; }
        ind_t & index(size_t i, size_t _I) { return _edges[i][_I].index; }
        const ind_t & index(size_t i, size_t _I) const { return _edges[i][_I].index; }
        Real & residual(size_t i, size_t _I) { return _edges[i][_I].residual; }
        const Real & residual(size_t i, size_t _I) const { return _edges[i][_I].residual; }

        void calcNewMessage( size_t i, size_t _I );
        void updateMessage( size_t i, size_t _I );
        void updateResidual( size_t i, size_t _I, Real r );
        void findMaxResidual( size_t &i, size_t &_I );
        /// Calculates unnormalized belief of variable
        void calcBeliefV( size_t i, Prob &p ) const;
        /// Calculates unnormalized belief of factor
        void calcBeliefF( size_t I, Prob &p ) const;

        void construct();
};


} // end of namespace dai


#endif
