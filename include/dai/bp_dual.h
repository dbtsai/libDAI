/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2009  Frederik Eaton [frederik at ofb dot net]
 */


/// \file
/// \brief Defines class BP_dual
/// \todo This replicates a large part of the functionality of BP; would it not be shorter to adapt BP instead?
/// \todo Improve documentation


#ifndef __defined_libdai_bp_dual_h
#define __defined_libdai_bp_dual_h


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>


namespace dai {


/// Calculates both types of BP messages and their normalizers from an InfAlg.
/** BP_dual calculates "dual" versions of BP messages (both messages from factors
 *  to variables and messages from variables to factors), and normalizers, given an InfAlg.
 *  These are computed from the variable and factor beliefs of the InfAlg.
 *  This class is used primarily by BBP.
 */
class BP_dual {

    protected:
        /// Convenience label for storing edge properties
        template<class T>
        struct _edges_t : public std::vector<std::vector<T> > {};

        /// Groups together the data structures for storing the two types of messages and their normalizers
        struct messages {
            /// Unnormalized variable->factor messages
            _edges_t<Prob> n;
            /// Normalizers of variable->factor messages
            _edges_t<Real> Zn;
            /// Unnormalized Factor->variable messages
            _edges_t<Prob> m;
            /// Normalizers of factor->variable messages
            _edges_t<Real> Zm;
        };
        /// Stores all messages
        messages _msgs;

        /// Groups together the data structures for storing the two types of beliefs and their normalizers
        struct beliefs {
            /// Unnormalized variable beliefs
            std::vector<Prob> b1;
            /// Normalizers of variable beliefs
            std::vector<Real> Zb1;
            /// Unnormalized factor beliefs
            std::vector<Prob> b2;
            /// Normalizers of factor beliefs
            std::vector<Real> Zb2;
        };
        /// Stores all beliefs
        beliefs _beliefs;

        /// Pointer to the InfAlg object
        const InfAlg *_ia;

        /// Does all necessary preprocessing
        void init();
        /// Allocates space for _msgs
        void regenerateMessages();
        /// Allocates space for _beliefs
        void regenerateBeliefs();

        /// Calculate all messages from InfAlg beliefs
        void calcMessages();
        /// Update factor->variable message (i->I)
        void calcNewM(size_t i, size_t _I);
        /// Update variable->factor message (I->i)
        void calcNewN(size_t i, size_t _I);

        /// Calculate all variable and factor beliefs from messages
        void calcBeliefs();
        /// Calculate variable belief
        void calcBeliefV(size_t i);
        /// Calculate factor belief
        void calcBeliefF(size_t I);

    public:
        /// Construct BP_dual object from (converged) InfAlg object's beliefs and factors.
        /*  A pointer to the the InfAlg object is stored,
         *  so the object must not be destroyed before the BP_dual is destroyed.
         */
        BP_dual( const InfAlg *ia ) : _ia(ia) { init(); }

        /// Returns the underlying FactorGraph
        const FactorGraph& fg() const { return _ia->fg(); }

        /// Returns reference to factor -> var message (I->i)
        Prob & msgM( size_t i, size_t _I ) { return _msgs.m[i][_I]; }
        /// Returns constant reference to factor -> var message (I->i)
        const Prob & msgM( size_t i, size_t _I ) const { return _msgs.m[i][_I]; }
        /// Returns reference to var -> factor message (i->I)
        Prob & msgN( size_t i, size_t _I ) { return _msgs.n[i][_I]; }
        /// Returns constant reference to var -> factor message (i->I)
        const Prob & msgN( size_t i, size_t _I ) const { return _msgs.n[i][_I]; }
        /// Returns reference to normalizer for msgM
        Real & zM( size_t i, size_t _I ) { return _msgs.Zm[i][_I]; }
        /// Returns constant reference to normalizer for msgM
        const Real & zM( size_t i, size_t _I ) const { return _msgs.Zm[i][_I]; }
        /// Returns reference to normalizer for msgN
        Real & zN( size_t i, size_t _I ) { return _msgs.Zn[i][_I]; }
        /// Returns constant reference to normalizer for msgN
        const Real & zN( size_t i, size_t _I ) const { return _msgs.Zn[i][_I]; }

        /// Returns variable belief
        Factor beliefV( size_t i ) const { return Factor( _ia->fg().var(i), _beliefs.b1[i] ); }
        /// Returns factor belief
        Factor beliefF( size_t I ) const { return Factor( _ia->fg().factor(I).vars(), _beliefs.b2[I] ); }

        /// Returns normalizer for variable belief
        Real beliefVZ( size_t i ) const { return _beliefs.Zb1[i]; }
        /// Returns normalizer for factor belief
        Real beliefFZ( size_t I ) const { return _beliefs.Zb2[I]; }
};


} // end of namespace dai


#endif
