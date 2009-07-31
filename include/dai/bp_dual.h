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


/// \file
/// \brief Defines class BP_dual
/// \todo Improve documentation
/// \todo Clean up


#ifndef __defined_libdai_bp_dual_h
#define __defined_libdai_bp_dual_h


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>


namespace dai {


/** Class to estimate "dual" versions of BP messages, and normalizers, given an InfAlg. 
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
            _edges_t<Prob> n;
            _edges_t<Real> Zn;
            _edges_t<Prob> m;
            _edges_t<Real> Zm;
        };
        messages _msgs;

        /// Groups together the data structures for storing the two types of beliefs and their normalizers
        struct beliefs {
            // indexed by node
            std::vector<Prob> b1;
            std::vector<Real> Zb1;
            // indexed by factor
            std::vector<Prob> b2;
            std::vector<Real> Zb2;
        };
        beliefs _beliefs;

        const InfAlg *_ia;
            
        void init();

        void regenerateMessages();
        void regenerateBeliefs();

        void calcMessages();
        void calcBeliefV(size_t i);
        void calcBeliefF(size_t I);
        void calcBeliefs();

        void calcNewM(size_t i, size_t _I);
        void calcNewN(size_t i, size_t _I);

    public:
        /// Construct BP_dual object from (converged) InfAlg object's beliefs and factors. 
        /*  A pointer to the the InfAlg object is stored, 
         *  so the object must not be destroyed before the BP_dual is destroyed.
         */
        BP_dual( const InfAlg *ia ) : _ia(ia) { init(); }

        const FactorGraph& fg() const { return _ia->fg(); }

        /// factor -> var message
        DAI_ACCMUT(Prob & msgM(size_t i, size_t _I), { return _msgs.m[i][_I]; });
        /// var -> factor message
        DAI_ACCMUT(Prob & msgN(size_t i, size_t _I), { return _msgs.n[i][_I]; });
        /// Normalizer for msgM
        DAI_ACCMUT(Real & zM(size_t i, size_t _I), { return _msgs.Zm[i][_I]; });
        /// Normalizer for msgN
        DAI_ACCMUT(Real & zN(size_t i, size_t _I), { return _msgs.Zn[i][_I]; });

        /// Variable belief
        Factor beliefV(size_t i) const { return Factor(_ia->fg().var(i), _beliefs.b1[i]); }
        /// Factor belief
        Factor beliefF(size_t I) const { return Factor(_ia->fg().factor(I).vars(), _beliefs.b2[I]); }

        /// Normalizer for variable belief
        Real beliefVZ(size_t i) const { return _beliefs.Zb1[i]; }
        /// Normalizer for factor belief
        Real beliefFZ(size_t I) const { return _beliefs.Zb2[I]; }

};


} // end of namespace dai


#endif
