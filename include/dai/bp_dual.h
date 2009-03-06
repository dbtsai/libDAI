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


#ifndef ____defined_libdai_bp_dual_h__
#define ____defined_libdai_bp_dual_h__


#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/enum.h>
#include <dai/bp.h>


namespace dai {

    
struct BP_dual_messages {
    // messages:
    // indexed by edge index (using VV2E)
    std::vector<Prob> n;
    std::vector<Real> Zn;
    std::vector<Prob> m;
    std::vector<Real> Zm;
};


struct BP_dual_beliefs {
    // beliefs:
    // indexed by node
    std::vector<Prob> b1;
    std::vector<Real> Zb1;
    // indexed by factor
    std::vector<Prob> b2;
    std::vector<Real> Zb2;
};


void _clamp( FactorGraph &g, const Var &n, const std::vector<size_t> &is );


/// Clamp a factor to have one of a set of values
void _clampFactor( FactorGraph &g, size_t I, const std::vector<size_t> &is );


class BP_dual : public DAIAlgFG {
    public:
        typedef std::vector<size_t>  _ind_t;

    protected:
        // indexed by edge index. for each edge i->I, contains a
        // vector whose entries correspond to those of I, and the
        // value of each entry is the corresponding entry of i
        std::vector<_ind_t>          _indices; 

        BP_dual_messages _msgs;
        BP_dual_messages _new_msgs;
    public:
        BP_dual_beliefs _beliefs;

        size_t _iters;
        double _maxdiff;

        struct Properties {
          typedef BP::Properties::UpdateType UpdateType;
          UpdateType updates;
          double tol;
          size_t maxiter;
          size_t verbose;
        } props;

        /// List of property names
        static const char *PropertyList[];
        /// Name of this inference algorithm
        static const char *Name;

    public:
        void Regenerate(); // used by constructor
        void RegenerateIndices();
        void RegenerateMessages();
        void RegenerateBeliefs();

        void CalcBelief1(size_t i);
        void CalcBelief2(size_t I);
        void CalcBeliefs(); // called after run()

        void calcNewM(size_t iI);
        void calcNewN(size_t iI);
        void upMsgM(size_t iI);
        void upMsgN(size_t iI);

        /* DAI_ENUM(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL) */
        typedef BP::Properties::UpdateType UpdateType;
        UpdateType Updates() const { return props.updates; }
        size_t Verbose() const { return props.verbose; }

        /// Default constructor
        BP_dual() {}

        /// construct BP_dual object from FactorGraph
        BP_dual(const FactorGraph & fg, const PropertySet &opts) : DAIAlgFG(fg) {
            setProperties(opts);
            Regenerate();
        }
        
        DAI_ACCMUT(Prob & msgM(size_t I, size_t i), { return _msgs.m[VV2E(i,I)]; });
        DAI_ACCMUT(Prob & msgN(size_t i, size_t I), { return _msgs.n[VV2E(i,I)]; });
        DAI_ACCMUT(Prob & msgM(size_t iI), { return _msgs.m[iI]; });
        DAI_ACCMUT(Prob & msgN(size_t iI), { return _msgs.n[iI]; });
        DAI_ACCMUT(Real & zM(size_t I, size_t i), { return _msgs.Zm[VV2E(i,I)]; });
        DAI_ACCMUT(Real & zN(size_t i, size_t I), { return _msgs.Zn[VV2E(i,I)]; });
        DAI_ACCMUT(Real & zM(size_t iI), { return _msgs.Zm[iI]; });
        DAI_ACCMUT(Real & zN(size_t iI), { return _msgs.Zn[iI]; });
        DAI_ACCMUT(Prob & newMsgM(size_t I, size_t i), { return _new_msgs.m[VV2E(i,I)]; });
        DAI_ACCMUT(Prob & newMsgN(size_t i, size_t I), { return _new_msgs.n[VV2E(i,I)]; });
        DAI_ACCMUT(Real & newZM(size_t I, size_t i), { return _new_msgs.Zm[VV2E(i,I)]; });
        DAI_ACCMUT(Real & newZN(size_t i, size_t I), { return _new_msgs.Zn[VV2E(i,I)]; });

        DAI_ACCMUT(_ind_t & index(size_t i, size_t I), { return( _indices[VV2E(i,I)] ); });

        Real belief1Z(size_t i) const { return _beliefs.Zb1[i]; }
        Real belief2Z(size_t I) const { return _beliefs.Zb2[I]; }

        size_t doneIters() const { return _iters; }


        /// @name General InfAlg interface
        //@{
        virtual BP_dual* clone() const { return new BP_dual(*this); }
        virtual BP_dual* create() const { return new BP_dual(); }
//        virtual BP_dual* create() const { return NULL; }
        virtual std::string identify() const;
        virtual Factor belief (const Var &n) const { return( belief1( findVar( n ) ) ); }
        virtual Factor belief (const VarSet &n) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual double run();
        virtual double maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        //@}

        void init(const std::vector<size_t>& state);
        Factor belief1 (size_t i) const { return Factor(var(i), _beliefs.b1[i]); }
        Factor belief2 (size_t I) const { return Factor(factor(I).vars(), _beliefs.b2[I]); }

        void updateMaxDiff( double maxdiff ) { if( maxdiff > _maxdiff ) _maxdiff = maxdiff; }

        /// Set Props according to the PropertySet opts, where the values can be stored as std::strings or as the type of the corresponding Props member
        void setProperties( const PropertySet &opts );
        PropertySet getProperties() const;
        std::string printProperties() const;
};


} // end of namespace dai


#endif
