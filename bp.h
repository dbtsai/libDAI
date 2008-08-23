/*  Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
    Radboud University Nijmegen, The Netherlands
    
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


#ifndef __BP_H__
#define __BP_H__


#include "daialg.h"
#include "factorgraph.h"
#include "enum.h"


using namespace std;


class BP : public DAIAlgFG {
    protected:
        typedef vector<size_t>  _ind_t;

        vector<_ind_t>          _indices;
        vector<Prob>            _messages, _newmessages;    

    public:
        ENUM4(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL)
        UpdateType Updates() const { return GetPropertyAs<UpdateType>("updates"); }

        // default constructor
        BP() : DAIAlgFG() {};
        // copy constructor
        BP(const BP & x) : DAIAlgFG(x), _indices(x._indices), _messages(x._messages), _newmessages(x._newmessages) {};
        BP* clone() const { return new BP(*this); }
        // construct BP object from FactorGraph
        BP(const FactorGraph & fg, const Properties &opts) : DAIAlgFG(fg, opts) {
            assert( checkProperties() );
            Regenerate();
        }
        // assignment operator
        BP & operator=(const BP & x) {
            if(this!=&x) {
                DAIAlgFG::operator=(x);
                _messages = x._messages;
                _newmessages = x._newmessages;
                _indices = x._indices;
            }
            return *this;
        }
        
        static const char *Name;

        Prob & message(size_t i1, size_t i2) { return( _messages[VV2E(i1,i2)] ); }  
        const Prob & message(size_t i1, size_t i2) const { return( _messages[VV2E(i1,i2)] ); }  
        Prob & newMessage(size_t i1, size_t i2) { return( _newmessages[VV2E(i1,i2)] ); }    
        const Prob & newMessage(size_t i1, size_t i2) const { return( _newmessages[VV2E(i1,i2)] ); }    
        _ind_t & index(size_t i1, size_t i2) { return( _indices[VV2E(i1,i2)] ); }
        const _ind_t & index(size_t i1, size_t i2) const { return( _indices[VV2E(i1,i2)] ); }

        string identify() const;
        void Regenerate();
        void init();
        void calcNewMessage(size_t iI);
        double run();
        Factor belief1 (size_t i) const;
        Factor belief2 (size_t I) const;
        Factor belief (const Var &n) const;
        Factor belief (const VarSet &n) const;
        vector<Factor> beliefs() const;
        Complex logZ() const;

        void init( const VarSet &ns );
        void undoProbs( const VarSet &ns ) { FactorGraph::undoProbs(ns); init(ns); }
        bool checkProperties();
};

#endif
