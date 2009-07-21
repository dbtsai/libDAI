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


/// \file
/// \brief Defines classes Evidence and Observation
/// \todo Improve documentation


#ifndef __defined_libdai_evidence_h
#define __defined_libdai_evidence_h


#include <istream>
#include <dai/daialg.h>


namespace dai {


/// Stores observed values of a subset of variables
class Observation {
    private:
        std::map<Var, size_t> _obs;

    public:
        /// Default constructor
        Observation() : _obs() {}
        /// Get all observations
        const std::map<Var, size_t>& observations() const { return _obs; }
        /// Add an observation
        void addObservation( Var node, size_t setting );
        /// Clamp variables to observed values
        void applyEvidence( InfAlg& alg ) const;
};


/// Stores multiple observations of variables
class Evidence {
    private:
        std::vector<Observation> _samples;

    public:
        /// Default constructor
        Evidence() : _samples() {}
      
        /** Read in tab-data from a stream. Each line contains one sample, and
         *  the first line is a header line with names.
         */
        void addEvidenceTabFile(std::istream& is, std::map<std::string, Var>& varMap);

        /** Read in tab-data from a stream. Each line contains one sample,
         *  and the first line is a header line with variable labels. 
         */
        void addEvidenceTabFile(std::istream& is, FactorGraph& fg);
      
        /// Returns total number of samples
        size_t nrSamples() const { return _samples.size(); }

        /// @name iterator interface
        //@{
        typedef std::vector<Observation>::iterator iterator;
        typedef std::vector<Observation>::const_iterator const_iterator;

        iterator begin() { return _samples.begin(); }
        const_iterator begin() const { return _samples.begin(); }
        iterator end() { return _samples.end(); }
        const_iterator end() const { return _samples.end(); }
        //@}
};


} // end of namespace dai


#endif
