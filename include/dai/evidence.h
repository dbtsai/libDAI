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

#ifndef __defined_libdai_evidence_h
#define __defined_libdai_evidence_h

#include <istream>

#include <dai/daialg.h>

namespace dai {

/// Store joint observations on a graphical model.
class SampleData {
private:
  std::string _name;
  std::map<Var, size_t> _obs;
public:
  /// Construct an empty object
  SampleData() : _name(), _obs() {}
  /// Set the name of the sample
  void name(const std::string& name) { _name = name; }
  /// Get the name of the sample
  const std::string& name() const { return _name; }
  /// Read from the observation map
  const std::map<Var, size_t>& observations() const { return _obs; }
  /// Add an observation
  void addObservation(Var node, size_t setting);
  /// Add evidence by clamping variables to observed values.
  void applyEvidence(InfAlg& alg) const;
};

/// Store observations from a graphical model.
class Evidence {
private:
  std::map< std::string, SampleData > _samples;
public:
  /// Start with empty obects, then fill with calls to addEvidenceTabFile()
 Evidence() : _samples() {}
  
  /** Read in tab-data from a stream. Each line contains one sample, and
   * the first line is a header line with names. The first column contains
   * names for each of the samples.
   */
  void addEvidenceTabFile(std::istream& is,
			  std::map< std::string, Var >& varMap);

  /** Read in tab-data from a stream. Each line contains one sample,
   * and the first line is a header line with variable IDs. The first
   * column contains names for each of the samples.
   */
  void addEvidenceTabFile(std::istream& is, FactorGraph& fg);
  
  /// Total number of samples in this evidence file
  size_t nrSamples() const { return _samples.size(); }

  /// @name iterator interface
  //@{
  typedef std::map< std::string, SampleData >::iterator iterator;
  typedef std::map< std::string, SampleData >::const_iterator const_iterator;
  iterator begin() { return _samples.begin(); }
  const_iterator begin() const { return _samples.begin(); }
  iterator end() { return _samples.end(); }
  const_iterator end() const { return _samples.end(); }
  //@}

};
  
}

#endif
