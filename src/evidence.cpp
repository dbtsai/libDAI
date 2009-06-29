#include <sstream>
#include <string>
#include <cstdlib>

#include <dai/util.h>

#include <dai/evidence.h>

namespace dai {

void SampleData::addObservation(Var node, size_t setting) {
    _obs[node] = setting;
}
  
void SampleData::applyEvidence(InfAlg& alg) const {
  std::map< Var, size_t>::const_iterator i = _obs.begin();
  for( ; i != _obs.end(); ++i) {
    alg.clamp(i->first, i->second);
  }
}
  
void Evidence::addEvidenceTabFile(std::istream& is, FactorGraph& fg) {
  std::map< std::string, Var > varMap;
  std::vector< Var >::const_iterator v = fg.vars().begin();
  for(; v != fg.vars().end(); ++v) {
    std::stringstream s;
    s << v->label();
    varMap[s.str()] = *v;
  }

  addEvidenceTabFile(is, varMap);
}

void Evidence::addEvidenceTabFile(std::istream& is, 
				  std::map< std::string, Var >& varMap) {
  
  std::vector< std::string > header_fields;
  std::vector< Var > vars;
  std::string line;
  getline(is, line);
  
  // Parse header
  tokenizeString(line, header_fields);
  std::vector< std::string >::const_iterator p_field = header_fields.begin();

  if (p_field == header_fields.end()) { DAI_THROW(INVALID_EVIDENCE_LINE); }

  ++p_field; // first column are sample labels
  for ( ; p_field != header_fields.end(); ++p_field) {
    std::map< std::string, Var >::iterator elem = varMap.find(*p_field);
    if (elem == varMap.end()) {
      DAI_THROW(INVALID_EVIDENCE_FILE);
    }
    vars.push_back(elem->second);
  }
  
  
  // Read samples
  while(getline(is, line)) {
    std::vector< std::string > fields;

    tokenizeString(line, fields);
    
    if (fields.size() != vars.size() + 1) { DAI_THROW(INVALID_EVIDENCE_LINE); }
    
    SampleData& sampleData = _samples[fields[0]];
    sampleData.name(fields[0]); // in case of a new sample
    for (size_t i = 0; i < vars.size(); ++i) {
      if (fields[i+1].size() > 0) { // skip if missing observation
	if (fields[i+1].find_first_not_of("0123456789") != std::string::npos) {
	  DAI_THROW(INVALID_EVIDENCE_OBSERVATION);
	}
	size_t state = atoi(fields[i+1].c_str());
	if (state >= vars[i].states()) {
	  DAI_THROW(INVALID_EVIDENCE_OBSERVATION);
	}
	sampleData.addObservation(vars[i], state);
      }
    }
  } // finished sample line
}

}
