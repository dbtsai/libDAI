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


#include <sstream>
#include <string>
#include <cstdlib>

#include <dai/util.h>
#include <dai/evidence.h>


namespace dai {


void Observation::addObservation( Var node, size_t setting ) {
    _obs[node] = setting;
}


void Observation::applyEvidence( InfAlg &alg ) const {
    for( std::map<Var, size_t>::const_iterator i = _obs.begin(); i != _obs.end(); ++i )
        alg.clamp( alg.fg().findVar(i->first), i->second );
}
  

void Evidence::addEvidenceTabFile( std::istream &is, FactorGraph &fg ) {
    std::map<std::string, Var> varMap;
    for( std::vector<Var>::const_iterator v = fg.vars().begin(); v != fg.vars().end(); ++v ) {
        std::stringstream s;
        s << v->label();
        varMap[s.str()] = *v;
    }

    addEvidenceTabFile( is, varMap );
}


void Evidence::addEvidenceTabFile( std::istream &is, std::map<std::string, Var> &varMap ) {
    std::string line;
    getline( is, line );
    
    // Parse header
    std::vector<std::string> header_fields;
    tokenizeString( line, header_fields );
    std::vector<std::string>::const_iterator p_field = header_fields.begin();
    if( p_field == header_fields.end() ) 
        DAI_THROW(INVALID_EVIDENCE_FILE);

    std::vector<Var> vars;
    for( ; p_field != header_fields.end(); ++p_field ) {
        std::map<std::string, Var>::iterator elem = varMap.find( *p_field );
        if( elem == varMap.end() )
            DAI_THROW(INVALID_EVIDENCE_FILE);
        vars.push_back( elem->second );
    }
    
    // Read samples
    while( getline(is, line) ) {
        std::vector<std::string> fields;
        tokenizeString( line, fields );
        if( fields.size() != vars.size() ) 
            DAI_THROW(INVALID_EVIDENCE_FILE);
        
        Observation sampleData;
        for( size_t i = 0; i < vars.size(); ++i ) {
            if( fields[i].size() > 0 ) { // skip if missing observation
                if( fields[i].find_first_not_of("0123456789") != std::string::npos )
                    DAI_THROW(INVALID_EVIDENCE_FILE);
                size_t state = atoi( fields[i].c_str() );
                if( state >= vars[i].states() )
                    DAI_THROW(INVALID_EVIDENCE_FILE);
                sampleData.addObservation( vars[i], state );
            }
        }
        _samples.push_back( sampleData );
    } // finished sample line
}


} // end of namespace dai
