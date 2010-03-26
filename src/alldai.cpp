/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <string>
#include <dai/alldai.h>
#include <dai/properties.h>
#include <dai/exceptions.h>


namespace dai {


using namespace std;


InfAlg *newInfAlg( const std::string &name, const FactorGraph &fg, const PropertySet &opts ) {
    if( name == ExactInf::Name )
        return new ExactInf (fg, opts);
#ifdef DAI_WITH_BP
    if( name == BP::Name )
        return new BP (fg, opts);
#endif
#ifdef DAI_WITH_FBP
    if( name == FBP::Name )
        return new FBP (fg, opts);
#endif
#ifdef DAI_WITH_TRWBP
    if( name == TRWBP::Name )
        return new TRWBP (fg, opts);
#endif
#ifdef DAI_WITH_MF
    if( name == MF::Name )
        return new MF (fg, opts);
#endif
#ifdef DAI_WITH_HAK
    if( name == HAK::Name )
        return new HAK (fg, opts);
#endif
#ifdef DAI_WITH_LC
    if( name == LC::Name )
        return new LC (fg, opts);
#endif
#ifdef DAI_WITH_TREEEP
    if( name == TreeEP::Name )
        return new TreeEP (fg, opts);
#endif
#ifdef DAI_WITH_JTREE
    if( name == JTree::Name )
        return new JTree (fg, opts);
#endif
#ifdef DAI_WITH_MR
    if( name == MR::Name )
        return new MR (fg, opts);
#endif
#ifdef DAI_WITH_GIBBS
    if( name == Gibbs::Name )
        return new Gibbs (fg, opts);
#endif
#ifdef DAI_WITH_CBP
    if( name == CBP::Name )
        return new CBP (fg, opts);
#endif
    DAI_THROWE(UNKNOWN_DAI_ALGORITHM,"Unknown libDAI algorithm: " + name);
}


InfAlg *newInfAlgFromString( const std::string &nameOpts, const FactorGraph &fg ) {
    pair<string,PropertySet> no = parseNameProperties( nameOpts );
    return newInfAlg( no.first, fg, no.second );
}


InfAlg *newInfAlgFromString( const std::string &nameOpts, const FactorGraph &fg, const std::map<std::string,std::string> &aliases ) {
    pair<string,PropertySet> no = parseNameProperties( nameOpts, aliases );
    return newInfAlg( no.first, fg, no.second );
}


std::pair<std::string, PropertySet> parseNameProperties( const std::string &s ) {
    string::size_type pos = s.find_first_of('[');
    string name;
    PropertySet opts;
    if( pos == string::npos ) {
        name = s;
    } else {
        name = s.substr(0,pos);

        stringstream ss;
        ss << s.substr(pos,s.length());
        ss >> opts;
    }
    return make_pair(name,opts);
}


std::pair<std::string, PropertySet> parseNameProperties( const std::string &s, const std::map<std::string,std::string> &aliases ) {
    // break string into method[properties]
    pair<string,PropertySet> ps = parseNameProperties(s);
    bool looped = false;

    // as long as 'method' is an alias, update:
    while( aliases.find(ps.first) != aliases.end() && !looped ) {
        string astr = aliases.find(ps.first)->second;
        pair<string,PropertySet> aps = parseNameProperties(astr);
        if( aps.first == ps.first )
            looped = true;
        // override aps properties by ps properties
        aps.second.set( ps.second );
        // replace ps by aps
        ps = aps;
        // repeat until method name == alias name ('looped'), or
        // there is no longer an alias 'method'
    }

    return ps;
}


std::map<std::string,std::string> readAliasesFile( const std::string &filename ) {
    // Read aliases
    map<string,string> result;
    ifstream infile;
    infile.open( filename.c_str() );
    if( infile.is_open() ) {
        while( true ) {
            string line;
            getline( infile,line );
            if( infile.fail() )
                break;
            if( (!line.empty()) && (line[0] != '#') ) {
                string::size_type pos = line.find(':',0);
                if( pos == string::npos )
                    DAI_THROWE(INVALID_ALIAS,"Invalid alias '" + line + "'");
                else {
                    string::size_type posl = line.substr(0, pos).find_last_not_of(" \t");
                    string key = line.substr(0, posl + 1);
                    string::size_type posr = line.substr(pos + 1, line.length()).find_first_not_of(" \t");
                    string val = line.substr(pos + 1 + posr, line.length());
                    result[key] = val;
                }
            }
        }
        infile.close();
    } else
        DAI_THROWE(CANNOT_READ_FILE,"Error opening aliases file " + filename);
    return result;
}


} // end of namespace dai
