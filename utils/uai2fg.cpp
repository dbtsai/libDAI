/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2008-2010  Joris Mooij  [joris dot mooij at libdai dot org]
 */


#include <iostream>
#include <fstream>
#include <dai/alldai.h>
#include <dai/util.h>
#include <dai/index.h>
#include <dai/jtree.h>


using namespace std;
using namespace dai;


/// Reads "evidence" (a mapping from observed variable labels to the observed values) from a UAI evidence file
map<size_t, size_t> ReadUAIEvidenceFile( char* filename ) {
    map<size_t, size_t> evid;

    // open file
    ifstream is;
    is.open( filename );
    if( is.is_open() ) {
        // read number of observed variables
        size_t nr_evid;
        is >> nr_evid;
        if( is.fail() )
            DAI_THROWE(INVALID_EVIDENCE_FILE,"Cannot read number of observed variables");

        // for each observation, read the variable label and the observed value
        for( size_t i = 0; i < nr_evid; i++ ) {
            size_t label, val;
            is >> label;
            if( is.fail() )
                DAI_THROWE(INVALID_EVIDENCE_FILE,"Cannot read label for " + toString(i) + "'th observed variable");
            is >> val;
            if( is.fail() )
                DAI_THROWE(INVALID_EVIDENCE_FILE,"Cannot read value of " + toString(i) + "'th observed variable");
            evid[label] = val;
        }

        // close file
        is.close();
    } else
        DAI_THROWE(CANNOT_READ_FILE,"Cannot read from file " + std::string(filename));

    return evid;
}


/// Reads factor graph (as a pair of a variable vector and factor vector) from a UAI factor graph file
void ReadUAIFGFile( const char *filename, size_t verbose, vector<Var>& vars, vector<Factor>& factors, vector<Permute>& permutations ) {
    vars.clear();
    factors.clear();
    permutations.clear();

    // open file
    ifstream is;
    is.open( filename );
    if( is.is_open() ) {
        size_t nrFacs, nrVars;
        string line;
        
        // read header line
        getline(is,line);
        if( is.fail() || (line != "BAYES" && line != "MARKOV" && line != "BAYES\r" && line != "MARKOV\r") )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"UAI factor graph file should start with \"BAYES\" or \"MARKOV\"");
        if( verbose >= 2 )
            cout << "Reading " << line << " network..." << endl;

        // read number of variables
        is >> nrVars;
        if( is.fail() )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of variables");
        if( verbose >= 2 )
            cout << "Reading " << nrVars << " variables..." << endl;

        // for each variable, read its number of states
        vars.reserve( nrVars );
        for( size_t i = 0; i < nrVars; i++ ) {
            size_t dim;
            is >> dim;
            if( is.fail() )
                DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of states of " + toString(i) + "'th variable");
            vars.push_back( Var( i, dim ) );
        }

        // read number of factors
        is >> nrFacs;
        if( is.fail() )
            DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of factors");
        if( verbose >= 2 )
            cout << "Reading " << nrFacs << " factors..." << endl;

        // for each factor, read the variables on which it depends
        vector<vector<Var> > factorVars;
        factors.reserve( nrFacs );
        factorVars.reserve( nrFacs );
        for( size_t I = 0; I < nrFacs; I++ ) {
            if( verbose >= 3 )
                cout << "Reading factor " << I << "..." << endl;

            // read number of variables for factor I
            size_t I_nrVars;
            is >> I_nrVars;
            if( is.fail() )
                DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of variables for " + toString(I) + "'th factor");
            if( verbose >= 3 )
                cout << "  which depends on " << I_nrVars << " variables" << endl;

            // read the variable labels
            vector<long> I_labels;
            vector<size_t> I_dims;
            I_labels.reserve( I_nrVars );
            I_dims.reserve( I_nrVars );
            factorVars[I].reserve( I_nrVars );
            for( size_t _i = 0; _i < I_nrVars; _i++ ) {
                long label;
                is >> label;
                if( is.fail() )
                    DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read variable labels for " + toString(I) + "'th factor");
                I_labels.push_back( label );
                I_dims.push_back( vars[label].states() );
                factorVars[I].push_back( vars[label] );
            }
            if( verbose >= 3 )
                cout << "  labels: " << I_labels << ", dimensions " << I_dims << endl;

            // add the factor and the labels
            factors.push_back( Factor( VarSet( factorVars[I].begin(), factorVars[I].end(), factorVars[I].size() ), (Real)0 ) );
        }

        // for each factor, read its values
        permutations.reserve( nrFacs );
        for( size_t I = 0; I < nrFacs; I++ ) {
            if( verbose >= 3 )
                cout << "Reading factor " << I << "..." << endl;

            // calculate permutation object, reversing the indexing in factorVars[I] first
            Permute permindex( factorVars[I], true );
            permutations.push_back( permindex );

            // read factor values
            size_t nrNonZeros;
            is >> nrNonZeros;
            if( is.fail() )
                DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read number of nonzero factor values for " + toString(I) + "'th factor");
            if( verbose >= 3 ) 
                cout << "  number of nonzero values: " << nrNonZeros << endl;
            DAI_ASSERT( nrNonZeros == factors[I].nrStates() );
            for( size_t li = 0; li < nrNonZeros; li++ ) {
                Real val;
                is >> val;
                if( is.fail() )
                    DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Cannot read factor values of " + toString(I) + "'th factor");
                // assign value after calculating its linear index corresponding to the permutation
                if( verbose >= 4 )
                    cout << "  " << li << "'th value " << val << " corresponds with index " << permindex.convertLinearIndex(li) << endl;
                factors[I].set( permindex.convertLinearIndex( li ), val );
            }
        }
        if( verbose >= 3 )
            cout << "variables:" << vars << endl;
        if( verbose >= 3 )
            cout << "factors:" << factors << endl;

        // close file
        is.close();
    } else
        DAI_THROWE(CANNOT_READ_FILE,"Cannot read from file " + std::string(filename));
}


int main( int argc, char *argv[] ) {
    if ( argc != 7 ) {
        cout << "This program is part of libDAI - http://www.libdai.org/" << endl << endl;
        cout << "Usage: ./uai2fg <filename.uai> <filename.uai.evid> <filename.fg> <type> <run_jtree> <verbose>" << endl << endl;
        cout << "Converts input files in the UAI 2008 approximate inference evaluation format" << endl;
        cout << "(see http://graphmod.ics.uci.edu/uai08/) to the libDAI factor graph format." << endl;
        cout << "Reads factor graph <filename.uai> and evidence <filename.uai.evid>" << endl;
        cout << "and writes the resulting clamped factor graph to <filename.fg>." << endl;
        cout << "If type==0, uses surgery (recommended), otherwise, uses just adds delta factors." << endl;
        cout << "If run_jtree!=0, runs a junction tree and reports the results in the UAI 2008 results file format." << endl;
        return 1;
    } else {
        long verbose = atoi( argv[6] );
        long type = atoi( argv[4] );
        bool run_jtree = atoi( argv[5] );

        // read factor graph
        vector<Var> vars;
        vector<Factor> facs;
        vector<Permute> permutations;
        ReadUAIFGFile( argv[1], verbose, vars, facs, permutations );

        // read evidence
        map<size_t,size_t> evid = ReadUAIEvidenceFile( argv[2] );

        // construct unclamped factor graph
        FactorGraph fg0( facs.begin(), facs.end(), vars.begin(), vars.end(), facs.size(), vars.size() );

        // change factor graph to reflect observed evidence
        if( type == 0 ) {
            // replace factors with clamped variables with slices
            for( size_t I = 0; I < facs.size(); I++ ) {
                for( map<size_t,size_t>::const_iterator e = evid.begin(); e != evid.end(); e++ ) {
                    if( facs[I].vars() >> vars[e->first] ) {
                        if( verbose >= 2 )
                            cout << "Clamping " << e->first << " to value " << e->second << " in factor " << I << " = " << facs[I].vars() << endl;
                        facs[I] = facs[I].slice( vars[e->first], e->second );
                        if( verbose >= 2 )
                            cout << "...remaining vars: " << facs[I].vars() << endl;
                    }
                }
            }
            // remove empty factors
            double logZcorr = 0.0;
            for( vector<Factor>::iterator I = facs.begin(); I != facs.end(); )
                if( I->vars().size() == 0 ) {
                    logZcorr += std::log( (Real)(*I)[0] );
                    I = facs.erase( I );
                } else
                    I++;
            // multiply with logZcorr constant
            if( facs.size() == 0 )
                facs.push_back( Factor( VarSet(), std::exp(logZcorr) ) );
            else
                facs.front() *= std::exp(logZcorr);
        }
        // add delta factors corresponding to observed variable values
        for( map<size_t,size_t>::const_iterator e = evid.begin(); e != evid.end(); e++ )
            facs.push_back( createFactorDelta( vars[e->first], e->second ) );

        // construct clamped factor graph
        FactorGraph fg( facs.begin(), facs.end(), vars.begin(), vars.end(), facs.size(), vars.size() );

        // write it to a file
        fg.WriteToFile( argv[3] );

        // if requested, perform various inference tasks
        if( run_jtree ) {
            // construct junction tree on unclamped factor graph
            JTree jt0( fg0, PropertySet()("updates",string("HUGIN")) );
            jt0.init();
            jt0.run();

            // construct junction tree on clamped factor graph
            JTree jt( fg, PropertySet()("updates",string("HUGIN")) );
            jt.init();
            jt.run();

            // output probability of evidence
            cout.precision( 8 );
            if( evid.size() )
                cout << "z " << (jt.logZ() - jt0.logZ()) / dai::log((Real)10.0) << endl;
            else
                cout << "z " << jt.logZ() / dai::log((Real)10.0) << endl;

            // output variable marginals
            cout << "m " << jt.nrVars() << " ";
            for( size_t i = 0; i < jt.nrVars(); i++ ) {
                cout << jt.var(i).states() << " ";
                for( size_t s = 0; s < jt.var(i).states(); s++ )
                    cout << jt.beliefV(i)[s] << " ";
            }
            cout << endl;

            // calculate MAP state
            jt.props.inference = JTree::Properties::InfType::MAXPROD;
            jt.init();
            jt.run();
            vector<size_t> MAP = jt.findMaximum();
            map<Var, size_t> state;
            for( size_t i = 0; i < MAP.size(); i++ )
                state[jt.var(i)] = MAP[i];
            double log_MAP_prob = 0.0;
            for( size_t I = 0; I < jt.nrFactors(); I++ )
                log_MAP_prob += dai::log( jt.factor(I)[calcLinearState( jt.factor(I).vars(), state )] );

            // output MAP state
            cout << "s ";
            if( evid.size() )
                cout << (log_MAP_prob - jt0.logZ()) / dai::log((Real)10.0) << " ";
            else
                cout << log_MAP_prob / dai::log((Real)10.0) << " ";
            cout << jt.nrVars() << " ";
            for( size_t i = 0; i < jt.nrVars(); i++ )
                cout << MAP[i] << " ";
            cout << endl;
        }
    }

    return 0;
}
