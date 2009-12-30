/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


#include <iostream>
#include <fstream>
#include <map>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <boost/program_options.hpp>
#include <dai/util.h>
#include <dai/alldai.h>


using namespace std;
using namespace dai;
namespace po = boost::program_options;


class TestDAI {
    protected:
        InfAlg          *obj;
        string          name;
        vector<Real>    err;

    public:
        vector<Factor>  varmargs;
        vector<Factor>  marginals;
        Real            logZ;
        Real            maxdiff;
        double          time;
        size_t          iters;
        bool            has_logZ;
        bool            has_maxdiff;
        bool            has_iters;

        TestDAI( const FactorGraph &fg, const string &_name, const PropertySet &opts ) : obj(NULL), name(_name), err(), varmargs(), marginals(), logZ(0.0), maxdiff(0.0), time(0), iters(0U), has_logZ(false), has_maxdiff(false), has_iters(false) {
            double tic = toc();
            if( name == "LDPC" ) {
                Prob zero(2,0.0);
                zero[0] = 1.0;
                for( size_t i = 0; i < fg.nrVars(); i++ )
                    varmargs.push_back( Factor(fg.var(i), zero) );
                marginals = varmargs;
                logZ = 0.0;
                maxdiff = 0.0;
                iters = 1;
                has_logZ = false;
                has_maxdiff = false;
                has_iters = false;
            } else
                obj = newInfAlg( name, fg, opts );
            time += toc() - tic;
        }

        ~TestDAI() {
            if( obj != NULL )
                delete obj;
        }

        string identify() const {
            if( obj != NULL )
                return obj->identify();
            else
                return "NULL";
        }

        void doDAI() {
            double tic = toc();
            if( obj != NULL ) {
                obj->init();
                obj->run();
                time += toc() - tic;

                try {
                    logZ = obj->logZ();
                    has_logZ = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_logZ = false;
                    else
                        throw;
                }

                try {
                    maxdiff = obj->maxDiff();
                    has_maxdiff = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_maxdiff = false;
                    else
                        throw;
                }

                try {
                    iters = obj->Iterations();
                    has_iters = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_iters = false;
                    else
                        throw;
                }

                varmargs.clear();
                for( size_t i = 0; i < obj->fg().nrVars(); i++ )
                    varmargs.push_back( obj->beliefV( i ) );

                marginals = obj->beliefs();
            };
        }

        void calcErrs( const TestDAI &x ) {
            err.clear();
            err.reserve( varmargs.size() );
            for( size_t i = 0; i < varmargs.size(); i++ )
                err.push_back( dist( varmargs[i], x.varmargs[i], Prob::DISTTV ) );
        }

        void calcErrs( const vector<Factor> &x ) {
            err.clear();
            err.reserve( varmargs.size() );
            for( size_t i = 0; i < varmargs.size(); i++ )
                err.push_back( dist( varmargs[i], x[i], Prob::DISTTV ) );
        }

        Real maxErr() {
            return( *max_element( err.begin(), err.end() ) );
        }

        Real avgErr() {
            return( accumulate( err.begin(), err.end(), 0.0 ) / err.size() );
        }
};


pair<string, PropertySet> parseMethodRaw( const string &s ) {
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


pair<string, PropertySet> parseMethod( const string &_s, const map<string,string> & aliases ) {
    // break string into method[properties]
    pair<string,PropertySet> ps = parseMethodRaw(_s);
    bool looped = false;

    // as long as 'method' is an alias, update:
    while( aliases.find(ps.first) != aliases.end() && !looped ) {
        string astr = aliases.find(ps.first)->second;
        pair<string,PropertySet> aps = parseMethodRaw(astr);
        if( aps.first == ps.first )
            looped = true;
        // override aps properties by ps properties
        aps.second.Set( ps.second );
        // replace ps by aps
        ps = aps;
        // repeat until method name == alias name ('looped'), or
        // there is no longer an alias 'method'
    }

    // check whether name is valid
    size_t n = 0;
    for( ; strlen( DAINames[n] ) != 0; n++ )
        if( ps.first == DAINames[n] )
            break;
    if( strlen( DAINames[n] ) == 0 && (ps.first != "LDPC") )
        DAI_THROWE(UNKNOWN_DAI_ALGORITHM,string("Unknown DAI algorithm \"") + ps.first + string("\" in \"") + _s + string("\""));

    return ps;
}


Real clipReal( Real x, Real minabs ) {
    if( abs(x) < minabs )
        return minabs;
    else
        return x;
}


DAI_ENUM(MarginalsOutputType,NONE,VAR,ALL);


int main( int argc, char *argv[] ) {
    string filename;
    string aliases;
    vector<string> methods;
    Real tol;
    size_t maxiter;
    size_t verbose;
    MarginalsOutputType marginals;
    bool report_iters = true;
    bool report_time = true;

    po::options_description opts_required("Required options");
    opts_required.add_options()
        ("filename", po::value< string >(&filename), "Filename of FactorGraph")
        ("methods", po::value< vector<string> >(&methods)->multitoken(), "DAI methods to test")
    ;

    po::options_description opts_optional("Allowed options");
    opts_optional.add_options()
        ("help", "produce help message")
        ("aliases", po::value< string >(&aliases), "Filename for aliases")
        ("tol", po::value< Real >(&tol), "Override tolerance")
        ("maxiter", po::value< size_t >(&maxiter), "Override maximum number of iterations")
        ("verbose", po::value< size_t >(&verbose), "Override verbosity")
        ("marginals", po::value< MarginalsOutputType >(&marginals), "Output marginals? (NONE,VAR,ALL)")
        ("report-time", po::value< bool >(&report_time), "Report calculation time")
        ("report-iters", po::value< bool >(&report_iters), "Report iterations needed")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(opts_required).add(opts_optional);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    if( vm.count("help") || !(vm.count("filename") && vm.count("methods")) ) {
        cout << "Reads factorgraph <filename.fg> and performs the approximate" << endl;
        cout << "inference algorithms <method*>, reporting calculation time, max and average" << endl;
        cout << "error and relative logZ error (comparing with the results of" << endl;
        cout << "<method0>, the base method, for which one can use JTREE_HUGIN)." << endl << endl;
        cout << opts_required << opts_optional << endl;
#ifdef DAI_DEBUG
        cout << "This is a debugging (unoptimised) build of libDAI." << endl;
#endif
        return 1;
    }

    try {
        // Read aliases
        map<string,string> Aliases;
        if( !aliases.empty() ) {
            ifstream infile;
            infile.open (aliases.c_str());
            if (infile.is_open()) {
                while( true ) {
                    string line;
                    getline(infile,line);
                    if( infile.fail() )
                        break;
                    if( (!line.empty()) && (line[0] != '#') ) {
                        string::size_type pos = line.find(':',0);
                        if( pos == string::npos )
                            DAI_THROWE(RUNTIME_ERROR,"Invalid alias");
                        else {
                            string::size_type posl = line.substr(0, pos).find_last_not_of(" \t");
                            string key = line.substr(0, posl + 1);
                            string::size_type posr = line.substr(pos + 1, line.length()).find_first_not_of(" \t");
                            string val = line.substr(pos + 1 + posr, line.length());
                            Aliases[key] = val;
                        }
                    }
                }
                infile.close();
            } else
                DAI_THROWE(RUNTIME_ERROR,"Error opening aliases file");
        }

        FactorGraph fg;
        fg.ReadFromFile( filename.c_str() );

        vector<Factor> varmargs0;
        Real logZ0 = 0.0;

        cout.setf( ios_base::scientific );
        cout.precision( 3 );

        cout << "# " << filename << endl;
        cout.width( 39 );
        cout << left << "# METHOD" << "\t";
        if( report_time )
            cout << right << "SECONDS  " << "\t";
        if( report_iters )
            cout << "ITERS" << "\t";
        cout << "MAX ERROR" << "\t";
        cout << "AVG ERROR" << "\t";
        cout << "LOGZ ERROR" << "\t";
        cout << "MAXDIFF" << "\t";
        cout << endl;

        for( size_t m = 0; m < methods.size(); m++ ) {
            pair<string, PropertySet> meth = parseMethod( methods[m], Aliases );

            if( vm.count("tol") )
                meth.second.Set("tol",tol);
            if( vm.count("maxiter") )
                meth.second.Set("maxiter",maxiter);
            if( vm.count("verbose") )
                meth.second.Set("verbose",verbose);
            TestDAI testdai(fg, meth.first, meth.second );
            testdai.doDAI();
            if( m == 0 ) {
                varmargs0 = testdai.varmargs;
                logZ0 = testdai.logZ;
            }
            testdai.calcErrs(varmargs0);

            cout.width( 39 );
            cout << left << methods[m] << "\t";
            if( report_time )
                cout << right << testdai.time << "\t";
            if( report_iters ) {
                if( testdai.has_iters ) {
                    cout << testdai.iters << "\t";
                } else {
                    cout << "N/A  \t";
                }
            }

            if( m > 0 ) {
                cout.setf( ios_base::scientific );
                cout.precision( 3 );

                Real me = clipReal( testdai.maxErr(), 1e-9 );
                cout << me << "\t";

                Real ae = clipReal( testdai.avgErr(), 1e-9 );
                cout << ae << "\t";

                if( testdai.has_logZ ) {
                    cout.setf( ios::showpos );
                    Real le = clipReal( testdai.logZ / logZ0 - 1.0, 1e-9 );
                    cout << le << "\t";
                    cout.unsetf( ios::showpos );
                } else
                    cout << "N/A       \t";

                if( testdai.has_maxdiff ) {
                    Real md = clipReal( testdai.maxdiff, 1e-9 );
                    if( isnan( me ) )
                        md = me;
                    if( isnan( ae ) )
                        md = ae;
                    if( md == INFINITY )
                        md = 1.0;
                    cout << md << "\t";
                } else
                    cout << "N/A    \t";
            }
            cout << endl;

            if( marginals == MarginalsOutputType::VAR ) {
                for( size_t i = 0; i < testdai.varmargs.size(); i++ )
                    cout << "# " << testdai.varmargs[i] << endl;
            } else if( marginals == MarginalsOutputType::ALL ) {
                for( size_t I = 0; I < testdai.marginals.size(); I++ )
                    cout << "# " << testdai.marginals[I] << endl;
            }
        }

        return 0;
    } catch( string &s ) {
        cerr << "Exception: " << s << endl;
        return 2;
    }
}
