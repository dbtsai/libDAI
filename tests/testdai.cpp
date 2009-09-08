/*  Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
    Radboud University Nijmegen, The Netherlands /
    Max Planck Institute for Biological Cybernetics, Germany

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
        vector<double>  err;

    public:
        vector<Factor>  q;
        double          logZ;
        double          maxdiff;
        double          time;
        size_t          iters;
        bool            has_logZ;
        bool            has_maxdiff;
        bool            has_iters;

        TestDAI( const FactorGraph &fg, const string &_name, const PropertySet &opts ) : obj(NULL), name(_name), err(), q(), logZ(0.0), maxdiff(0.0), time(0), iters(0U), has_logZ(false), has_maxdiff(false), has_iters(false) {
            double tic = toc();
            if( name == "LDPC" ) {
                double zero[2] = {1.0, 0.0};
                q.clear();
                for( size_t i = 0; i < fg.nrVars(); i++ )
                    q.push_back( Factor(Var(i,2), zero) );
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

        vector<Factor> allBeliefs() {
            vector<Factor> result;
            for( size_t i = 0; i < obj->fg().nrVars(); i++ )
                result.push_back( obj->belief( obj->fg().var(i) ) );
            return result;
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

                q = allBeliefs();
            };
        }

        void calcErrs( const TestDAI &x ) {
            err.clear();
            err.reserve( q.size() );
            for( size_t i = 0; i < q.size(); i++ )
                err.push_back( dist( q[i], x.q[i], Prob::DISTTV ) );
        }

        void calcErrs( const vector<Factor> &x ) {
            err.clear();
            err.reserve( q.size() );
            for( size_t i = 0; i < q.size(); i++ )
                err.push_back( dist( q[i], x[i], Prob::DISTTV ) );
        }

        double maxErr() {
            return( *max_element( err.begin(), err.end() ) );
        }

        double avgErr() {
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


double clipdouble( double x, double minabs ) {
    if( fabs(x) < minabs )
        return minabs;
    else
        return x;
}


int main( int argc, char *argv[] ) {
    string filename;
    string aliases;
    vector<string> methods;
    double tol;
    size_t maxiter;
    size_t verbose;
    bool marginals = false;
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
        ("tol", po::value< double >(&tol), "Override tolerance")
        ("maxiter", po::value< size_t >(&maxiter), "Override maximum number of iterations")
        ("verbose", po::value< size_t >(&verbose), "Override verbosity")
        ("marginals", po::value< bool >(&marginals), "Output single node marginals?")
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

        vector<Factor> q0;
        double logZ0 = 0.0;

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
            TestDAI piet(fg, meth.first, meth.second );
            piet.doDAI();
            if( m == 0 ) {
                q0 = piet.q;
                logZ0 = piet.logZ;
            }
            piet.calcErrs(q0);

            cout.width( 39 );
            cout << left << methods[m] << "\t";
            if( report_time )
                cout << right << piet.time << "\t";
            if( report_iters ) {
                if( piet.has_iters ) {
                    cout << piet.iters << "\t";
                } else {
                    cout << "N/A  \t";
                }
            }

            if( m > 0 ) {
                cout.setf( ios_base::scientific );
                cout.precision( 3 );

                double me = clipdouble( piet.maxErr(), 1e-9 );
                cout << me << "\t";

                double ae = clipdouble( piet.avgErr(), 1e-9 );
                cout << ae << "\t";

                if( piet.has_logZ ) {
                    cout.setf( ios::showpos );
                    double le = clipdouble( piet.logZ / logZ0 - 1.0, 1e-9 );
                    cout << le << "\t";
                    cout.unsetf( ios::showpos );
                } else
                    cout << "N/A       \t";

                if( piet.has_maxdiff ) {
                    double md = clipdouble( piet.maxdiff, 1e-9 );
                    if( isnan( me ) )
                        md = me;
                    if( isnan( ae ) )
                        md = ae;
                    cout << md << "\t";
                } else
                    cout << "N/A    \t";
            }
            cout << endl;

            if( marginals ) {
                for( size_t i = 0; i < piet.q.size(); i++ )
                    cout << "# " << piet.q[i] << endl;
            }
        }

        return 0;
    } catch( string &s ) {
        cerr << "Exception: " << s << endl;
        return 2;
    }
}
