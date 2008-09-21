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


#include <iostream>
#include <fstream>
#include <map>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <boost/program_options.hpp>
#include <dai/util.h>
#include <dai/alldai.h>


using namespace std;
using namespace dai;
namespace po = boost::program_options;


class TestAI {
    protected:
        InfAlg          *obj;
        string          name;
        vector<double>  err;

    public:
        vector<Factor>  q;
        double          logZ;
        double          maxdiff;
        double          time;

        TestAI( const FactorGraph &fg, const string &_name, const PropertySet &opts ) : obj(NULL), name(_name), err(), q(), logZ(0.0), maxdiff(0.0), time(0) {
            double tic = toc();
            obj = newInfAlg( name, fg, opts );
            time += toc() - tic;
/*
            } else if( method.substr(0,5) == "EXACT" ) { // EXACT
                // Look if the network is small enough to do brute-force exact method
                bool toolarge = false;
                size_t total_statespace = 1;
                for( size_t i = 0; i < fg.nrVars(); i++ ) {
                    total_statespace *= fg.var(i).states();
                    if( total_statespace > (1UL << 16) )
                        toolarge = true;
                }

                if( !toolarge ) {
                    Factor piet;
                    for( size_t I = 0; I < fg.nrFactors(); I++ )
                        piet *= fg.factor( I );
                    for( size_t i = 0; i < fg.nrVars(); i++ )
                        q.push_back(piet.marginal(fg.var(i)));
                    time += toc() - tic;
                    logZ = fg.ExactlogZ();
                } else
                    throw "Network too large for EXACT method";
            }
*/
        }

        ~TestAI() { 
            if( obj != NULL )
                delete obj;
        }

        string identify() { 
            if( obj != NULL )
                return obj->identify(); 
            else
                return "NULL";
        }

        vector<Factor> allBeliefs() {
            vector<Factor> result;
            for( size_t i = 0; i < obj->nrVars(); i++ )
                result.push_back( obj->belief( obj->var(i) ) );
            return result;
        }

        void doAI() {
            double tic = toc();
//            if( name == "EXACT" ) {
//                // calculation has already been done
//            } 
            if( obj != NULL ) {
                obj->init();
                obj->run();
                time += toc() - tic;
                logZ = real(obj->logZ());
                maxdiff = obj->maxDiff();
                q = allBeliefs();
            };
        }

        void calcErrs( const TestAI &x ) {
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


pair<string, PropertySet> parseMethod( const string &_s, const map<string,string> & aliases ) {
    string s = _s;
    if( aliases.find(_s) != aliases.end() )
        s = aliases.find(_s)->second;

    pair<string, PropertySet> result;
    string & name = result.first;
    PropertySet & opts = result.second;

    string::size_type pos = s.find_first_of('[');
    name = s.substr( 0, pos );
    if( pos == string::npos )
        throw "Malformed method";
    size_t n = 0;
    for( ; n < sizeof(DAINames) / sizeof(string); n++ )
        if( name == DAINames[n] )
            break;
    if( n == sizeof(DAINames) / sizeof(string) )
        throw "Unknown inference algorithm";

    stringstream ss;
    ss << s.substr(pos,s.length());
    ss >> opts;
    
    return result;
}


double clipdouble( double x, double minabs ) {
    if( fabs(x) < minabs )
        return minabs;
    else
        return x;
}


int main( int argc, char *argv[] ) {
    try {
        string filename;
        string aliases;
        vector<string> methods;
        double tol;
        size_t maxiter;
        size_t verbose;

        po::options_description opts_required("Required options");
        opts_required.add_options()
            ("filename", po::value< string >(&filename), "Filename of FactorGraph")
            ("methods", po::value< vector<string> >(&methods)->multitoken(), "AI methods to test")
        ;

        po::options_description opts_optional("Allowed options");
        opts_optional.add_options()
            ("help", "produce help message")
            ("aliases", po::value< string >(&aliases), "Filename for aliases")
            ("tol", po::value< double >(&tol), "Override tolerance")
            ("maxiter", po::value< size_t >(&maxiter), "Override maximum number of iterations")
            ("verbose", po::value< size_t >(&verbose), "Override verbosity")
        ;

        po::options_description cmdline_options;
        cmdline_options.add(opts_required).add(opts_optional);

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
        po::notify(vm);

        if( vm.count("help") || !(vm.count("filename") && vm.count("methods")) ) {
            cout << "Reads factorgraph <filename.fg> and performs the approximate" << endl;
            cout << "inference algorithms <method*>, reporting clocks, max and average" << endl;
            cout << "error and relative logZ error (comparing with the results of" << endl;
            cout << "<method0>, the base method, for which one can use JTREE_HUGIN)." << endl << endl;
            cout << opts_required << opts_optional << endl;
            return 1;
        }

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
                            throw "Invalid alias";
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
                throw "Error opening aliases file";
        }

        FactorGraph fg;
        if( fg.ReadFromFile(filename.c_str()) ) {
            cout << "Error reading " << filename << endl;
            return 2;
        } else {
            vector<Factor> q0;
            double logZ0 = 0.0;

            cout << "# " << filename << endl;
            cout.width( 40 );
            cout << left << "# METHOD" << "  ";
            cout.width( 10 );
            cout << right << "SECONDS" << "   ";
            cout.width( 10 );
            cout << "MAX ERROR" << "  ";
            cout.width( 10 );
            cout << "AVG ERROR" << "  ";
            cout.width( 10 );
            cout << "LOGZ ERROR" << "  ";
            cout.width( 10 );
            cout << "MAXDIFF" << endl;

            for( size_t m = 0; m < methods.size(); m++ ) {
                pair<string, PropertySet> meth = parseMethod( methods[m], Aliases );

                if( vm.count("tol") )
                    meth.second.Set("tol",tol);
                if( vm.count("maxiter") )
                    meth.second.Set("maxiter",maxiter);
                if( vm.count("verbose") )
                    meth.second.Set("verbose",verbose);
                TestAI piet(fg, meth.first, meth.second );
                piet.doAI();
                if( m == 0 ) {
                    q0 = piet.q;
                    logZ0 = piet.logZ;
                }
                piet.calcErrs(q0);

                cout.width( 40 );
//                cout << left << piet.identify() << "  ";
                cout << left << methods[m] << "  ";
                cout.width( 10 );
                cout << right << piet.time << "    ";

                if( m > 0 ) {
                    cout.setf( ios_base::scientific );
                    cout.precision( 3 );
                    cout.width( 10 ); 
                    double me = clipdouble( piet.maxErr(), 1e-9 );
                    cout << me << "  ";
                    cout.width( 10 );
                    double ae = clipdouble( piet.avgErr(), 1e-9 );
                    cout << ae << "  ";
                    cout.width( 10 );
                    double le = clipdouble( piet.logZ / logZ0 - 1.0, 1e-9 );
                    cout << le << "  ";
                    cout.width( 10 );
                    double md = clipdouble( piet.maxdiff, 1e-9 );
                    if( isnan( me ) )
                        md = me;
                    if( isnan( ae ) )
                        md = ae;
                    cout << md << endl;
                } else
                    cout << endl;
            }
        }
    } catch(const char *e) {
        cerr << "Exception: " << e << endl;
        return 1;
    } catch(exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
    }

    return 0;
}
