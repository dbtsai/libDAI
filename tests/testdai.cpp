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

        string identify() { 
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
                    has_logZ = false;
                }
                try {
                    maxdiff = obj->maxDiff();
                    has_maxdiff = true;
                } catch( Exception &e ) {
                    has_maxdiff = false;
                }
                try {
                    iters = obj->Iterations();
                    has_iters = true;
                } catch( Exception &e ) {
                    has_iters = false;
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


pair<string, PropertySet> parseMethod( const string &_s, const map<string,string> & aliases ) {
    // s = first part of _s, until '['
    string::size_type pos = _s.find_first_of('[');
    string s;
    if( pos == string::npos )
        s = _s;
    else
        s = _s.substr(0,pos);

    // if the first part is an alias, substitute
    if( aliases.find(s) != aliases.end() )
        s = aliases.find(s)->second;

    // attach second part, merging properties if necessary
    if( pos != string::npos ) {
        if( s.at(s.length()-1) == ']' ) {
            s = s.erase(s.length()-1,1) + ',' + _s.substr(pos+1);
        } else
            s = s + _s.substr(pos);
    }

    pair<string, PropertySet> result;
    string & name = result.first;
    PropertySet & opts = result.second;

    pos = s.find_first_of('[');
    if( pos == string::npos )
        throw "Malformed method";
    name = s.substr( 0, pos );
    size_t n = 0;
    for( ; strlen( DAINames[n] ) != 0; n++ )
        if( name == DAINames[n] )
            break;
    if( strlen( DAINames[n] ) == 0 && (name != "LDPC") )
        DAI_THROW(UNKNOWN_DAI_ALGORITHM);

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
        fg.ReadFromFile( filename.c_str() );

        vector<Factor> q0;
        double logZ0 = 0.0;

        cout.setf( ios_base::scientific );
        cout.precision( 3 );

        cout << "# " << filename << endl;
        cout.width( 40 );
        cout << left << "# METHOD" << "  ";
        if( report_time ) {
            cout.width( 10 );
            cout << right << "SECONDS" << "   ";
        }
        if( report_iters ) {
            cout.width( 10 );
            cout << "ITERS" << "  ";
        }
        cout.width( 10 );
        cout << "MAX ERROR" << "  ";
        cout.width( 10 );
        cout << "AVG ERROR" << "  ";
        cout.width( 10 );
        cout << "LOGZ ERROR" << "  ";
        cout.width( 10 );
        cout << "MAXDIFF" << "  ";
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

            cout.width( 40 );
            cout << left << methods[m] << "  ";
            if( report_time ) {
                cout.width( 10 );
                cout << right << piet.time << "    ";
            }
            if( report_iters ) {
                cout.width( 10 );
                if( piet.has_iters ) {
                    cout << piet.iters << "  ";
                } else {
                    cout << "N/A         ";
                }
            }

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
                if( piet.has_logZ ) {
                    double le = clipdouble( piet.logZ / logZ0 - 1.0, 1e-9 );
                    cout << le << "  ";
                } else
                    cout << "N/A         ";

                cout.width( 10 );
                if( piet.has_maxdiff ) {
                    double md = clipdouble( piet.maxdiff, 1e-9 );
                    if( isnan( me ) )
                        md = me;
                    if( isnan( ae ) )
                        md = ae;
                    cout << md << "  ";
                } else
                    cout << "N/A         ";
            }
            cout << endl;

            if( marginals ) {
                for( size_t i = 0; i < piet.q.size(); i++ )
                    cout << "# " << piet.q[i] << endl;
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
