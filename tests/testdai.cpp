/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
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


/// Wrapper class for DAI approximate inference algorithms
class TestDAI {
    protected:
        /// Stores a pointer to an InfAlg object, managed by this class
        InfAlg          *obj;
        /// Stores the name of the InfAlg algorithm
        string          name;
        /// Stores the total variation distances of the variable marginals
        vector<Real>    err;

    public:
        /// Stores the variable marginals
        vector<Factor>  varMarginals;
        /// Stores all marginals
        vector<Factor>  allMarginals;
        /// Stores the logarithm of the partition sum
        Real            logZ;
        /// Stores the maximum difference in the last iteration
        Real            maxdiff;
        /// Stores the computation time (in seconds)
        double          time;
        /// Stores the number of iterations needed
        size_t          iters;
        /// Does the InfAlg support logZ()?
        bool            has_logZ;
        /// Does the InfAlg support maxDiff()?
        bool            has_maxdiff;
        /// Does the InfAlg support Iterations()?
        bool            has_iters;

        /// Construct from factor graph \a fg, name \a _name, and set of properties \a opts
        TestDAI( const FactorGraph &fg, const string &_name, const PropertySet &opts ) : obj(NULL), name(_name), err(), varMarginals(), allMarginals(), logZ(0.0), maxdiff(0.0), time(0), iters(0U), has_logZ(false), has_maxdiff(false), has_iters(false) {
            double tic = toc();

            if( name == "LDPC" ) {
                // special case: simulating a Low Density Parity Check code
                Real zero[2] = {1.0, 0.0};
                for( size_t i = 0; i < fg.nrVars(); i++ )
                    varMarginals.push_back( Factor(fg.var(i), zero) );
                allMarginals = varMarginals;
                logZ = 0.0;
                maxdiff = 0.0;
                iters = 1;
                has_logZ = false;
                has_maxdiff = false;
                has_iters = false;
            } else
                // create a corresponding InfAlg object
                obj = newInfAlg( name, fg, opts );

            // Add the time needed to create the object
            time += toc() - tic;
        }

        /// Destructor
        ~TestDAI() {
            if( obj != NULL )
                delete obj;
        }

        /// Identify
        string identify() const {
            if( obj != NULL )
                return obj->identify();
            else
                return "NULL";
        }

        /// Run the algorithm and store its results
        void doDAI() {
            double tic = toc();
            if( obj != NULL ) {
                // Initialize
                obj->init();
                // Run
                obj->run();
                // Record the time
                time += toc() - tic;

                // Store logarithm of the partition sum (if supported)
                try {
                    logZ = obj->logZ();
                    has_logZ = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_logZ = false;
                    else
                        throw;
                }

                // Store maximum difference encountered in last iteration (if supported)
                try {
                    maxdiff = obj->maxDiff();
                    has_maxdiff = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_maxdiff = false;
                    else
                        throw;
                }

                // Store number of iterations needed (if supported)
                try {
                    iters = obj->Iterations();
                    has_iters = true;
                } catch( Exception &e ) {
                    if( e.code() == Exception::NOT_IMPLEMENTED )
                        has_iters = false;
                    else
                        throw;
                }

                // Store variable marginals
                varMarginals.clear();
                for( size_t i = 0; i < obj->fg().nrVars(); i++ )
                    varMarginals.push_back( obj->beliefV( i ) );

                // Store all marginals calculated by the method
                allMarginals = obj->beliefs();
            };
        }

        /// Calculate total variation distance of variable marginals with respect to those in \a x
        void calcErrs( const TestDAI &x ) {
            err.clear();
            err.reserve( varMarginals.size() );
            for( size_t i = 0; i < varMarginals.size(); i++ )
                err.push_back( dist( varMarginals[i], x.varMarginals[i], Prob::DISTTV ) );
        }

        /// Calculate total variation distance of variable marginals with respect to those in \a x
        void calcErrs( const vector<Factor> &x ) {
            err.clear();
            err.reserve( varMarginals.size() );
            for( size_t i = 0; i < varMarginals.size(); i++ )
                err.push_back( dist( varMarginals[i], x[i], Prob::DISTTV ) );
        }

        /// Return maximum error
        Real maxErr() {
            return( *max_element( err.begin(), err.end() ) );
        }

        /// Return average error
        Real avgErr() {
            return( accumulate( err.begin(), err.end(), 0.0 ) / err.size() );
        }
};


/// Clips a real number: if the absolute value of \a x is less than \a minabs, return \a minabs, else return \a x
Real clipReal( Real x, Real minabs ) {
    if( abs(x) < minabs )
        return minabs;
    else
        return x;
}


/// Whether to output no marginals, only variable marginals, or all calculated marginals
DAI_ENUM(MarginalsOutputType,NONE,VAR,ALL);


/// Main function
int main( int argc, char *argv[] ) {
    // Variables to store command line options
    // Filename of factor graph
    string filename;
    // Filename for aliases
    string aliases;
    // Approximate Inference methods to use
    vector<string> methods;
    // Which marginals to output
    MarginalsOutputType marginals;
    // Output number of iterations?
    bool report_iters = true;
    // Output calculation time?
    bool report_time = true;

    // Define required command line options
    po::options_description opts_required("Required options");
    opts_required.add_options()
        ("filename", po::value< string >(&filename), "Filename of factor graph")
        ("methods", po::value< vector<string> >(&methods)->multitoken(), "DAI methods to perform")
    ;

    // Define allowed command line options
    po::options_description opts_optional("Allowed options");
    opts_optional.add_options()
        ("help", "Produce help message")
        ("aliases", po::value< string >(&aliases), "Filename for aliases")
        ("marginals", po::value< MarginalsOutputType >(&marginals), "Output marginals? (NONE/VAR/ALL, default=NONE)")
        ("report-time", po::value< bool >(&report_time), "Output calculation time (default==1)?")
        ("report-iters", po::value< bool >(&report_iters), "Output iterations needed (default==1)?")
    ;

    // Define all command line options
    po::options_description cmdline_options;
    cmdline_options.add(opts_required).add(opts_optional);

    // Parse command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    // Display help message if necessary
    if( vm.count("help") || !(vm.count("filename") && vm.count("methods")) ) {
        cout << "This program is part of libDAI - http://www.libdai.org/" << endl << endl;
        cout << "Usage: ./testdai --filename <filename.fg> --methods <method1> [<method2> <method3> ...]" << endl << endl;
        cout << "Reads factor graph <filename.fg> and performs the approximate inference algorithms" << endl;
        cout << "<method*>, reporting for each method:" << endl;
        cout << "  o the calculation time needed, in seconds (if report-time == 1);" << endl;
        cout << "  o the number of iterations needed (if report-iters == 1);" << endl;
        cout << "  o the maximum (over all variables) total variation error in the variable marginals;" << endl;
        cout << "  o the average (over all variables) total variation error in the variable marginals;" << endl;
        cout << "  o the error (difference) of the logarithm of the partition sums;" << endl << endl;
        cout << "All errors are calculated by comparing the results of the current method with" << endl; 
        cout << "the results of the first method (the base method). If marginals==VAR, additional" << endl;
        cout << "output consists of the variable marginals, and if marginals==ALL, all marginals" << endl;
        cout << "calculated by the method are reported." << endl << endl;
        cout << "<method*> should be a list of one or more methods, seperated by spaces, in the format:" << endl << endl;
        cout << "    name[key1=val1,key2=val2,key3=val3,...,keyn=valn]" << endl << endl;
        cout << "where name should be the name of an algorithm in libDAI (or an alias, if an alias" << endl;
        cout << "filename is provided), followed by a list of properties (surrounded by rectangular" << endl;
        cout << "brackets), where each property consists of a key=value pair and the properties are" << endl;
        cout << "seperated by commas. If an alias file is specified, alias substitution is performed." << endl;
        cout << "This is done by looking up the name in the alias file and substituting the alias" << endl;
        cout << "by its corresponding method as defined in the alias file. Properties are parsed from" << endl;
        cout << "left to right, so if a property occurs repeatedly, the right-most value is used." << endl << endl;
        cout << opts_required << opts_optional << endl;
#ifdef DAI_DEBUG
        cout << "Note: this is a debugging build of libDAI." << endl << endl;
#endif
        cout << "Example:  ./testdai --filename testfast.fg --aliases aliases.conf --methods JTREE_HUGIN BP_SEQFIX BP_PARALL[maxiter=5]" << endl;
        return 1;
    }

    try {
        // Read aliases
        map<string,string> Aliases;
        if( !aliases.empty() )
            Aliases = readAliasesFile( aliases );

        // Read factor graph
        FactorGraph fg;
        fg.ReadFromFile( filename.c_str() );

        // Declare variables used for storing variable marginals and log partition sum of base method
        vector<Factor> varMarginals0;
        Real logZ0 = 0.0;

        // Output header
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

        // For each method...
        for( size_t m = 0; m < methods.size(); m++ ) {
            // Parse method
            pair<string, PropertySet> meth = parseNameProperties( methods[m], Aliases );

            // Check whether name is valid
            size_t n = 0;
            for( ; strlen( DAINames[n] ) != 0; n++ )
                if( meth.first == DAINames[n] )
                    break;
            if( strlen( DAINames[n] ) == 0 )
                DAI_THROWE(UNKNOWN_DAI_ALGORITHM,string("Unknown DAI algorithm \"") + meth.first + string("\" in \"") + methods[m] + string("\""));

            // Construct object for running the method
            TestDAI testdai(fg, meth.first, meth.second );

            // Run the method
            testdai.doDAI();

            // For the base method, store its variable marginals and logarithm of the partition sum
            if( m == 0 ) {
                varMarginals0 = testdai.varMarginals;
                logZ0 = testdai.logZ;
            }

            // Calculate errors relative to base method
            testdai.calcErrs( varMarginals0 );

            // Output method name
            cout.width( 39 );
            cout << left << methods[m] << "\t";
            // Output calculation time, if requested
            if( report_time )
                cout << right << testdai.time << "\t";
            // Output number of iterations, if requested
            if( report_iters ) {
                if( testdai.has_iters ) {
                    cout << testdai.iters << "\t";
                } else {
                    cout << "N/A  \t";
                }
            }

            // If this is not the base method
            if( m > 0 ) {
                cout.setf( ios_base::scientific );
                cout.precision( 3 );

                // Output maximum error in variable marginals
                Real me = clipReal( testdai.maxErr(), 1e-9 );
                cout << me << "\t";

                // Output average error in variable marginals
                Real ae = clipReal( testdai.avgErr(), 1e-9 );
                cout << ae << "\t";

                // Output error in log partition sum
                if( testdai.has_logZ ) {
                    cout.setf( ios::showpos );
                    Real le = clipReal( testdai.logZ - logZ0, 1e-9 );
                    cout << le << "\t";
                    cout.unsetf( ios::showpos );
                } else
                    cout << "N/A       \t";

                // Output maximum difference in last iteration
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

            // Output marginals, if requested
            if( marginals == MarginalsOutputType::VAR ) {
                for( size_t i = 0; i < testdai.varMarginals.size(); i++ )
                    cout << "# " << testdai.varMarginals[i] << endl;
            } else if( marginals == MarginalsOutputType::ALL ) {
                for( size_t I = 0; I < testdai.allMarginals.size(); I++ )
                    cout << "# " << testdai.allMarginals[I] << endl;
            }
        }

        return 0;
    } catch( string &s ) {
        // Abort with error message
        cerr << "Exception: " << s << endl;
        return 2;
    }
}
