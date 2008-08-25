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
#include <iterator>
#include "factorgraph.h"
#include "weightedgraph.h"
#include "util.h"
#include <boost/program_options.hpp>


using namespace std;
using namespace dai;
namespace po = boost::program_options;


void MakeHOIFG( size_t N, size_t M, size_t k, double sigma, FactorGraph &fg ) {
	vector<Var> vars;
	vector<Factor> factors;

	for( size_t i = 0; i < N; i++ )
		vars.push_back(Var(i,2));

	for( size_t I = 0; I < M; I++ ) {
		VarSet vars;
		while( vars.size() < k ) {
			do {
				size_t newind = (size_t)(N * rnd_uniform());
				Var newvar = Var(newind, 2);
				if( !(vars && newvar) ) {
					vars |= newvar;
					break;
				}
			} while( 1 );
		}
		Factor newfac(vars);
		for( size_t t = 0; t < newfac.stateSpace(); t++ )
			newfac[t] = exp(rnd_stdnormal() * sigma);
		factors.push_back(newfac);
	}

    fg = FactorGraph(factors);
};


void MakeFullFG( size_t N, double sigma_w, double sigma_th, string type, FactorGraph &fg ) {
    vector<Var> vars;
    vector<Factor> factors;

    double w[N][N];
    double th[N];
    double buf[4];

    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i,2));
    
    for( size_t i = 0; i < N; i++ )
        for( size_t j = 0; j < N; j++ )
            w[i][j] = 0.0;

    for( size_t i = 0; i < N; i++ )
        for( size_t j = i+1; j < N; j++ ) {
            w[i][j] = rnd_stdnormal() * sigma_w;
            if( type == "fe" )
                w[i][j] = fabs(w[i][j]);
            else if( type == "af" )
                w[i][j] = -fabs(w[i][j]);
            w[j][i] = w[i][j];
            buf[0] = (buf[3] = exp(w[i][j]));
            buf[1] = (buf[2] = exp(-w[i][j]));
            factors.push_back(Factor(VarSet(vars[i],vars[j]),buf));
        }

    for( size_t i = 0; i < N; i++ ) {
        th[i] = rnd_stdnormal() * sigma_th;
        buf[0] = exp(th[i]);
        buf[1] = exp(-th[i]);
        factors.push_back(Factor(vars[i],buf));
    }

    fg = FactorGraph(factors);
};


void MakeGridFG( long periodic, long n, double sigma_w, double sigma_th, string type, FactorGraph &fg ) {
    vector<Var> vars;
    vector<Factor> factors;

    long N = n*n;

    double w[N][N];
    double th[N];
    double buf[4];

    for( long i = 0; i < N; i++ )
        vars.push_back(Var(i,2));
    
    for( long i = 0; i < N; i++ )
        for( long j = 0; j < N; j++ )
            w[i][j] = 0.0;

    for( long i = 0; i < n; i++ )
        for( long j = 0; j < n; j++ ) {
            if( i+1 < n || periodic )
                w[i*n+j][((i+1)%n)*n+j] = 1.0;
            if( i > 0 || periodic )
                w[i*n+j][((i+n-1)%n)*n+j] = 1.0;
            if( j+1 < n || periodic )
                w[i*n+j][i*n+((j+1)%n)] = 1.0;
            if( j > 0 || periodic )
                w[i*n+j][i*n+((j+n-1)%n)] = 1.0;
        }

    for( long i = 0; i < N; i++ )
        for( long j = i+1; j < N; j++ ) 
            if( w[i][j] ) {
                w[i][j] = rnd_stdnormal() * sigma_w;
                if( type == "fe" )
                    w[i][j] = fabs(w[i][j]);
                else if( type == "af" )
                    w[i][j] = -fabs(w[i][j]);
                w[j][i] = w[i][j];
                buf[0] = (buf[3] = exp(w[i][j]));
                buf[1] = (buf[2] = exp(-w[i][j]));
                factors.push_back(Factor(VarSet(vars[i],vars[j]),buf));
            }

    for( long i = 0; i < N; i++ ) {
        th[i] = rnd_stdnormal() * sigma_th;
        buf[0] = exp(th[i]);
        buf[1] = exp(-th[i]);
        factors.push_back(Factor(vars[i],buf));
    }

    fg = FactorGraph(factors);
};


void MakeDRegFG( size_t N, size_t d, double sigma_w, double sigma_th, string type, FactorGraph &fg ) {
    vector<Var> vars;
    vector<Factor> factors;

    double w[N][N];
    double th[N];
    double buf[4];

    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i,2));
    
    for( size_t i = 0; i < N; i++ )
        for( size_t j = 0; j < N; j++ )
            w[i][j] = 0.0;

    UEdgeVec g = RandomDRegularGraph( N, d );
    for( size_t i = 0; i < g.size(); i++ ) {
        w[g[i].n1][g[i].n2] = 1.0;
        w[g[i].n2][g[i].n1] = 1.0;
    }
    
    for( size_t i = 0; i < N; i++ )
        for( size_t j = i+1; j < N; j++ ) 
            if( w[i][j] ) {
                w[i][j] = rnd_stdnormal() * sigma_w;
                if( type == "fe" )
                    w[i][j] = fabs(w[i][j]);
                else if( type == "af" )
                    w[i][j] = -fabs(w[i][j]);
                w[j][i] = w[i][j];
                buf[0] = (buf[3] = exp(w[i][j]));
                buf[1] = (buf[2] = exp(-w[i][j]));
                factors.push_back(Factor(VarSet(vars[i],vars[j]),buf));
            }

    for( size_t i = 0; i < N; i++ ) {
        th[i] = rnd_stdnormal() * sigma_th;
        buf[0] = exp(th[i]);
        buf[1] = exp(-th[i]);
        factors.push_back(Factor(vars[i],buf));
    }

    fg = FactorGraph(factors);
};


const char *HOITYPE = "hoi";
const char *FULLTYPE = "full";
const char *GRIDTYPE = "grid";
const char *DREGTYPE = "dreg";


// Old usages:
// create_full_fg <N> <sigma_w> <sigma_th> <subtype>
// create_grid_fg <periodic> <n> <sigma_w> <sigma_th> <subtype>
// create_dreg_fg <d> <N> <sigma_w> <sigma_th> <subtype>


int main( int argc, char *argv[] ) {
    try {
		size_t N, M, k, d;
        size_t periodic;
        size_t seed;
		double beta, sigma_w, sigma_th;
        string type, subtype;

        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help",     "produce help message")
            ("type",     po::value<string>(&type),     "factor graph type:\n\t'full', 'grid', 'dreg' or 'hoi'")
            ("seed",     po::value<size_t>(&seed),     "random number seed")
            ("subtype",  po::value<string>(&subtype),  "interactions type:\n\t'sg', 'fe' or 'af'\n\t(ignored for type=='hoi')")
            ("N",        po::value<size_t>(&N),        "number of (binary) variables")
            ("M",        po::value<size_t>(&M),        "number of factors\n\t(only for type=='hoi')")
            ("k",        po::value<size_t>(&k),        "connectivity of the factors\n\t(only for type=='hoi')")
            ("d",        po::value<size_t>(&d),        "variable connectivity\n\t(only for type=='dreg')")
            ("beta",     po::value<double>(&beta),     "stddev of log-factor entries\n\t(only for type=='hoi')")
            ("sigma_w",  po::value<double>(&sigma_w),  "stddev of pairwise interactions w_{ij}\n\t(ignored for type=='hoi')")
            ("sigma_th", po::value<double>(&sigma_th), "stddev of singleton interactions th_i\n\t(ignored for type=='hoi')")
            ("periodic", po::value<size_t>(&periodic), "0/1 corresponding to nonperiodic/periodic grid\n\t(only for type=='grid')")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if( vm.count("help") || !vm.count("type") ) {
            if( vm.count("type") ) {
                if( type == HOITYPE ) {
                    cout << "Creates a random factor graph of <N> binary variables and" << endl;
                    cout << "<M> factors, each factor being an interaction of <k> variables." << endl;
                    cout << "The entries of the factors are exponentials of i.i.d. Gaussian" << endl;
                    cout << "variables with mean 0 and standard deviation <beta>." << endl;
                } else if( type == FULLTYPE ) {
                    cout << "Creates fully connected pairwise graphical model of <N> variables;" << endl;
                } else if( type == GRIDTYPE ) {
                    cout << "Creates 2D Ising grid (periodic if <periodic>!=0) of (approx.) <N> variables;" << endl;
                } else if( type == DREGTYPE ) {
                    cout << "Creates random d-regular graph of <N> nodes with uniform degree <d>" << endl;
                    cout << "(where <d><N> should be even);" << endl;
                } else
                    cerr << "Unknown type (should be one of 'full', 'grid', 'dreg' or 'hoi')" << endl;
                
                if( type == FULLTYPE || type == GRIDTYPE || type == DREGTYPE ) {
                    cout << "singleton interactions are Gaussian with mean 0 and standard" << endl;
                    cout << "deviation <sigma_th>; pairwise interactions are Gaussian with mean 0" << endl;
                    cout << "and standard deviation <sigma_w> if <subtype>=='sg', absolute value" << endl;
                    cout << "is taken if <subtype>=='fe' and a minus sign is added if <subtype>=='af'." << endl;
                }
            }
            cout << endl << desc << endl;
            return 1;
        }

        if( !vm.count("seed") )
            throw "Please specify random number seed.";
        rnd_seed( seed );
//      srand( gsl_rng_default_seed );

        FactorGraph fg;

        cout << "# Factor graph made by " << argv[0] << endl;
        cout << "# type = " << type << endl;

        if( type == HOITYPE ) {
            if( !vm.count("N") || !vm.count("M") || !vm.count("k") || !vm.count("beta") )
                throw "Please specify all required arguments";
            do {
                MakeHOIFG( N, M, k, beta, fg );
            } while( !fg.isConnected() );

            cout << "# N = " << N << endl;
            cout << "# M = " << M << endl;
            cout << "# k = " << k << endl;
            cout << "# beta = " << beta << endl;
        } else if( type == FULLTYPE ) {
            if( !vm.count("N") || !vm.count("sigma_w") || !vm.count("sigma_th") || !vm.count("subtype") )
                throw "Please specify all required arguments";
            MakeFullFG( N, sigma_w, sigma_th, subtype, fg );

            cout << "# N = " << N << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
            cout << "# subtype = " << subtype << endl;
        } else if( type == GRIDTYPE ) {
            if( !vm.count("N") || !vm.count("sigma_w") || !vm.count("sigma_th") || !vm.count("subtype") || !vm.count("periodic") )
                throw "Please specify all required arguments";

            size_t n = (size_t)sqrt((long double)N);
            N = n * n;

            MakeGridFG( periodic, n, sigma_w, sigma_th, subtype, fg );
            
            cout << "# periodic = " << periodic << endl;
            cout << "# n = " << n << endl;
            cout << "# N = " << N << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
            cout << "# subtype = " << subtype << endl;
        } else if( type == DREGTYPE ) {
            if( !vm.count("N") || !vm.count("sigma_w") || !vm.count("sigma_th") || !vm.count("subtype") || !vm.count("d") )
                throw "Please specify all required arguments";

            MakeDRegFG( N, d, sigma_w, sigma_th, subtype, fg );

            cout << "# N = " << N << endl;
            cout << "# d = " << d << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
            cout << "# subtype = " << subtype << endl;
        }

        cout << "# seed = " << seed << endl;
        cout << fg;
    }
    catch(exception& e) {
        cerr << "Error: " << e.what() << endl; 
        return 1;
    } 
    catch(const char * e) {
        cerr << "Error: " << e << endl;
        return 1;
    } 
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
    }

    return 0;
}
