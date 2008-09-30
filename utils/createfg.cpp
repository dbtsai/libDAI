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
#include <vector>
#include <iterator>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <dai/factorgraph.h>
#include <dai/weightedgraph.h>
#include <dai/util.h>


using namespace std;
using namespace dai;
namespace po = boost::program_options;
typedef boost::numeric::ublas::compressed_matrix<double> matrix;
typedef matrix::value_array_type::const_iterator matrix_vcit;
typedef matrix::index_array_type::const_iterator matrix_icit;


Factor BinaryFactor( const Var &n, double field ) {
    assert( n.states() == 2 );
    double buf[2];
    buf[0] = exp(-field); 
    buf[1] = exp(field);
    return Factor(n, &buf[0]);
}


Factor BinaryFactor( const Var &n1, const Var &n2, double coupling ) {
    assert( n1.states() == 2 );
    assert( n2.states() == 2 );
    assert( n1 != n2 );
    double buf[4];
    buf[0] = (buf[3] = exp(coupling));
    buf[1] = (buf[2] = exp(-coupling));
    return Factor( VarSet(n1, n2), &buf[0] );
}


Factor RandomFactor( const VarSet &ns, double beta ) {
    Factor fac( ns );
    for( size_t t = 0; t < fac.states(); t++ )
        fac[t] = exp(rnd_stdnormal() * beta);
    return fac;
}


Factor PottsFactor( const Var &n1, const Var &n2, double beta ) {
    Factor fac( VarSet( n1, n2 ), 1.0 );
    assert( n1.states() == n2.states() );
    for( size_t s = 0; s < n1.states(); s++ )
        fac[ s * (n1.states() + 1) ] = exp(beta);
    return fac;
}


void MakeHOIFG( size_t N, size_t M, size_t k, double sigma, FactorGraph &fg ) {
	vector<Var> vars;
	vector<Factor> factors;

    vars.reserve(N);
	for( size_t i = 0; i < N; i++ )
		vars.push_back(Var(i,2));

	for( size_t I = 0; I < M; I++ ) {
		VarSet vars;
		while( vars.size() < k ) {
			do {
				size_t newind = (size_t)(N * rnd_uniform());
				Var newvar = Var(newind, 2);
				if( !vars.contains( newvar ) ) {
					vars |= newvar;
					break;
				}
			} while( 1 );
		}
        factors.push_back( RandomFactor( vars, sigma ) );
	}

    fg = FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


// w should be upper triangular or lower triangular
void WTh2FG( const matrix &w, const vector<double> &th, FactorGraph &fg ) {
    vector<Var>    vars;
    vector<Factor> factors;

    size_t N = th.size();
    assert( (w.size1() == N) && (w.size2() == N) );

    vars.reserve(N);
    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i,2));

    factors.reserve( w.nnz() + N );
    // walk through the sparse array structure
    // this is similar to matlab sparse arrays
    // index2 gives the column index (similar to ir in matlab)
    // index1 gives the starting indices for each row (similar to jc in matlab)
    size_t i = 0;
    for( size_t pos = 0; pos < w.nnz(); pos++ ) {
        while( pos == w.index1_data()[i+1] )
            i++;
        size_t j = w.index2_data()[pos];
        double w_ij = w.value_data()[pos];
        factors.push_back( BinaryFactor( vars[i], vars[j], w_ij ) );
    }
    for( size_t i = 0; i < N; i++ )
        factors.push_back( BinaryFactor( vars[i], th[i] ) );

    fg = FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


void MakeFullFG( size_t N, double mean_w, double mean_th, double sigma_w, double sigma_th, FactorGraph &fg ) {
    matrix w(N,N,N*(N-1)/2);
    vector<double> th(N,0.0);

    for( size_t i = 0; i < N; i++ ) {
        for( size_t j = i+1; j < N; j++ )
            w(i,j) = rnd_stdnormal() * sigma_w + mean_w;
        th[i] = rnd_stdnormal() * sigma_th + mean_th;
    }

    WTh2FG( w, th, fg );
}


void Make3DPotts( size_t n1, size_t n2, size_t n3, size_t states, double beta, FactorGraph &fg ) {
    vector<Var> vars;
    vector<Factor> factors;
    
    for( size_t i1 = 0; i1 < n1; i1++ )
        for( size_t i2 = 0; i2 < n2; i2++ )
            for( size_t i3 = 0; i3 < n3; i3++ ) {
                vars.push_back( Var( i1*n2*n3 + i2*n3 + i3, states ) );
                if( i1 )
                    factors.push_back( PottsFactor( vars.back(), vars[ (i1-1)*n2*n3 + i2*n3 + i3 ], beta ) );
                if( i2 )
                    factors.push_back( PottsFactor( vars.back(), vars[ i1*n2*n3 + (i2-1)*n3 + i3 ], beta ) );
                if( i3 )
                    factors.push_back( PottsFactor( vars.back(), vars[ i1*n2*n3 + i2*n3 + (i3-1) ], beta ) );
            }

    fg = FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


void MakeGridFG( long periodic, size_t n, double mean_w, double mean_th, double sigma_w, double sigma_th, FactorGraph &fg ) {
    size_t N = n*n;

    matrix w(N,N,2*N);
    vector<double> th(N,0.0);

    for( size_t i = 0; i < n; i++ )
        for( size_t j = 0; j < n; j++ ) {
            if( i+1 < n || periodic )
                w(i*n+j, ((i+1)%n)*n+j) = rnd_stdnormal() * sigma_w + mean_w;
            if( j+1 < n || periodic )
                w(i*n+j, i*n+((j+1)%n)) = rnd_stdnormal() * sigma_w + mean_w;
            th[i*n+j] = rnd_stdnormal() * sigma_th + mean_th;
        }

    WTh2FG( w, th, fg );
}

            
void MakeGridNonbinaryFG( bool periodic, size_t n, size_t states, double beta, FactorGraph &fg ) {
    size_t N = n*n;

    vector<Var>    vars;
    vector<Factor> factors;

    vars.reserve(N);
    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i, states));

    factors.reserve( 2 * N );
    for( size_t i = 0; i < n; i++ ) {
        for( size_t j = 0; j < n; j++ ) {
            if( i+1 < n || periodic )
                factors.push_back( RandomFactor( VarSet( vars[i*n+j], vars[((i+1)%n)*n+j] ), beta ) );
            if( j+1 < n || periodic )
                factors.push_back( RandomFactor( VarSet( vars[i*n+j], vars[i*n+((j+1)%n)] ), beta ) );
        }
    }

    fg = FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


void MakeLoopFG( size_t N, double mean_w, double mean_th, double sigma_w, double sigma_th, FactorGraph &fg ) {
    matrix w(N,N,N);
    vector<double> th(N,0.0);

    for( size_t i = 0; i < N; i++ ) {
        w(i, (i+1)%N) = rnd_stdnormal() * sigma_w + mean_w;
        th[i] = rnd_stdnormal() * sigma_th + mean_th;
    }

    WTh2FG( w, th, fg );
}


void MakeLoopNonbinaryFG( size_t N, size_t states, double beta, FactorGraph &fg ) {
    vector<Var>    vars;
    vector<Factor> factors;

    vars.reserve(N);
    for( size_t i = 0; i < N; i++ )
        vars.push_back(Var(i, states));

    factors.reserve( N );
    for( size_t i = 0; i < N; i++ ) {
        factors.push_back( RandomFactor( VarSet( vars[i], vars[(i+1)%N] ), beta ) );
    }

    fg = FactorGraph( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
}


void MakeTreeFG( size_t N, double mean_w, double mean_th, double sigma_w, double sigma_th, FactorGraph &fg ) {
    matrix w(N,N,N-1);
    vector<double> th(N,0.0);

    for( size_t i = 0; i < N; i++ ) {
        th[i] = rnd_stdnormal() * sigma_th + mean_th;
        if( i > 0 ) {
            size_t j = rnd_int( 0, i-1 );
            w(i,j) = rnd_stdnormal() * sigma_w + mean_w;
        }
    }

    WTh2FG( w, th, fg );
}


void MakeDRegFG( size_t N, size_t d, double mean_w, double mean_th, double sigma_w, double sigma_th, FactorGraph &fg ) {
    matrix w(N,N,(d*N)/2);
    vector<double> th(N,0.0);

    UEdgeVec g = RandomDRegularGraph( N, d );
    for( size_t i = 0; i < g.size(); i++ )
        w(g[i].n1, g[i].n2) = rnd_stdnormal() * sigma_w + mean_w;
    
    for( size_t i = 0; i < N; i++ )
        th[i] = rnd_stdnormal() * sigma_th + mean_th;

    WTh2FG( w, th, fg );
}


// N = number of variables
// n = size of variable neighborhoods
// K = number of factors
// k = size of factor neighborhoods
// asserts: N * n == K * k
BipartiteGraph CreateRandomBipartiteGraph( size_t N, size_t K, size_t n, size_t k ) {
    BipartiteGraph G;

    assert( N * n == K * k );

    // build lists of degree-repeated vertex numbers
    std::vector<size_t> stubs1(N*n,0);
    for( size_t i = 0; i < N; i++ )
        for( size_t t = 0; t < n; t++ )
            stubs1[i*n + t] = i;

    // build lists of degree-repeated vertex numbers
    std::vector<size_t> stubs2(K*k,0);
    for( size_t I = 0; I < K; I++ )
        for( size_t t = 0; t < k; t++ )
            stubs2[I*k + t] = I;

    // shuffle lists
    random_shuffle( stubs1.begin(), stubs1.end() );
    random_shuffle( stubs2.begin(), stubs2.end() );

    // add edges
    vector<BipartiteGraph::Edge> edges;
    edges.reserve( N*n );
    for( size_t e = 0; e < N*n; e++ )
        edges.push_back( BipartiteGraph::Edge(stubs1[e], stubs2[e]) );

    // finish construction
    G.construct( N, K, edges.begin(), edges.end() );

    return G;
}


// Returns x**n % p, assuming p is prime
int powmod (int x, int n, int p) {
    int y = 1;
    for( int m = 0; m < n; m++ )
        y = (x * y) % p;
    return y;
}


// Returns order of x in GF(p) with p prime
size_t order (int x, int p) {
    x = x % p;
    assert( x != 0 );
    size_t n = 0;
    size_t prod = 1;
    do {
        prod = (prod * x) % p;
        n++;
    } while( prod != 1 ); 
    return n;
}


// Returns whether n is a prime number
bool isPrime (size_t n) {
    bool result = true;
    for( size_t k = 2; (k < n) && result; k++ )
        if( n % k == 0 )
            result = false;
    return result;
}


// Make a regular LDPC graph with N=6, j=2, K=4, k=3
BipartiteGraph CreateSmallLDPCGraph() {
    BipartiteGraph G;
    size_t N=4, j=3, K=4; // k=3;

    typedef BipartiteGraph::Edge Edge;
    vector<Edge> edges;
    edges.reserve( N*j );
    edges.push_back( Edge(0,0) ); edges.push_back( Edge(1,0) ); edges.push_back( Edge(2,0) );
    edges.push_back( Edge(0,1) ); edges.push_back( Edge(1,1) ); edges.push_back( Edge(3,1) );
    edges.push_back( Edge(0,2) ); edges.push_back( Edge(2,2) ); edges.push_back( Edge(3,2) );
    edges.push_back( Edge(1,3) ); edges.push_back( Edge(2,3) ); edges.push_back( Edge(3,3) );

    // finish construction
    G.construct( N, K, edges.begin(), edges.end() );

    return G;
}


// Use construction described in "A Class of Group-Structured LDPC Codes"
// by R. M. Tanner, D. Sridhara and T. Fuja
// Proceedings of ICSTA, 2001
//
// Example parameters: (p,j,k) = (31,3,5)
// j and k must be divisors of p-1
BipartiteGraph CreateGroupStructuredLDPCGraph( size_t p, size_t j, size_t k ) {
    BipartiteGraph G;

    size_t n = j;
    size_t N = p * k;
    size_t K = p * j;

    size_t a, b;
    for( a = 2; a < p; a++ )
        if( order(a,p) == k )
            break;
    assert( a != p );
    for( b = 2; b < p; b++ )
        if( order(b,p) == j )
            break;
    assert( b != p );
//    cout << "# order(a=" << a << ") = " << order(a,p) << endl;
//    cout << "# order(b=" << b << ") = " << order(b,p) << endl;

    assert( N * n == K * k );

    typedef BipartiteGraph::Edge Edge;
    vector<Edge> edges;
    edges.reserve( N * n );

    for( size_t s = 0; s < j; s++ )
        for( size_t t = 0; t < k; t++ ) {
            size_t P = (powmod(b,s,p) * powmod(a,t,p)) % p;
            for( size_t m = 0; m < p; m++ )
                edges.push_back( Edge(t*p + m, s*p + ((m + P) % p)) );
        }

    // finish construction
    G.construct( N, K, edges.begin(), edges.end() );

    return G;
}


// Make parity check table
void MakeParityCheck( double *result, size_t n, double eps ) {
    size_t N = 1 << n;
    for( size_t i = 0; i < N; i++ ) {
        size_t c = 0;
        for( size_t t = 0; t < n; t++ )
            if( i & (1 << t) )
                c ^= 1;
        if( c )
            result[i] = eps;
        else
            result[i] = 1.0 - eps;
    }
    return;
}


const char *FULL_TYPE = "full";
const char *GRID_TYPE = "grid";
const char *GRID_TORUS_TYPE = "grid_torus";
const char *DREG_TYPE = "dreg";
const char *HOI_TYPE = "hoi";
const char *LDPC_RANDOM_TYPE = "ldpc_random";
const char *LDPC_GROUP_TYPE = "ldpc_group";
const char *LDPC_SMALL_TYPE = "ldpc_small";
const char *LOOP_TYPE = "loop";
const char *TREE_TYPE = "tree";
const char *POTTS3D_TYPE = "potts3d";


int main( int argc, char *argv[] ) {
    try {
		size_t N, K, k, d, j, n1, n2, n3;
        size_t prime;
        size_t seed;
		double beta, sigma_w, sigma_th, noise, mean_w, mean_th;
        string type;
        size_t states = 2;

        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help",     "produce help message")
            ("type",     po::value<string>(&type),     "factor graph type:\n\t'full', 'grid', 'grid_torus', 'dreg', 'loop', 'tree', 'hoi', 'ldpc_random', 'ldpc_group', 'ldpc_small', 'potts3d'")
            ("seed",     po::value<size_t>(&seed),     "random number seed (tries to read from /dev/urandom if not specified)")
            ("N",        po::value<size_t>(&N),        "number of variables (not for type=='ldpc_small')")
            ("n1",       po::value<size_t>(&n1),       "width of 3D grid (only for type=='potts3d')")
            ("n2",       po::value<size_t>(&n2),       "height of 3D grid (only for type=='potts3d')")
            ("n3",       po::value<size_t>(&n3),       "length of 3D grid (only for type=='potts3d')")
            ("K",        po::value<size_t>(&K),        "number of factors\n\t(only for type=='hoi' and 'type=='ldpc_{random,group}')")
            ("k",        po::value<size_t>(&k),        "number of variables per factor\n\t(only for type=='hoi' and type=='ldpc_{random,group}')")
            ("d",        po::value<size_t>(&d),        "variable connectivity\n\t(only for type=='dreg')")
            ("j",        po::value<size_t>(&j),        "number of parity checks per bit\n\t(only for type=='ldpc_{random,group}')")
            ("prime",    po::value<size_t>(&prime),    "prime number for construction of LDPC code\n\t(only for type=='ldpc_group')")
            ("beta",     po::value<double>(&beta),     "stddev of log-factor entries\n\t(only for type=='hoi', 'potts3d', 'grid' if states>2)")
            ("mean_w",   po::value<double>(&mean_w),   "mean of pairwise interactions w_{ij}\n\t(not for type=='hoi', 'ldpc_*', 'potts3d')")
            ("mean_th",  po::value<double>(&mean_th),  "mean of singleton interactions th_i\n\t(not for type=='hoi', 'ldpc_*', 'potts3d')")
            ("sigma_w",  po::value<double>(&sigma_w),  "stddev of pairwise interactions w_{ij}\n\t(not for type=='hoi', 'ldpc_*', 'potts3d')")
            ("sigma_th", po::value<double>(&sigma_th), "stddev of singleton interactions th_i\n\t(not for type=='hoi', 'ldpc_*', 'potts3d'")
            ("noise",    po::value<double>(&noise),    "bitflip probability for binary symmetric channel (only for type=='ldpc')")
            ("states",   po::value<size_t>(&states),   "number of states of each variable (should be 2 for all but type=='grid', 'grid_torus', 'loop', 'potts3d')")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if( vm.count("help") || !vm.count("type") ) {
            if( vm.count("type") ) {
                if( type == FULL_TYPE ) {
                    cout << "Creates fully connected pairwise graphical model of <N> binary variables;" << endl;
                } else if( type == GRID_TYPE ) {
                    cout << "Creates (non-periodic) 2D Ising grid of (approx.) <N> variables (which need not be binary);" << endl;
                } else if( type == GRID_TORUS_TYPE ) {
                    cout << "Creates periodic 2D Ising grid of (approx.) <N> variables (which need not be binary);" << endl;
                } else if( type == DREG_TYPE ) {
                    cout << "Creates random d-regular graph of <N> binary variables with uniform degree <d>" << endl;
                    cout << "(where <d><N> should be even);" << endl;
                } else if( type == LOOP_TYPE ) {
                    cout << "Creates a pairwise graphical model consisting of a single loop of" << endl;
                    cout << "<N> variables (which need not be binary);" << endl;
                } else if( type == TREE_TYPE ) {
                    cout << "Creates a pairwise, connected graphical model without cycles (i.e., a tree)" << endl;
                    cout << "of <N> binary variables;" << endl;
                } else if( type == HOI_TYPE ) {
                    cout << "Creates a random factor graph of <N> binary variables and" << endl;
                    cout << "<K> factors, each factor being an interaction of <k> variables." << endl;
                    cout << "The entries of the factors are exponentials of i.i.d. Gaussian" << endl;
                    cout << "variables with mean 0 and standard deviation <beta>." << endl;
                } else if( type == LDPC_RANDOM_TYPE ) {
                    cout << "Simulates LDPC decoding problem, using a LDPC code of <N> bits and <K> parity" << endl;
                    cout << "checks, with <k> bits per check and <j> checks per bit, transmitted on a binary" << endl;
                    cout << "symmetric channel with probability <noise> of flipping a bit. The transmitted" << endl;
                    cout << "codeword has all bits set to zero. The LDPC code is randomly generated." << endl;
                } else if( type == LDPC_GROUP_TYPE ) {
                    cout << "Simulates LDPC decoding problem, using a LDPC code of <N> bits and <K> parity" << endl;
                    cout << "checks, with <k> bits per check and <j> checks per bit, transmitted on a binary" << endl;
                    cout << "symmetric channel with probability <noise> of flipping a bit. The transmitted" << endl;
                    cout << "codeword has all bits set to zero. The LDPC code is constructed (using group" << endl;
                    cout << "theory) using a parameter <prime>; <j> and <k> should both be divisors of <prime>-1." << endl;
                } else if( type == LDPC_SMALL_TYPE ) {
                    cout << "Simulates LDPC decoding problem, using a LDPC code of 4 bits and 4 parity" << endl;
                    cout << "checks, with 3 bits per check and 3 checks per bit, transmitted on a binary" << endl;
                    cout << "symmetric channel with probability <noise> of flipping a bit. The transmitted" << endl;
                    cout << "codeword has all bits set to zero. The LDPC code is fixed." << endl;
                } else if( type == POTTS3D_TYPE ) {
                    cout << "Builds 3D Potts model of size <n1>x<n2>x<n3> with nearest-neighbour Potts" << endl;
                    cout << "interactions with <states> states and inverse temperature <beta>." << endl;
                } else
                    cerr << "Unknown type (should be one of 'full', 'grid', 'grid_torus', 'dreg', 'loop', 'tree', 'hoi', 'ldpc_random', 'ldpc_group', 'ldpc_small', 'potts3d')" << endl;
                
                if( type == FULL_TYPE || type == GRID_TYPE || type == GRID_TORUS_TYPE || type == DREG_TYPE || type == LOOP_TYPE || type == TREE_TYPE ) {
                    if( type == GRID_TYPE || type == GRID_TORUS_TYPE || type == LOOP_TYPE ) {
                        cout << "if <states> > 2: factor entries are exponents of Gaussians with mean 0 and standard deviation beta; otherwise," << endl;
                    }
                    cout << "singleton interactions are Gaussian with mean <mean_th> and standard" << endl;
                    cout << "deviation <sigma_th>; pairwise interactions are Gaussian with mean" << endl;
                    cout << "<mean_w> and standard deviation <sigma_w>." << endl;
                } 
            }
            cout << endl << desc << endl;
            return 1;
        }

        if( !vm.count("states") )
            states = 2;

        if( !vm.count("seed") ) {
            ifstream infile;
            bool success;
            infile.open( "/dev/urandom" );
            success = infile.is_open();
            if( success ) {
                infile.read( (char *)&seed, sizeof(size_t) / sizeof(char) );
                success = infile.good();
                infile.close();
            }
            if( !success )
                throw "Please specify random number seed.";
        }
        rnd_seed( seed );

        FactorGraph fg;

        cout << "# Factor graph made by " << argv[0] << endl;
        cout << "# type = " << type << endl;

        if( type == FULL_TYPE ) {
            if( !vm.count("N") || !vm.count("mean_w") || !vm.count("mean_th") || !vm.count("sigma_w") || !vm.count("sigma_th") )
                throw "Please specify all required arguments";
            MakeFullFG( N, mean_w, mean_th, sigma_w, sigma_th, fg );

            cout << "# N = " << N << endl;
            cout << "# mean_w = " << mean_w << endl;
            cout << "# mean_th = " << mean_th << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
        } else if( type == GRID_TYPE || type == GRID_TORUS_TYPE ) {
            if( states > 2 ) {
                if( !vm.count("N") || !vm.count("beta") )
                    throw "Please specify all required arguments";
            } else {
                if( !vm.count("N") || !vm.count("mean_w") || !vm.count("mean_th") || !vm.count("sigma_w") || !vm.count("sigma_th") )
                    throw "Please specify all required arguments";
            }

            size_t n = (size_t)sqrt((long double)N);
            N = n * n;

            bool periodic = false;
            if( type == GRID_TYPE )
                periodic = false;
            else
                periodic = true;

            if( states > 2 )
                MakeGridNonbinaryFG( periodic, n, states, beta, fg );
            else
                MakeGridFG( periodic, n, mean_w, mean_th, sigma_w, sigma_th, fg );
                
            cout << "# n = " << n << endl;
            cout << "# N = " << N << endl;
            
            if( states > 2 )
                cout << "# beta = " << beta << endl;
            else {
                cout << "# mean_w = " << mean_w << endl;
                cout << "# mean_th = " << mean_th << endl;
                cout << "# sigma_w = " << sigma_w << endl;
                cout << "# sigma_th = " << sigma_th << endl;
            }
        } else if( type == DREG_TYPE ) {
            if( !vm.count("N") || !vm.count("mean_w") || !vm.count("mean_th") || !vm.count("sigma_w") || !vm.count("sigma_th") || !vm.count("d") )
                throw "Please specify all required arguments";

            MakeDRegFG( N, d, mean_w, mean_th, sigma_w, sigma_th, fg );

            cout << "# N = " << N << endl;
            cout << "# d = " << d << endl;
            cout << "# mean_w = " << mean_w << endl;
            cout << "# mean_th = " << mean_th << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
        } else if( type == LOOP_TYPE ) {
            if( states > 2 ) {
                if( !vm.count("N") || !vm.count("beta") )
                    throw "Please specify all required arguments";
            } else {
                if( !vm.count("N") || !vm.count("mean_w") || !vm.count("mean_th") || !vm.count("sigma_w") || !vm.count("sigma_th") )
                    throw "Please specify all required arguments";
            }
            if( states > 2 )
                MakeLoopNonbinaryFG( N, states, beta, fg );
            else
                MakeLoopFG( N, mean_w, mean_th, sigma_w, sigma_th, fg );

            cout << "# N = " << N << endl;

            if( states > 2 )
                cout << "# beta = " << beta << endl;
            else {
                cout << "# mean_w = " << mean_w << endl;
                cout << "# mean_th = " << mean_th << endl;
                cout << "# sigma_w = " << sigma_w << endl;
                cout << "# sigma_th = " << sigma_th << endl;
            }
        } else if( type == TREE_TYPE ) {
            if( !vm.count("N") || !vm.count("mean_w") || !vm.count("mean_th") || !vm.count("sigma_w") || !vm.count("sigma_th") )
                throw "Please specify all required arguments";
            MakeTreeFG( N, mean_w, mean_th, sigma_w, sigma_th, fg );

            cout << "# N = " << N << endl;
            cout << "# mean_w = " << mean_w << endl;
            cout << "# mean_th = " << mean_th << endl;
            cout << "# sigma_w = " << sigma_w << endl;
            cout << "# sigma_th = " << sigma_th << endl;
        } else if( type == HOI_TYPE ) {
            if( !vm.count("N") || !vm.count("K") || !vm.count("k") || !vm.count("beta") )
                throw "Please specify all required arguments";
            do {
                MakeHOIFG( N, K, k, beta, fg );
            } while( !fg.isConnected() );

            cout << "# N = " << N << endl;
            cout << "# K = " << K << endl;
            cout << "# k = " << k << endl;
            cout << "# beta = " << beta << endl;
        } else if( type == LDPC_RANDOM_TYPE || type == LDPC_GROUP_TYPE || type == LDPC_SMALL_TYPE ) {
            if( !vm.count("noise") )
                throw "Please specify all required arguments";

            if( type == LDPC_RANDOM_TYPE ) {
                if( !vm.count("N") || !vm.count("K") || !vm.count("j") || !vm.count("k") )
                    throw "Please specify all required arguments";

                if( N * j != K * k )
                    throw "Parameters should satisfy N * j == K * k";
            } else if( type == LDPC_GROUP_TYPE ) {
                if( !vm.count("prime") || !vm.count("j") || !vm.count("k") )
                    throw "Please specify all required arguments";

                if( !isPrime(prime) )
                    throw "Parameter <prime> should be prime";
                if( !((prime-1) % j == 0 ) )
                    throw "Parameters should satisfy (prime-1) % j == 0";
                if( !((prime-1) % k == 0 ) )
                    throw "Parameters should satisfy (prime-1) % k == 0";

                N = prime * k;
                K = prime * j;
            } else if( type == LDPC_SMALL_TYPE ) {
                N = 4;
                K = 4;
                j = 3;
                k = 3;
            }

            cout << "# N = " << N << endl;
            cout << "# K = " << K << endl;
            cout << "# j = " << j << endl;
            cout << "# k = " << k << endl;
            if( type == LDPC_GROUP_TYPE )
                cout << "# prime = " << prime << endl;
            cout << "# noise = " << noise << endl;

            // p = 31, j = 3, k = 5
            // p = 37, j = 3, k = 4
            // p = 7 , j = 2, k = 3
            // p = 29, j = 2, k = 4

            // Construct likelihood and paritycheck factors
            double likelihood[4] = {1.0 - noise, noise, noise, 1.0 - noise};
            double *paritycheck = new double[1 << k];
            MakeParityCheck(paritycheck, k, 0.0);

            // Create LDPC structure
            BipartiteGraph ldpcG;
            bool regular;
            do {
                if( type == LDPC_GROUP_TYPE )
                    ldpcG = CreateGroupStructuredLDPCGraph( prime, j, k );
                else if( type == LDPC_RANDOM_TYPE )
                    ldpcG = CreateRandomBipartiteGraph( N, K, j, k );
                else if( type == LDPC_SMALL_TYPE )
                    ldpcG = CreateSmallLDPCGraph();

                regular = true;
                for( size_t i = 0; i < N; i++ )
                    if( ldpcG.nb1(i).size() != j )
                        regular = false;
                for( size_t I = 0; I < K; I++ )
                    if( ldpcG.nb2(I).size() != k )
                        regular = false;
            } while( !regular && !ldpcG.isConnected() );

            // Convert to FactorGraph
            vector<Factor> factors;
            for( size_t I = 0; I < K; I++ ) {
                VarSet vs;
                for( size_t _i = 0; _i < k; _i++ ) {
                    size_t i = ldpcG.nb2(I)[_i];
                    vs |= Var( i, 2 );
                }
                factors.push_back( Factor( vs, paritycheck ) );
            }
            delete paritycheck;

            // Generate noise vector
            vector<char> noisebits(N,0);
            size_t bitflips = 0;
            for( size_t i = 0; i < N; i++ ) {
                if( rnd_uniform() < noise ) {
                    noisebits[i] = 1;
                    bitflips++;
                }
            }
            cout << "# bitflips = " << bitflips << endl;

            // Simulate transmission of all-zero codeword
            vector<char> input(N,0);
            vector<char> output(N,0);
            for( size_t i = 0; i < N; i++ )
                output[i] = (input[i] + noisebits[i]) & 1;

            // Add likelihoods
            for( size_t i = 0; i < N; i++ )
               factors.push_back( Factor(Var(i,2), likelihood + output[i]*2) ); 

            // Construct Factor Graph
            fg = FactorGraph( factors );
        } else if( type == POTTS3D_TYPE ) {
            if( !vm.count("n1") || !vm.count("n2") || !vm.count("n3") || !vm.count("beta") || !vm.count("states") )
                throw "Please specify all required arguments";
            Make3DPotts( n1, n2, n3, states, beta, fg );

            cout << "# N = " << n1*n2*n3 << endl;
            cout << "# n1 = " << n1 << endl;
            cout << "# n2 = " << n2 << endl;
            cout << "# n3 = " << n3 << endl;
            cout << "# beta = " << beta << endl;
            cout << "# states = " << states << endl;
        } else {
            throw "Invalid type";
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
