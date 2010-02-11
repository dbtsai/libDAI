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
#include <cstdlib>
#include <dai/factorgraph.h>
#include <dai/jtree.h>


using namespace std;
using namespace dai;


void findLoopClusters( const FactorGraph & fg, std::set<VarSet> &allcl, VarSet newcl, const Var & root, size_t length, VarSet vars ) {
    for( VarSet::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        VarSet ind = fg.delta( *in );
        if( (newcl.size()) >= 2 && ind.contains( root ) ) {
            allcl.insert( newcl | *in );
        }
        else if( length > 1 )
            findLoopClusters( fg, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


size_t countLoops( const FactorGraph & fg, size_t loopdepth ) {
    set<VarSet> loops;
    for( vector<Var>::const_iterator i0 = fg.vars().begin(); i0 != fg.vars().end(); i0++ ) {
        VarSet i0d = fg.delta(*i0);
        if( loopdepth > 1 )
            findLoopClusters( fg, loops, *i0, *i0, loopdepth - 1, i0d );
    }
    return loops.size();
}


bool hasShortLoops( const std::vector<Factor> &P ) {
    bool found = false;
    vector<Factor>::const_iterator I, J;
    for( I = P.begin(); I != P.end(); I++ ) {
        J = I;
        J++;
        for( ; J != P.end(); J++ )
            if( (I->vars() & J->vars()).size() >= 2 ) {
                found = true;
                break;
            }
        if( found )
            break;
    }
    return found;
}


bool hasNegatives( const std::vector<Factor> &P ) {
    bool found = false;
    for( size_t I = 0; I < P.size(); I++ )
        if( P[I].hasNegatives() ) {
            found = true;
            break;
        }
    return found;
}


int main( int argc, char *argv[] ) {
    if( argc != 3 ) {
        // Display help message if number of command line arguments is incorrect
        cout << "This program is part of libDAI - http://www.libdai.org/" << endl << endl;
        cout << "Usage: ./fginfo <in.fg> <tw>" << endl << endl;
        cout << "Reports some detailed information about the factor graph <in.fg>." << endl;
        cout << "Also calculates treewidth (which may take some time) unless <tw> == 0." << endl;
        return 1;
    } else {
        // Read factorgraph
        FactorGraph fg;
        char *infile = argv[1];
        int calc_tw = atoi(argv[2]);
        fg.ReadFromFile( infile );

        // Output various statistics
        cout << "Number of variables:   " << fg.nrVars() << endl;
        cout << "Number of factors:     " << fg.nrFactors() << endl;
        cout << "Connected:             " << fg.isConnected() << endl;
        cout << "Tree:                  " << fg.isTree() << endl;
        cout << "Has short loops:       " << hasShortLoops(fg.factors()) << endl;
        cout << "Has negatives:         " << hasNegatives(fg.factors()) << endl;
        cout << "Binary variables?      " << fg.isBinary() << endl;
        cout << "Pairwise interactions? " << fg.isPairwise() << endl;
        // Calculate treewidth using various heuristics, if requested
        if( calc_tw ) {
            std::pair<size_t,size_t> tw;
            tw = boundTreewidth(fg, &eliminationCost_MinNeighbors);
            cout << "Treewidth (MinNeighbors):     " << tw.first << " (" << tw.second << " states)" << endl;
            tw = boundTreewidth(fg, &eliminationCost_MinWeight);
            cout << "Treewidth (MinWeight):        " << tw.first << " (" << tw.second << " states)" << endl;
            tw = boundTreewidth(fg, &eliminationCost_MinFill);
            cout << "Treewidth (MinFill):          " << tw.first << " (" << tw.second << " states)" << endl;
            tw = boundTreewidth(fg, &eliminationCost_WeightedMinFill);
            cout << "Treewidth (WeightedMinFill):  " << tw.first << " (" << tw.second << " states)" << endl;
        }
        // Calculate total state space
        long double stsp = 1.0;
        for( size_t i = 0; i < fg.nrVars(); i++ )
            stsp *= fg.var(i).states();
        cout << "Total state space:   " << stsp << endl;
        // Output type of factor graph
        cout << "Type: " << (fg.isPairwise() ? "pairwise" : "higher order") << " interactions, " << (fg.isBinary() ? "binary" : "nonbinary") << " variables" << endl;

        // Calculate complexity for LCBP
        long double cavsum_lcbp = 0.0;
        long double cavsum_lcbp2 = 0.0;
        size_t max_Delta_size = 0;
        map<size_t,size_t> cavsizes;
        for( size_t i = 0; i < fg.nrVars(); i++ ) {
            VarSet di = fg.delta(i);
            if( cavsizes.count(di.size()) )
                cavsizes[di.size()]++;
            else
                cavsizes[di.size()] = 1;
            size_t Ds = fg.Delta(i).nrStates();
            if( Ds > max_Delta_size )
                max_Delta_size = Ds;
            cavsum_lcbp += di.nrStates();
            for( VarSet::const_iterator j = di.begin(); j != di.end(); j++ )
                cavsum_lcbp2 += j->states();
        }
        cout << "Maximum pancake has " << max_Delta_size << " states" << endl;
        cout << "LCBP with full cavities needs " << cavsum_lcbp << " BP runs" << endl;
        cout << "LCBP with only pairinteractions needs " << cavsum_lcbp2 << " BP runs" << endl;
        cout << "Cavity sizes: ";
        for( map<size_t,size_t>::const_iterator it = cavsizes.begin(); it != cavsizes.end(); it++ )
            cout << it->first << "(" << it->second << ") ";
        cout << endl;

        // Calculate girth and length of loops
        if( fg.isPairwise() ) {
            bool girth_reached = false;
            size_t loopdepth;
            for( loopdepth = 2; loopdepth <= fg.nrVars() && !girth_reached; loopdepth++ ) {
                size_t nr_loops = countLoops( fg, loopdepth );
                cout << "Loops up to " << loopdepth << " variables: " << nr_loops << endl;
                if( nr_loops > 0 )
                    girth_reached = true;
            }
            if( girth_reached )
                cout << "Girth: " << loopdepth-1 << endl;
            else
                cout << "Girth: infinity" << endl;
        }

        // Output factor state spaces
        map<size_t,size_t> facsizes;
        for( size_t I = 0; I < fg.nrFactors(); I++ ) {
            size_t Isize = fg.factor(I).vars().size();
            if( facsizes.count( Isize ) )
                facsizes[Isize]++;
            else
                facsizes[Isize] = 1;
        }
        cout << "Factor sizes: ";
        for( map<size_t,size_t>::const_iterator it = facsizes.begin(); it != facsizes.end(); it++ )
            cout << it->first << "(" << it->second << ") ";
        cout << endl;

        return 0;
    }
}
