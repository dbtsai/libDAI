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
#include <dai/alldai.h>


using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if ( argc != 2 ) {
        cout << "Usage: " << argv[0] << " <filename.fg>" << endl << endl;
        cout << "Reads factor graph <filename.fg> and runs" << endl;
        cout << "Belief Propagation and JunctionTree on it." << endl << endl;
        return 1;
    } else {
        FactorGraph fg;
        fg.ReadFromFile(argv[1]);

        size_t  maxiter = 10000;
        double  tol = 1e-9;
        size_t  verb = 1;

        PropertySet opts;
        opts.Set("maxiter",maxiter);
        opts.Set("tol",tol);
        opts.Set("verbose",verb);

        JTree jt( fg, opts("updates",string("HUGIN")) );
        jt.init();
        jt.run();

        BP bp(fg, opts("updates",string("SEQFIX"))("logdomain",false));
        bp.init();
        bp.run();

        cout << "Exact single node marginals:" << endl;
        for( size_t i = 0; i < fg.nrVars(); i++ )
            cout << jt.belief(fg.var(i)) << endl;

        cout << "Approximate (loopy belief propagation) single node marginals:" << endl;
        for( size_t i = 0; i < fg.nrVars(); i++ )
            cout << bp.belief(fg.var(i)) << endl;

        cout << "Exact factor marginals:" << endl;
        for( size_t I = 0; I < fg.nrFactors(); I++ )
            cout << jt.belief(fg.factor(I).vars()) << endl;

        cout << "Approximate (loopy belief propagation) factor marginals:" << endl;
        for( size_t I = 0; I < fg.nrFactors(); I++ )
            cout << bp.belief(fg.factor(I).vars()) << "=" << bp.beliefF(I) << endl;

        cout << "Exact log partition sum: " << jt.logZ() << endl;
        cout << "Approximate (loopy belief propagation) log partition sum: " << bp.logZ() << endl;
    }

    return 0;
}
