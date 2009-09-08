/*  Copyright (C) 2009  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
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
#include <dai/bbp.h>


using namespace dai;
using namespace std;


int main( int argc, char *argv[] ) {
    if ( argc != 2 ) {
        cout << "Usage: " << argv[0] << " <filename.fg>" << endl << endl;
        cout << "Reads factor graph <filename.fg> and verifies" << endl;
        cout << "whether BBP works correctly on it." << endl << endl;
        return 1;
    } else {
        // Read FactorGraph from the file specified by the first command line argument
        FactorGraph fg;
        fg.ReadFromFile(argv[1]);

        // Set some constants
        size_t verbose = 0;
        double tol = 1.0e-9;
        size_t maxiter = 10000;
        double damping = 0.0;
        BBP::Properties::UpdateType updates = BBP::Properties::UpdateType::PAR;

        // Store the constants in a PropertySet object
        PropertySet opts;
        opts.Set("verbose",verbose);  // Verbosity (amount of output generated)
        opts.Set("tol",tol);          // Tolerance for convergence
        opts.Set("maxiter",maxiter);  // Maximum number of iterations
        opts.Set("damping",damping);  // Amount of damping applied

        // Construct a BP (belief propagation) object from the FactorGraph fg
        BP bp(fg, opts("updates",string("SEQFIX"))("logdomain",false));
        bp.recordSentMessages = true;
        bp.init();
        bp.run();

        vector<size_t> state( fg.nrVars(), 0 );

        for( size_t t = 0; t < 45; t++ ) {
            BBP::Properties::UpdateType updates;
            switch( t % 5 ) {
                case BBP::Properties::UpdateType::SEQ_FIX:
                    updates = BBP::Properties::UpdateType::SEQ_FIX;
                    break;
                case BBP::Properties::UpdateType::SEQ_MAX:
                    updates = BBP::Properties::UpdateType::SEQ_MAX;
                    break;
                case BBP::Properties::UpdateType::SEQ_BP_REV:
                    updates = BBP::Properties::UpdateType::SEQ_BP_REV;
                    break;
                case BBP::Properties::UpdateType::SEQ_BP_FWD:
                    updates = BBP::Properties::UpdateType::SEQ_BP_FWD;
                    break;
                case BBP::Properties::UpdateType::PAR:
                    updates = BBP::Properties::UpdateType::PAR;
                    break;
            }
            bbp_cfn_t cfn;
            switch( (t / 5) % 9 ) {
                case 0:
                    cfn = bbp_cfn_t::CFN_GIBBS_B;
                    break;
                case 1:
                    cfn = bbp_cfn_t::CFN_GIBBS_B2;
                    break;
                case 2:
                    cfn = bbp_cfn_t::CFN_GIBBS_EXP;
                    break;
                case 3:
                    cfn = bbp_cfn_t::CFN_GIBBS_B_FACTOR;
                    break;
                case 4:
                    cfn = bbp_cfn_t::CFN_GIBBS_B2_FACTOR;
                    break;
                case 5:
                    cfn = bbp_cfn_t::CFN_GIBBS_EXP_FACTOR;
                    break;
                case 6:
                    cfn = bbp_cfn_t::CFN_VAR_ENT;
                    break;
                case 7:
                    cfn = bbp_cfn_t::CFN_FACTOR_ENT;
                    break;
                case 8:
                    cfn = bbp_cfn_t::CFN_BETHE_ENT;
                    break;
            }

            double h = 1e-4;
            double result = numericBBPTest( bp, &state, opts("updates",updates), cfn, h );
            cout << "result: " << result << ",\tupdates=" << updates << ", cfn=" << cfn << endl;
        }
    }

    return 0;
}
