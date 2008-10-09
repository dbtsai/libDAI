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


/// \file
/// \brief Main libDAI header file
/// \todo Improve documentation


#ifndef __defined_libdai_alldai_h
#define __defined_libdai_alldai_h


#include <string>
#include <dai/daialg.h>
#include <dai/properties.h>
#include <dai/exactinf.h>
#ifdef DAI_WITH_BP
    #include <dai/bp.h>
#endif
#ifdef DAI_WITH_MF
    #include <dai/mf.h>
#endif
#ifdef DAI_WITH_HAK
    #include <dai/hak.h>
#endif
#ifdef DAI_WITH_LC
    #include <dai/lc.h>
#endif
#ifdef DAI_WITH_TREEEP
    #include <dai/treeep.h>
#endif
#ifdef DAI_WITH_JTREE
    #include <dai/jtree.h>
#endif
#ifdef DAI_WITH_MR
    #include <dai/mr.h>
#endif


/// Namespace for libDAI
namespace dai {


/// Constructs a new approximate inference algorithm.
/** \param name The name of the approximate inference algorithm (should be one of the names in DAINames).
 *  \param fg The FactorGraph that the algorithm should be applied to.
 *  \param opts A PropertySet specifying the options for the algorithm.
 *  \return Returns a pointer to the new InfAlg object; it is the responsibility of the caller to delete it later.
 */
InfAlg *newInfAlg( const std::string &name, const FactorGraph &fg, const PropertySet &opts );


/// Contains the names of all approximate inference algorithms compiled into libDAI.
static const char* DAINames[] = {
    ExactInf::Name,
#ifdef DAI_WITH_BP
    BP::Name, 
#endif
#ifdef DAI_WITH_MF
    MF::Name,
#endif
#ifdef DAI_WITH_HAK
    HAK::Name,
#endif
#ifdef DAI_WITH_LC
    LC::Name,
#endif
#ifdef DAI_WITH_TREEEP
    TreeEP::Name,
#endif
#ifdef DAI_WITH_JTREE
    JTree::Name,
#endif
#ifdef DAI_WITH_MR
    MR::Name,
#endif
    ""
};


} // end of namespace dai


/** \mainpage libDAI reference manual
 *  \author Joris Mooij
 *  \version git HEAD
 *  \date October 8, 2008
 *
 *  \section about About libDAI
 *  libDAI is a free/open source C++ library (licensed under GPL) that provides
 *  implementations of various (approximate) inference methods for discrete 
 *  graphical models. libDAI supports arbitrary factor graphs with discrete 
 *  variables; this includes discrete Markov Random Fields and Bayesian 
 *  Networks.
 *
 *  The library is targeted at researchers; to be able to use the library, a 
 *  good understanding of graphical models is needed. 
 *
 *  \section limitations Limitations
 *  libDAI is not intended to be a complete package for approximate inference. 
 *  Instead, it should be considered as an "inference engine", providing 
 *  various inference methods. In particular, it contains no GUI, currently 
 *  only supports its own file format for input and output (although support 
 *  for standard file formats may be added later), and provides very limited
 *  visualization functionalities.
 *
 *  \section features Features
 *  Currently, libDAI supports the following (approximate) inference methods:
 *  - Exact inference by brute force enumeration;
 *  - Exact inference by junction-tree methods;
 *  - Mean Field;
 *  - Loopy Belief Propagation [\ref KFL01];
 *  - Tree Expectation Propagation [\ref MiQ04];
 *  - Generalized Belief Propagation [\ref YFW05];
 *  - Double-loop GBP [\ref HAK03];
 *  - Various variants of Loop Corrected Belief Propagation
 *    [\ref MoK07, \ref MoR05].
 *
 *  \section language Why C++?
 *  Because libDAI is implemented in C++, it is very fast compared with
 *  implementations in MatLab (a factor 1000 faster is not uncommon).
 *  libDAI does provide a MatLab interface for easy integration with MatLab. 
 *
 *  \section quickstart Quick start
 *  An example program illustrating basic usage of libDAI is given in examples/example.cpp.
 */

/// \example example.cpp

/** \page Bibliography
 *  \section Bibliograpy
 *  \anchor KFL01 \ref KFL01
 *  F. R. Kschischang and B. J. Frey and H.-A. Loeliger (2001):
 *  "Factor Graphs and the Sum-Product Algorithm",
 *  <em>IEEE Transactions on Information Theory</em> 47(2):498-519.
 *  http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=910572
 *
 *  \anchor MiQ04 \ref MiQ04
 *  T. Minka and Y. Qi (2004):
 *  "Tree-structured Approximations by Expectation Propagation",
 *  <em>Advances in Neural Information Processing Systems</em> (NIPS) 16.
 *  http://books.nips.cc/papers/files/nips16/NIPS2003_AA25.pdf
 *
 *  \anchor MoR05 \ref MoR05
 *  A. Montanari and T. Rizzo (2005):
 *  "How to Compute Loop Corrections to the Bethe Approximation",
 *  <em>Journal of Statistical Mechanics: Theory and Experiment</em>
 *  2005(10)-P10011.
 *  http://stacks.iop.org/1742-5468/2005/P10011
 *
 *  \anchor YFW05 \ref YFW05
 *  J. S. Yedidia and W. T. Freeman and Y. Weiss (2005):
 *  "Constructing Free-Energy Approximations and Generalized Belief Propagation Algorithms",
 *  <em>IEEE Transactions on Information Theory</em>
 *  51(7):2282-2312.
 *  http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=1459044
 *
 *  \anchor HAK03 \ref HAK03
 *  T. Heskes and C. A. Albers and H. J. Kappen (2003):
 *  "Approximate Inference and Constrained Optimization",
 *  <em>Proceedings of the 19th Annual Conference on Uncertainty in Artificial Intelligence (UAI-03)</em> pp. 313-320.
 *  http://www.snn.ru.nl/reports/Heskes.uai2003.ps.gz
 *
 *  \anchor MoK07 \ref MoK07
 *  J. M. Mooij and H. J. Kappen (2007):
 *  "Loop Corrections for Approximate Inference on Factor Graphs",
 *  <em>Journal of Machine Learning Research</em> 8:1113-1143.
 *  http://www.jmlr.org/papers/volume8/mooij07a/mooij07a.pdf
 *
 *  \anchor MoK07b \ref MoK07b
 *  J. M. Mooij and H. J. Kappen (2007):
 *  "Sufficient Conditions for Convergence of the Sum-Product Algorithm",
 *  <em>IEEE Transactions on Information Theory</em> 53(12):4422-4437.
 *  http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=4385778
 */


#endif
