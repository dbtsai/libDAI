/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2008-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 */


/** \file
 *  \brief Contains additional doxygen documentation
 *
 *  \todo Improve documentation
 *
 *  \todo Merge COPYING into doxygen documentation
 *  \todo Merge README into doxygen documentation
 *  \todo Document examples, tests and utils
 *
 *  \todo Adapt (part of the) guidelines in http://www.boost.org/development/requirements.html#Design_and_Programming
 *
 *  \todo Use "gcc -MM" to generate dependencies for targets: http://make.paulandlesley.org/autodep.html
 *  \todo Investigate whether switching to cmake as cross-platform build system would be a good idea.
 *
 *  \todo Replace VarSets by SmallSet<size_t> where appropriate, in order to minimize the use of FactorGraph::findVar().
 *
 *  \idea Disentangle structures. In particular, ensure that graphical properties are not
 *  entangled with probabilistic properties. For example, a FactorGraph contains several
 *  components:
 *  - a BipartiteGraph
 *  - an array of variable labels
 *  - an array of variable state space sizes
 *  - an array of pointers to factor value vectors
 *  In this way, each factor could be implemented differently, e.g., we could have
 *  some sparse factors, some noisy-OR factors, some dense factors, some arbitrary
 *  precision factors, etc.
 *
 *  \idea Use Boost::uBLAS framework to deal with matrices, especially, with 2D sparse matrices.
 *  See http://www.boost.org/libs/numeric/ublas/doc/matrix_sparse.htm
 *  I read somewhere that boost::uBLAS concentrates more on correct implementation than on performance.
 *
 *  \idea Introduce naming scheme:
 *  - all Vars should be named v_..., e.g. v_i instead of i
 *  - all VarSets should be named vs_..., e.g. v_i instead of i
 *  - all Factors should be named f_..., e.g. f_I instead of I
 *  - all indices should be named _..., e.g. _k instead of k
 *  - all iterators should be named i_, e.g. i_i is an iterator to i
 *  - all const_iterators should be named ci_, e.g. ci_i is an iterator to i
 **/


/** \page discussion Discussion of possible improvements
 *  \section discuss_extendedgraphs Extended factorgraphs/regiongraphs
 *
 *  A FactorGraph and a RegionGraph are often equipped with
 *  additional properties for nodes and edges. The code to initialize those
 *  is often quite similar. Maybe one could abstract this, e.g.:
 *  \code
 *  template <typename Node1Properties, typename Node2Properties, typename EdgeProperties>
 *  class ExtFactorGraph : public FactorGraph {
 *      public:
 *          std::vector<Node1Properties>              node1Props;
 *          std::vector<Node2Properties>              node2Props;
 *          std::vector<std::vector<EdgeProperties> > edgeProps;
 *         // ...
 *  }
 *  \endcode
 *
 *  Advantages:
 *  - Less code duplication.
 *  - Easier maintainability.
 *  - Easier to write new inference algorithms.
 *
 *  Disadvantages:
 *  - Cachability may be worse.
 *  - A problem is the case where there are no properties for either type of nodes or for edges.
 *    Maybe this can be solved using specializations, or using variadac template arguments?
 *    Another possible solution would be to define a "class Empty {}", and add some code
 *    that checks for the typeid, comparing it with Empty, and doing something special in that case
 *    (e.g., not allocating memory).
 *  - The main disadvantage of this approach seems to be that it leads to even more entanglement.
 *    Therefore this is probably a bad idea.
 *
 *  \section discuss_templates Polymorphism by template parameterization
 *  Instead of polymorphism by inheritance, use polymorphism by template parameterization.
 *  For example, the real reason for introducing the complicated inheritance scheme of dai::InfAlg
 *  was for functions like dai::calcMarginal. Instead, one could use a template function:
 *  \code
 *  template<typename InfAlg>
 *  Factor calcMarginal( const InfAlg &obj, const VarSet &ns, bool reInit );
 *  \endcode
 *  This would assume that the type InfAlg supports certain methods. Ideally, one would use
 *  concepts to define different classes of inference algorithms with different capabilities,
 *  for example the ability to calculate logZ, the ability to calculate marginals, the ability to
 *  calculate bounds, the ability to calculate MAP states, etc. Then, one would use traits
 *  classes in order to be able to query the capabilities of the model. For example, one would be
 *  able to query whether the inference algorithm supports calculation of logZ.  Unfortunately,
 *  this is compile-time polymorphism, whereas tests/testdai needs runtime polymorphism.
 *  Therefore this is probably a bad idea.
 */


/** \mainpage libDAI reference manual
 *  \author Joris Mooij
 *  \version git HEAD
 *  \date October 10, 2008
 *
 *  \section about About libDAI
 *  libDAI is a free/open source C++ library (licensed under GPLv2+) that provides
 *  implementations of various (approximate) inference methods for discrete
 *  graphical models. libDAI supports arbitrary factor graphs with discrete
 *  variables; this includes discrete Markov Random Fields and Bayesian
 *  Networks.
 *
 *  The library is targeted at researchers. To be able to use the library, a
 *  good understanding of graphical models is needed.
 *
 *  \section limitations Limitations
 *  libDAI is not intended to be a complete package for approximate inference.
 *  Instead, it should be considered as an "inference engine", providing
 *  various inference methods. In particular, it contains no GUI, currently
 *  only supports its own file format for input and output (although support
 *  for standard file formats may be added later), and provides very limited
 *  visualization functionalities. The only learning method supported currently
 *  is EM for learning factor parameters.
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
 *    [\ref MoK07, \ref MoR05];
 *  - Gibbs sampler;
 *  - Conditioned BP [\ref EaG09];
 *
 *  In addition, libDAI supports parameter learning of conditional probability
 *  tables by Expectation Maximization.
 *
 *  \section language Why C++?
 *  Because libDAI is implemented in C++, it is very fast compared with
 *  implementations in MatLab (a factor 1000 faster is not uncommon).
 *  libDAI does provide a MatLab interface for easy integration with MatLab.
 */


/** \example example.cpp
 */


/** \page quickstart Quick start
 *  An example program illustrating basic usage of libDAI is given in examples/example.cpp.
 */


/** \page bibliography Bibliography
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
 *
 *  \anchor EaG09 \ref EaG09
 *  F. Eaton and Z. Ghahramani (2009):
 *  "Choosing a Variable to Clamp",
 *  <em>Proceedings of the Twelfth International Conference on Artificial Intelligence and Statistics (AISTATS 2009)</em> 5:145-152
 *  http://jmlr.csail.mit.edu/proceedings/papers/v5/eaton09a/eaton09a.pdf
 *
 *  \anchor StW99 \ref StW99
 *  A. Steger and N. C. Wormald (1999):
 *  "Generating Random Regular Graphs Quickly",
 *  <em>Combinatorics, Probability and Computing</em> Vol 8, Issue 4, pp. 377-396
 *  http://www.math.uwaterloo.ca/~nwormald/papers/randgen.pdf
 *
 *  \anchor EMK06 \ref EMK06
 *  G. Elidan and I. McGraw and D. Koller (2006):
 *  "Residual Belief Propagation: Informed Scheduling for Asynchronous Message Passing",
 *  <em>Proceedings of the 22nd Annual Conference on Uncertainty in Artificial Intelligence (UAI-06)</em>
 *  http://uai.sis.pitt.edu/papers/06/UAI2006_0091.pdf
 */


/** \page fileformat libDAI factorgraph file format
 *
 *  This page describes the .fg fileformat used in libDAI to store factor graphs.
 *  Markov Random Fields are special cases of factor graphs, as are Bayesian
 *  networks. A factor graph can be specified as follows: for each factor, one has
 *  to specify which variables occur in the factor, what their respective
 *  cardinalities (i.e., number of possible values) are, and a table listing all
 *  the values of that factor for all possible configurations of these variables.
 *  A .fg file is not much more than that. It starts with a line containing the
 *  number of factors in that graph, followed by an empty line. Then all factors
 *  are specified, one block for each factor, where the blocks are seperated by
 *  empty lines. Each variable occurring in the factor graph has a unique
 *  identifier, its index (which should be a nonnegative integer). Comment lines
 *  start with #.
 *
 *  Each block starts with a line containing the number of variables in that
 *  factor. The second line contains the indices of these variables, seperated by
 *  spaces (indices are nonnegative integers and to avoid confusion, it is
 *  suggested to start counting at 0). The third line contains the number of
 *  possible values of each of these variables, also seperated by spaces. Note that
 *  there is some redundancy here, since if a variable appears in more than one
 *  factor, the cardinality of that variable appears several times in the .fg file.
 *  The fourth line contains the number of nonzero entries in the factor table.
 *  The rest of the lines contain these nonzero entries; each entry consists of a
 *  table index, followed by white-space, followed by the value corresponding to
 *  that table index. The most difficult part is getting the indexing right. The
 *  convention that is used is that the left-most variables cycle through their
 *  values the fastest (similar to MATLAB indexing of multidimensional arrays). An
 *  example block describing one factor is:
 *
 *  3\n
 *  4 8 7\n
 *  3 2 2\n
 *  11\n
 *  0 0.1\n
 *  1 3.5\n
 *  2 2.8\n
 *  3 6.3\n
 *  4 8.4\n
 *  6 7.4\n
 *  7 2.4\n
 *  8 8.9\n
 *  9 1.3\n
 *  10 1.6\n
 *  12 6.4\n
 *  11 2.6\n
 *
 *  which corresponds to the following factor:
 *
 *  \f[
 *  \begin{array}{ccc|c}
 *  x_4 & x_8 & x_7 & \mbox{value}\\
 *  \hline
 *   0 & 0 & 0  &  0.1\\
 *   1 & 0 & 0  &  3.5\\
 *   2 & 0 & 0  &  2.8\\
 *   0 & 1 & 0  &  6.3\\
 *   1 & 1 & 0  &  8.4\\
 *   2 & 1 & 0  &  0.0\\
 *   0 & 0 & 1  &  7.4\\
 *   1 & 0 & 1  &  2.4\\
 *   2 & 0 & 1  &  8.9\\
 *   0 & 1 & 1  &  1.3\\
 *   1 & 1 & 1  &  1.6\\
 *   2 & 1 & 1  &  2.6
 *  \end{array}
 *  \f]
 *
 *  Note that the value of x_4 changes fastest, followed by that of x_8, and x_7
 *  varies the slowest, corresponding to the second line of the block ("4 8 7").
 *  Further, x_4 can take on three values, and x_8 and x_7 each have two possible
 *  values, as described in the third line of the block ("3 2 2"). The table
 *  contains 11 non-zero entries (all except for the fifth entry). Note that the
 *  eleventh and twelveth entries are interchanged.
 *
 *  A final note: the internal representation in libDAI of the factor above is
 *  different, because the variables are ordered according to their indices
 *  (i.e., the ordering would be x_4 x_7 x_8) and the values of the table are
 *  stored accordingly, with the variable having the smallest index changing
 *  fastest:
 *
 *  \f[
 *  \begin{array}{ccc|c}
 *  x_4 & x_7 & x_8 & \mbox{value}\\
 *  \hline
 *   0 & 0 & 0  &  0.1\\
 *   1 & 0 & 0  &  3.5\\
 *   2 & 0 & 0  &  2.8\\
 *   0 & 1 & 0  &  7.4\\
 *   1 & 1 & 0  &  2.4\\
 *   2 & 1 & 0  &  8.9\\
 *   0 & 0 & 1  &  6.3\\
 *   1 & 0 & 1  &  8.4\\
 *   2 & 0 & 1  &  0.0\\
 *   0 & 1 & 1  &  1.3\\
 *   1 & 1 & 1  &  1.6\\
 *   2 & 1 & 1  &  2.6
 *  \end{array}
 *  \f]
 */


 /** \page license License
  *  \verbinclude COPYING
  */
