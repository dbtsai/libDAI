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


/** \page fileformats libDAI file formats
 *
 *  \section fileformats-factorgraph Factor graph (.fg) file format
 *
 *  This section describes the .fg file format used in libDAI to store factor graphs.
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
 *  <pre>
 *  3
 *  4 8 7
 *  3 2 2
 *  11
 *  0 0.1
 *  1 3.5
 *  2 2.8
 *  3 6.3
 *  4 8.4
 *  6 7.4
 *  7 2.4
 *  8 8.9
 *  9 1.3
 *  10 1.6
 *  12 6.4
 *  11 2.6
 *  </pre>
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
 *  Note that the value of \f$x_4\f$ changes fastest, followed by that of \f$x_8\f$, and \f$x_7\f$
 *  varies the slowest, corresponding to the second line of the block ("4 8 7").
 *  Further, \f$x_4\f$ can take on three values, and \f$x_8\f$ and \f$x_7\f$ each have two possible
 *  values, as described in the third line of the block ("3 2 2"). The table
 *  contains 11 non-zero entries (all except for the fifth entry). Note that the
 *  eleventh and twelveth entries are interchanged.
 *
 *  A final note: the internal representation in libDAI of the factor above is
 *  different, because the variables are ordered according to their indices
 *  (i.e., the ordering would be \f$x_4 x_7 x_8\f$) and the values of the table are
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
 *
 *
 *  \section fileformats-evidence Evidence (.tab) file format
 *
 *  This page describes the .tab fileformat used in libDAI to store "evidence",
 *  i.e., a data set consisting of multiple samples, where each sample is the 
 *  observed joint state of some variables.
 *
 *  A .tab file is a tabular data file, consisting of a header line followed by
 *  one line for each data sample. Each line should have the same number of columns,
 *  where columns are separated by one tab character. Each column corresponds to 
 *  a variable. The header line consists of the variable labels (corresponding to 
 *  Var::label()). The other lines are observed joint states of the variables, i.e.,
 *  each line corresponds to a joint observation of the variables, and each column
 *  of a line contains the state of the variable associated with that column.
 *  Missing data is handled simply by having two consecutive tab characters, 
 *  without any characters in between.
 *
 *  \par  Example:
 *
 *  <pre>
 *  1       3       2
 *  0       0       1
 *  1       0       1
 *  1               1
 *  </pre>
 *
 *  This would correspond to a data set consisting of three observations concerning
 *  the variables with labels 1, 3 and 2; the first observation being
 *  \f$x_1 = 0, x_3 = 0, x_2 = 1\f$, the second observation being
 *  \f$x_1 = 1, x_3 = 0, x_2 = 1\f$, and the third observation being
 *  \f$x_1 = 1, x_2 = 1\f$ (where the state of \f$x_3\f$ is missing).
 */

/** \page license License

<b>libDAI is licensed under the GNU General Public License version 2, or
(at your option) any later version. The complete license text is
included below.</b>

\htmlonly
<pre>
		    GNU GENERAL PUBLIC LICENSE
		       Version 2, June 1991

 Copyright (C) 1989, 1991 Free Software Foundation, Inc.
                       51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

			    Preamble

  The licenses for most software are designed to take away your
freedom to share and change it.  By contrast, the GNU General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  This
General Public License applies to most of the Free Software
Foundation's software and to any other program whose authors commit to
using it.  (Some other Free Software Foundation software is covered by
the GNU Library General Public License instead.)  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
this service if you wish), that you receive source code or can get it
if you want it, that you can change the software or use pieces of it
in new free programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must show them these terms so they know their
rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  Finally, any free program is threatened constantly by software
patents.  We wish to avoid the danger that redistributors of a free
program will individually obtain patent licenses, in effect making the
program proprietary.  To prevent this, we have made it clear that any
patent must be licensed for everyone's free use or not licensed at all.

  The precise terms and conditions for copying, distribution and
modification follow.

		    GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License applies to any program or other work which contains
a notice placed by the copyright holder saying it may be distributed
under the terms of this General Public License.  The "Program", below,
refers to any such program or work, and a "work based on the Program"
means either the Program or any derivative work under copyright law:
that is to say, a work containing the Program or a portion of it,
either verbatim or with modifications and/or translated into another
language.  (Hereinafter, translation is included without limitation in
the term "modification".)  Each licensee is addressed as "you".

Activities other than copying, distribution and modification are not
covered by this License; they are outside its scope.  The act of
running the Program is not restricted, and the output from the Program
is covered only if its contents constitute a work based on the
Program (independent of having been made by running the Program).
Whether that is true depends on what the Program does.

  1. You may copy and distribute verbatim copies of the Program's
source code as you receive it, in any medium, provided that you
conspicuously and appropriately publish on each copy an appropriate
copyright notice and disclaimer of warranty; keep intact all the
notices that refer to this License and to the absence of any warranty;
and give any other recipients of the Program a copy of this License
along with the Program.

You may charge a fee for the physical act of transferring a copy, and
you may at your option offer warranty protection in exchange for a fee.

  2. You may modify your copy or copies of the Program or any portion
of it, thus forming a work based on the Program, and copy and
distribute such modifications or work under the terms of Section 1
above, provided that you also meet all of these conditions:

    a) You must cause the modified files to carry prominent notices
    stating that you changed the files and the date of any change.

    b) You must cause any work that you distribute or publish, that in
    whole or in part contains or is derived from the Program or any
    part thereof, to be licensed as a whole at no charge to all third
    parties under the terms of this License.

    c) If the modified program normally reads commands interactively
    when run, you must cause it, when started running for such
    interactive use in the most ordinary way, to print or display an
    announcement including an appropriate copyright notice and a
    notice that there is no warranty (or else, saying that you provide
    a warranty) and that users may redistribute the program under
    these conditions, and telling the user how to view a copy of this
    License.  (Exception: if the Program itself is interactive but
    does not normally print such an announcement, your work based on
    the Program is not required to print an announcement.)

These requirements apply to the modified work as a whole.  If
identifiable sections of that work are not derived from the Program,
and can be reasonably considered independent and separate works in
themselves, then this License, and its terms, do not apply to those
sections when you distribute them as separate works.  But when you
distribute the same sections as part of a whole which is a work based
on the Program, the distribution of the whole must be on the terms of
this License, whose permissions for other licensees extend to the
entire whole, and thus to each and every part regardless of who wrote it.

Thus, it is not the intent of this section to claim rights or contest
your rights to work written entirely by you; rather, the intent is to
exercise the right to control the distribution of derivative or
collective works based on the Program.

In addition, mere aggregation of another work not based on the Program
with the Program (or with a work based on the Program) on a volume of
a storage or distribution medium does not bring the other work under
the scope of this License.

  3. You may copy and distribute the Program (or a work based on it,
under Section 2) in object code or executable form under the terms of
Sections 1 and 2 above provided that you also do one of the following:

    a) Accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of Sections
    1 and 2 above on a medium customarily used for software interchange; or,

    b) Accompany it with a written offer, valid for at least three
    years, to give any third party, for a charge no more than your
    cost of physically performing source distribution, a complete
    machine-readable copy of the corresponding source code, to be
    distributed under the terms of Sections 1 and 2 above on a medium
    customarily used for software interchange; or,

    c) Accompany it with the information you received as to the offer
    to distribute corresponding source code.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form with such
    an offer, in accord with Subsection b above.)

The source code for a work means the preferred form of the work for
making modifications to it.  For an executable work, complete source
code means all the source code for all modules it contains, plus any
associated interface definition files, plus the scripts used to
control compilation and installation of the executable.  However, as a
special exception, the source code distributed need not include
anything that is normally distributed (in either source or binary
form) with the major components (compiler, kernel, and so on) of the
operating system on which the executable runs, unless that component
itself accompanies the executable.

If distribution of executable or object code is made by offering
access to copy from a designated place, then offering equivalent
access to copy the source code from the same place counts as
distribution of the source code, even though third parties are not
compelled to copy the source along with the object code.

  4. You may not copy, modify, sublicense, or distribute the Program
except as expressly provided under this License.  Any attempt
otherwise to copy, modify, sublicense or distribute the Program is
void, and will automatically terminate your rights under this License.
However, parties who have received copies, or rights, from you under
this License will not have their licenses terminated so long as such
parties remain in full compliance.

  5. You are not required to accept this License, since you have not
signed it.  However, nothing else grants you permission to modify or
distribute the Program or its derivative works.  These actions are
prohibited by law if you do not accept this License.  Therefore, by
modifying or distributing the Program (or any work based on the
Program), you indicate your acceptance of this License to do so, and
all its terms and conditions for copying, distributing or modifying
the Program or works based on it.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the
original licensor to copy, distribute or modify the Program subject to
these terms and conditions.  You may not impose any further
restrictions on the recipients' exercise of the rights granted herein.
You are not responsible for enforcing compliance by third parties to
this License.

  7. If, as a consequence of a court judgment or allegation of patent
infringement or for any other reason (not limited to patent issues),
conditions are imposed on you (whether by court order, agreement or
otherwise) that contradict the conditions of this License, they do not
excuse you from the conditions of this License.  If you cannot
distribute so as to satisfy simultaneously your obligations under this
License and any other pertinent obligations, then as a consequence you
may not distribute the Program at all.  For example, if a patent
license would not permit royalty-free redistribution of the Program by
all those who receive copies directly or indirectly through you, then
the only way you could satisfy both it and this License would be to
refrain entirely from distribution of the Program.

If any portion of this section is held invalid or unenforceable under
any particular circumstance, the balance of the section is intended to
apply and the section as a whole is intended to apply in other
circumstances.

It is not the purpose of this section to induce you to infringe any
patents or other property right claims or to contest validity of any
such claims; this section has the sole purpose of protecting the
integrity of the free software distribution system, which is
implemented by public license practices.  Many people have made
generous contributions to the wide range of software distributed
through that system in reliance on consistent application of that
system; it is up to the author/donor to decide if he or she is willing
to distribute software through any other system and a licensee cannot
impose that choice.

This section is intended to make thoroughly clear what is believed to
be a consequence of the rest of this License.

  8. If the distribution and/or use of the Program is restricted in
certain countries either by patents or by copyrighted interfaces, the
original copyright holder who places the Program under this License
may add an explicit geographical distribution limitation excluding
those countries, so that distribution is permitted only in or among
countries not thus excluded.  In such case, this License incorporates
the limitation as if written in the body of this License.

  9. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of this License which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
this License, you may choose any version ever published by the Free Software
Foundation.

  10. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

			    NO WARRANTY

  11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

		     END OF TERMS AND CONDITIONS

	    How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to the public, the best way to achieve this is to make it
free software which everyone can redistribute and change under these terms.

  To do so, attach the following notices to the program.  It is safest
to attach them to the start of each source file to most effectively
convey the exclusion of warranty; and each file should have at least
the "copyright" line and a pointer to where the full notice is found.

    &lt;one line to give the program's name and a brief idea of what it does.&gt;
    Copyright (C) &lt;year&gt;  &lt;name of author&gt;

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


Also add information on how to contact you by electronic and paper mail.

If the program is interactive, make it output a short notice like this
when it starts in an interactive mode:

    Gnomovision version 69, Copyright (C) year name of author
    Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the appropriate
parts of the General Public License.  Of course, the commands you use may
be called something other than `show w' and `show c'; they could even be
mouse-clicks or menu items--whatever suits your program.

You should also get your employer (if you work as a programmer) or your
school, if any, to sign a "copyright disclaimer" for the program, if
necessary.  Here is a sample; alter the names:

  Yoyodyne, Inc., hereby disclaims all copyright interest in the program
  `Gnomovision' (which makes passes at compilers) written by James Hacker.

  &lt;signature of Ty Coon&gt;, 1 April 1989
  Ty Coon, President of Vice

This General Public License does not permit incorporating your program into
proprietary programs.  If your program is a subroutine library, you may
consider it more useful to permit linking proprietary applications with the
library.  If this is what you want to do, use the GNU Library General
Public License instead of this License.</pre>
\endhtmlonly
  */
