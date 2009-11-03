/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  libDAI is licensed under the terms of the GNU General Public License version
 *  2, or (at your option) any later version. libDAI is distributed without any
 *  warranty. See the file COPYING for more details.
 *
 *  Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at libdai dot org]
 *  Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands
 */


/// \file
/// \brief Main libDAI header file. It \#includes all other libDAI headers.


#ifndef __defined_libdai_alldai_h
#define __defined_libdai_alldai_h


#include <string>
#include <dai/daialg.h>
#include <dai/properties.h>
#include <dai/exactinf.h>
#include <dai/evidence.h>
#include <dai/emalg.h>
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
#ifdef DAI_WITH_GIBBS
    #include <dai/gibbs.h>
#endif
#ifdef DAI_WITH_CBP
    #include <dai/cbp.h>
#endif


/// Namespace for libDAI
namespace dai {


/// Constructs a new inference algorithm.
/** \param name The name of the inference algorithm (should be one of the names in DAINames).
 *  \param fg The FactorGraph that the algorithm should be applied to.
 *  \param opts A PropertySet specifying the options for the algorithm.
 *  \return Returns a pointer to the new InfAlg object; it is the responsibility of the caller to delete it later.
 */
InfAlg *newInfAlg( const std::string &name, const FactorGraph &fg, const PropertySet &opts );


/// Constructs a new inference algorithm.
/** \param nameOpts The name and options of the inference algorithm (should be in the format "name[key1=val1,key2=val2,...,keyn=valn]").
 *  \param fg The FactorGraph that the algorithm should be applied to.
 *  \return Returns a pointer to the new InfAlg object; it is the responsibility of the caller to delete it later.
 *  \todo Support aliases like in testdai
 */
InfAlg *newInfAlgFromString( const std::string &nameOpts, const FactorGraph &fg );


/// Contains the names of all inference algorithms compiled into libDAI.
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
#ifdef DAI_WITH_GIBBS
    Gibbs::Name,
#endif
#ifdef DAI_WITH_CBP
    CBP::Name,
#endif
    ""
};


} // end of namespace dai


#endif
