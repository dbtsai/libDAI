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


#ifndef __defined_libdai_alldai_h
#define __defined_libdai_alldai_h


#include <string>
#include <dai/daialg.h>
#ifdef WITH_BP
    #include <dai/bp.h>
#endif
#ifdef WITH_MF
    #include <dai/mf.h>
#endif
#ifdef WITH_HAK
    #include <dai/hak.h>
#endif
#ifdef WITH_LC
    #include <dai/lc.h>
#endif
#ifdef WITH_TREEEP
    #include <dai/treeep.h>
#endif
#ifdef WITH_JTREE
    #include <dai/jtree.h>
#endif
#ifdef WITH_MR
    #include <dai/mr.h>
#endif


namespace dai {


/// newInfAlg constructs a new approximate inference algorithm named name for the
/// FactorGraph fg with optionts opts and returns a pointer to the new object.
/// The caller needs to delete it (maybe some sort of smart_ptr might be useful here).
InfAlg *newInfAlg( const std::string &name, const FactorGraph &fg, const Properties &opts );


/// DAINames contains the names of all approximate inference algorithms

static const char* DAINames[] = {
#ifdef WITH_BP
    BP::Name, 
#endif
#ifdef WITH_MF
    MF::Name,
#endif
#ifdef WITH_HAK
    HAK::Name,
#endif
#ifdef WITH_LC
    LC::Name,
#endif
#ifdef WITH_TREEEP
    TreeEP::Name,
#endif
#ifdef WITH_JTREE
    JTree::Name,
#endif
#ifdef WITH_MR
    MR::Name,
#endif
    ""
};


} // end of namespace dai


#endif
