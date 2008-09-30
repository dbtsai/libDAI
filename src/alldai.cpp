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


#include <string>
#include <dai/alldai.h>
#include <dai/properties.h>
#include <dai/exceptions.h>


namespace dai {


using namespace std;


InfAlg *newInfAlg( const string &name, const FactorGraph &fg, const PropertySet &opts ) {
    if( name == ExactInf::Name )
        return new ExactInf (fg, opts);
#ifdef DAI_WITH_BP
    if( name == BP::Name ) 
        return new BP (fg, opts);
#endif
#ifdef DAI_WITH_MF
    if( name == MF::Name ) 
        return new MF (fg, opts);
#endif
#ifdef DAI_WITH_HAK
    if( name == HAK::Name ) 
        return new HAK (fg, opts);
#endif
#ifdef DAI_WITH_LC
    if( name == LC::Name )
        return new LC (fg, opts);
#endif
#ifdef DAI_WITH_TREEEP
    if( name == TreeEP::Name )
        return new TreeEP (fg, opts);
#endif
#ifdef DAI_WITH_JTREE
    if( name == JTree::Name )
        return new JTree (fg, opts);
#endif
#ifdef DAI_WITH_MR
    if( name == MR::Name )
        return new MR (fg, opts);
#endif
    DAI_THROW(UNKNOWN_DAI_ALGORITHM);
}


} // end of namespace dai
