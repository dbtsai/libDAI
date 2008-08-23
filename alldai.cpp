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


#include "alldai.h"


InfAlg *newInfAlg( const string &name, const FactorGraph &fg, const Properties &opts ) {
    if( name == BP::Name ) 
        return new BP (fg, opts);
    else if( name == MF::Name ) 
        return new MF (fg, opts);
    else if( name == HAK::Name ) 
        return new HAK (fg, opts);
    else if( name == LC::Name )
        return new LC (fg, opts);
    else if( name == TreeEP::Name )
        return new TreeEP (fg, opts);
    else if( name == MR::Name )
        return new MR (fg, opts);
    else if( name == JTree::Name )
        return new JTree (fg, opts);
    else
        throw "Unknown inference algorithm";
}
