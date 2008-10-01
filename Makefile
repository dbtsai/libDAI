#   Copyright (C) 2006-2008  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
#   Radboud University Nijmegen, The Netherlands /
#   Max Planck Institute for Biological Cybernetics, Germany
#   
#   This file is part of libDAI.
#
#   libDAI is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   libDAI is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with libDAI; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


# Enable/disable various approximate inference methods
WITH_BP=true
WITH_MF=true
WITH_HAK=true
WITH_LC=true
WITH_TREEEP=true
WITH_JTREE=true
WITH_MR=true
# Build with debug info?
DEBUG=true
# Build matlab interface?
WITH_MATLAB=
# New/old matlab version?
NEW_MATLAB=true
# Windows or linux (default)?
WINDOWS=

# Directories
INC=include/dai
SRC=src
LIB=lib

# Extensions (library, object, executable, MEX file extensions)
LE=.a
OE=.o
EE=
ME=.mexglx

# Libraries
LIBS=-ldai

# We use the BOOST Program Options library
BOOSTLIBS=-lboost_program_options

# Compile using GNU C++ Compiler
CC=g++
# Output filename option
CCO=-o

# Flags for the C++ compiler
CCFLAGS=-Wno-deprecated -Wall -W -Wextra -fpic -I./include -Llib -O3 #-pg
ifdef DEBUG
CCFLAGS:=$(CCFLAGS) -g -DDAI_DEBUG
endif

ifdef WINDOWS
CCFLAGS:=$(CCFLAGS) -DWINDOWS
endif

OBJECTS:=exactinf$(OE)
ifdef WITH_BP
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_BP
OBJECTS:=$(OBJECTS) bp$(OE)
endif
ifdef WITH_MF
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_MF
OBJECTS:=$(OBJECTS) mf$(OE)
endif
ifdef WITH_HAK
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_HAK
OBJECTS:=$(OBJECTS) hak$(OE)
endif
ifdef WITH_LC
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_LC
OBJECTS:=$(OBJECTS) lc$(OE)
endif
ifdef WITH_TREEEP
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_TREEEP
OBJECTS:=$(OBJECTS) treeep$(OE)
endif
ifdef WITH_JTREE
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_JTREE
OBJECTS:=$(OBJECTS) jtree$(OE)
endif
ifdef WITH_MR
CCFLAGS:=$(CCFLAGS) -DDAI_WITH_MR
OBJECTS:=$(OBJECTS) mr$(OE)
endif

ifdef WITH_MATLAB
# Replace the following by the directory where Matlab has been installed
MATLABDIR=/opt/matlab/bin
# Replace the following with the extension of compiled MEX files on this platform, e.g. .mexglx for x86
MEX=$(MATLABDIR)/mex
MEXFLAGS=-I.
ifdef DEBUG
MEXFLAGS:=$(MEXFLAGS) -g -DDAI_DEBUG
endif
ifdef NEW_MATLAB
MEXFLAGS:=$(MEXFLAGS) -largeArrayDims
else
MEXFLAGS:=$(MEXFLAGS) -DSMALLMEM
endif
endif

HEADERS=$(INC)/bipgraph.h $(INC)/index.h $(INC)/var.h $(INC)/factor.h $(INC)/varset.h $(INC)/smallset.h $(INC)/prob.h $(INC)/daialg.h $(INC)/properties.h $(INC)/alldai.h $(INC)/enum.h $(INC)/exceptions.h

TARGETS=tests utils $(LIB)/libdai$(LE) example$(EE) testregression doc
ifdef WITH_MATLAB
TARGETS:=$(TARGETS) matlabs
endif
all : $(TARGETS)
	echo -e "\a"

matlabs : matlab/dai$(ME) matlab/dai_readfg$(ME) matlab/dai_writefg$(ME) matlab/dai_potstrength$(ME)

$(LIB)/libdai$(LE) : bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
	-mkdir -p lib
	ar rcus $(LIB)/libdai$(LE) bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)

tests : tests/testdai$(EE)

utils : utils/createfg$(EE) utils/fg2dot$(EE) utils/fginfo$(EE)

testregression : tests/testdai
	@echo Starting regression test...this can take a minute or so!
	cd tests; time ./testregression; cd ..

doc : $(INC)/*.h $(SRC)/*.cpp doxygen.conf
	-mkdir -p doc
	doxygen doxygen.conf

.PHONY : clean
clean :
	-rm *$(OE) 
	-rm matlab/*$(ME) matlab/*$(OE) 
	-rm example$(EE) tests/testdai$(EE) utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE)
	-rm -R doc
	-rm -R lib

include Makefile.shared
