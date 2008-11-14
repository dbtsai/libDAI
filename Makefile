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
WITH_GIBBS=true
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

# Extensions (library, object, executable, matlab compiled MEX file)
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
CCFLAGS=-O3 -Wno-deprecated -Wall -W -Wextra -fpic -Iinclude -Llib
CCDEBUGFLAGS=-g -DDAI_DEBUG

# Build targets
TARGETS=tests utils $(LIB)/libdai$(LE) examples testregression doc

ifdef WITH_MATLAB
	# Replace the following by the directory where Matlab has been installed
	MATLABDIR=/agbs/share/sw/matlab
	MEX=$(MATLABDIR)/bin/mex
	MEXFLAGS=-Iinclude CXX\#$(CC) CXXFLAGS\#'-O3 -Wno-deprecated -Wall -W -Wextra -fpic'
endif


include Makefile.shared


$(LIB)/libdai$(LE) : bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
	-mkdir -p lib
	ar rcus $(LIB)/libdai$(LE) bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)

testregression : tests/testdai$(EE)
	@echo Starting regression test...this can take a minute or so!
	cd tests && ./testregression && cd ..

doc : $(INC)/*.h $(SRC)/*.cpp examples/*.cpp doxygen.conf
	doxygen doxygen.conf

.PHONY : clean
clean :
	-rm *$(OE)
	-rm matlab/*$(ME)
	-rm examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE)
	-rm tests/testdai$(EE)
	-rm utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE)
	-rm -R doc
	-rm -R lib
