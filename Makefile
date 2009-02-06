# Copyright (C) 2006-2009  Joris Mooij  [joris dot mooij at tuebingen dot mpg dot de]
# Radboud University Nijmegen, The Netherlands /
# Max Planck Institute for Biological Cybernetics, Germany
#   
# This file is part of libDAI.
#
# libDAI is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# libDAI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libDAI; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


# Choose OS from {LINUX, WINDOWS, CYGWIN, MACOSX}
# LINUX:   GNU/Linux and other UNIX variants
# WINDOWS: Visual C++ with GNU Make
# CYGWIN:  CygWin
# MACOSX:  MacOSX
OS=LINUX

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

# Directories
#   Location libDAI headers
INC=include/dai
#   Location of libDAI source files
SRC=src
#   Destination directory of libDAI library
LIB=lib
#   Additional iclude paths for C compiler
CCINC=-Iinclude

# Extensions (library, object, executable, matlab compiled MEX file)
ifneq ($(OS),WINDOWS)
  LE=.a
  OE=.o
  EE=
  ME=.mexglx
else
  LE=.lib
  OE=.obj
  EE=.exe
  ME=.mexglx
endif

# Libraries
ifneq ($(OS),WINDOWS)
  LIBS=-ldai
  # Additional library paths for linker
  CCLIB=-Llib
else
  # For some reason, we have to add the VC library path, although it is in the environment
  LIBS=/link $(LIB)/libdai$(LE) /LIBPATH:"C:\Program Files\Microsoft Visual Studio 9.0\VC\ATLMFC\LIB" /LIBPATH:"C:\Program Files\Microsoft Visual Studio 9.0\VC\LIB" /LIBPATH:"C:\Program Files\Microsoft SDKs\Windows\v6.0A\lib"
endif

# Tell the linker to link with the BOOST Program Options library
ifeq ($(OS),CYGWIN)
  BOOSTLIBS=-lboost_program_options-gcc34-mt
endif
ifeq ($(OS),MACOSX)
  BOOSTLIBS=-lboost_program_options-mt
endif
ifeq ($(OS),LINUX)
  BOOSTLIBS=-lboost_program_options
endif
ifeq ($(OS),WINDOWS)
  BOOSTLIBS=/LIBPATH:C:\boost_1_36_0\stage\lib
endif

# Compiler specific options
ifneq ($(OS),WINDOWS)
  # Compile using GNU C++ Compiler
  CC=g++
  # Output filename option of the compiler
  CCO=-o 
else
  # Compile using Visual C++ Compiler
  CC=cl
  # Output filename option
  CCO=/Fe
endif

# Flags for the C++ compiler
ifneq ($(OS),WINDOWS)
  CCFLAGS=-O3 -Wno-deprecated -Wall -W -Wextra -fpic
  CCDEBUGFLAGS=-g -DDAI_DEBUG
else
  CCFLAGS=/Iinclude /IC:\boost_1_36_0 /EHsc /Ox
  CCDEBUGFLAGS=/Zi -DDAI_DEBUG
endif

ifeq ($(OS),CYGWIN)
  CCINC:=$(CCINC) -I/usr/local/include/boost-1_37
  # dynamic linking of Boost libraries seems not to work on Cygwin
  CCFLAGS:=$(CCFLAGS) -DCYGWIN -static
endif
ifeq ($(OS),MACOSX)
  # indicate where your boost headers and libraries are (likely where macports installs libraries)
  CCINC:=$(CCINC) -I/opt/local/include
  CCLIB:=$(CCLIB) -L/opt/local/lib
endif

# Build targets
TARGETS=tests utils $(LIB)/libdai$(LE) examples testregression
ifneq ($(OS),WINDOWS)
  TARGETS:=$(TARGETS) doc 
endif

ifdef WITH_MATLAB
  ifneq ($(OS),WINDOWS)
    # Replace the following by the directory where Matlab has been installed
    MATLABDIR=/agbs/share/sw/matlab
    MEX=$(MATLABDIR)/bin/mex
    MEXFLAGS=CXX\#$(CC) CXXFLAGS\#'$(CCFLAGS)'
  else
    # Replace the following by the directory where Matlab has been installed
    MATLABDIR=c:\matlab
    MEX=$(MATLABDIR)\bin\mex
    MEXFLAGS=-Iinclude CXX\#$(CC) CXXFLAGS\#"/EHsc /Ox"
  endif
endif


ifdef DEBUG
  CCFLAGS:=$(CCFLAGS) $(CCDEBUGFLAGS)
endif
ifeq ($(OS),WINDOWS)
  CCFLAGS:=$(CCFLAGS) -DWINDOWS
endif
ifdef WITH_MATLAB
  TARGETS:=$(TARGETS) matlabs
  ifdef NEW_MATLAB
    MEXFLAGS:=$(MEXFLAGS) -largeArrayDims
  else
    MEXFLAGS:=$(MEXFLAGS) -DSMALLMEM
  endif
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
ifdef WITH_GIBBS
  CCFLAGS:=$(CCFLAGS) -DDAI_WITH_GIBBS
  OBJECTS:=$(OBJECTS) gibbs$(OE)
endif


HEADERS=$(INC)/bipgraph.h $(INC)/index.h $(INC)/var.h $(INC)/factor.h $(INC)/varset.h $(INC)/smallset.h $(INC)/prob.h $(INC)/daialg.h $(INC)/properties.h $(INC)/alldai.h $(INC)/enum.h $(INC)/exceptions.h


CC:=$(CC) $(CCINC) $(CCLIB) $(CCFLAGS)
MEX:=$(MEX) $(CCLIB) $(CCINC) $(MEXFLAGS)


# META TARGETS
###############

all : $(TARGETS)

examples : examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE) examples/example_sprinkler$(EE)

matlabs : matlab/dai$(ME) matlab/dai_readfg$(ME) matlab/dai_writefg$(ME) matlab/dai_potstrength$(ME)

tests : tests/testdai$(EE)

utils : utils/createfg$(EE) utils/fg2dot$(EE) utils/fginfo$(EE)


# OBJECTS
##########

bipgraph$(OE) : $(SRC)/bipgraph.cpp $(HEADERS)
	$(CC) -c $(SRC)/bipgraph.cpp

daialg$(OE) : $(SRC)/daialg.cpp $(HEADERS)
	$(CC) -c $(SRC)/daialg.cpp

exactinf$(OE) : $(SRC)/exactinf.cpp $(INC)/exactinf.h $(HEADERS)
	$(CC) -c $(SRC)/exactinf.cpp

bp$(OE) : $(SRC)/bp.cpp $(INC)/bp.h $(HEADERS)
	$(CC) -c $(SRC)/bp.cpp

lc$(OE) : $(SRC)/lc.cpp $(INC)/lc.h $(HEADERS)
	$(CC) -c $(SRC)/lc.cpp

mf$(OE) : $(SRC)/mf.cpp $(INC)/mf.h $(HEADERS)
	$(CC) -c $(SRC)/mf.cpp

factorgraph$(OE) : $(SRC)/factorgraph.cpp $(INC)/factorgraph.h $(HEADERS)
	$(CC) -c $(SRC)/factorgraph.cpp

util$(OE) : $(SRC)/util.cpp $(INC)/util.h $(HEADERS)
	$(CC) -c $(SRC)/util.cpp

regiongraph$(OE) : $(SRC)/regiongraph.cpp $(INC)/regiongraph.h $(HEADERS)
	$(CC) -c $(SRC)/regiongraph.cpp

hak$(OE) : $(SRC)/hak.cpp $(INC)/hak.h $(HEADERS) $(INC)/regiongraph.h
	$(CC) -c $(SRC)/hak.cpp

clustergraph$(OE) : $(SRC)/clustergraph.cpp $(INC)/clustergraph.h $(HEADERS)
	$(CC) -c $(SRC)/clustergraph.cpp

jtree$(OE) : $(SRC)/jtree.cpp $(INC)/jtree.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h
	$(CC) -c $(SRC)/jtree.cpp

treeep$(OE) : $(SRC)/treeep.cpp $(INC)/treeep.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h $(INC)/jtree.h
	$(CC) -c $(SRC)/treeep.cpp

weightedgraph$(OE) : $(SRC)/weightedgraph.cpp $(INC)/weightedgraph.h $(HEADERS)
	$(CC) -c $(SRC)/weightedgraph.cpp

mr$(OE) : $(SRC)/mr.cpp $(INC)/mr.h $(HEADERS)
	$(CC) -c $(SRC)/mr.cpp

gibbs$(OE) : $(SRC)/gibbs.cpp $(INC)/gibbs.h $(HEADERS)
	$(CC) -c $(SRC)/gibbs.cpp

properties$(OE) : $(SRC)/properties.cpp $(HEADERS)
	$(CC) -c $(SRC)/properties.cpp

exceptions$(OE) : $(SRC)/exceptions.cpp $(HEADERS)
	$(CC) -c $(SRC)/exceptions.cpp

alldai$(OE) : $(SRC)/alldai.cpp $(HEADERS)
	$(CC) -c $(SRC)/alldai.cpp


# EXAMPLES
###########

examples/example$(EE) : examples/example.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)examples/example$(EE) examples/example.cpp $(LIBS)

examples/example_bipgraph$(EE) : examples/example_bipgraph.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)examples/example_bipgraph$(EE) examples/example_bipgraph.cpp $(LIBS)

examples/example_varset$(EE) : examples/example_varset.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)examples/example_varset$(EE) examples/example_varset.cpp $(LIBS)

examples/example_sprinkler$(EE) : examples/example_sprinkler.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)examples/example_sprinkler$(EE) examples/example_sprinkler.cpp $(LIBS)


# TESTS
########

tests/testdai$(EE) : tests/testdai.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)tests/testdai$(EE) tests/testdai.cpp $(LIBS) $(BOOSTLIBS)


# MATLAB INTERFACE
###################

matlab/dai$(ME) : $(SRC)/matlab/dai.cpp $(HEADERS) matlab$(OE) $(LIB)/libdai$(LE)
	$(MEX) -o matlab/dai $(SRC)/matlab/dai.cpp matlab$(OE) $(LIB)/libdai$(LE)

matlab/dai_readfg$(ME) : $(SRC)/matlab/dai_readfg.cpp $(HEADERS) factorgraph$(OE) matlab$(OE) exceptions$(OE)
	$(MEX) -o matlab/dai_readfg $(SRC)/matlab/dai_readfg.cpp factorgraph$(OE) matlab$(OE) exceptions$(OE)

matlab/dai_writefg$(ME) : $(SRC)/matlab/dai_writefg.cpp $(HEADERS) factorgraph$(OE) matlab$(OE) exceptions$(OE)
	$(MEX) -o matlab/dai_writefg $(SRC)/matlab/dai_writefg.cpp factorgraph$(OE) matlab$(OE) exceptions$(OE)

matlab/dai_potstrength$(ME) : $(SRC)/matlab/dai_potstrength.cpp $(HEADERS) matlab$(OE) exceptions$(OE)
	$(MEX) -o matlab/dai_potstrength $(SRC)/matlab/dai_potstrength.cpp matlab$(OE) exceptions$(OE)

matlab$(OE) : $(SRC)/matlab/matlab.cpp $(INC)/matlab/matlab.h $(HEADERS)
	$(MEX) -c $(SRC)/matlab/matlab.cpp


# UTILS
########

utils/createfg$(EE) : utils/createfg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)utils/createfg$(EE) utils/createfg.cpp $(LIBS) $(BOOSTLIBS)

utils/fg2dot$(EE) : utils/fg2dot.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)utils/fg2dot$(EE) utils/fg2dot.cpp $(LIBS)

utils/fginfo$(EE) : utils/fginfo.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)utils/fginfo$(EE) utils/fginfo.cpp $(LIBS)


# LIBRARY
##########

ifneq ($(OS),WINDOWS)
$(LIB)/libdai$(LE) : bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
	-mkdir -p lib
	ar rcus $(LIB)/libdai$(LE) bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
else
$(LIB)/libdai$(LE) : bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
	-mkdir lib
	lib /out:$(LIB)/libdai$(LE) bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS)
endif


# REGRESSION TESTS
###################

ifneq ($(OS),WINDOWS)
testregression : tests/testdai$(EE)
	@echo Starting regression test...this can take a minute or so!
	cd tests && ./testregression && cd ..
else
testregression : tests/testdai$(EE)
	@echo Starting regression test...this can take a minute or so!
	cd tests && testregression.bat && cd ..
endif


# DOCUMENTATION
################

doc : $(INC)/*.h $(SRC)/*.cpp examples/*.cpp doxygen.conf
	doxygen doxygen.conf


# CLEAN
########

ifneq ($(OS),WINDOWS)
.PHONY : clean
clean :
	-rm *$(OE)
	-rm matlab/*$(ME)
	-rm examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE) examples/example_sprinkler$(EE)
	-rm tests/testdai$(EE)
	-rm utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE)
	-rm -R doc
	-rm -R lib
else
.PHONY : clean
clean :
	-del *$(OE) *.ilk *.pdb *$(EE) matlab\*$(ME) examples\*$(EE) examples\*.ilk examples\*.pdb tests\testdai$(EE) tests\*.pdb tests\*.ilk utils\*$(EE) utils\*.pdb utils\*.ilk $(LIB)\libdai$(LE)
endif
