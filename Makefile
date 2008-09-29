#   Copyright (C) 2006-2008  Joris Mooij  [j dot mooij at science dot ru dot nl]
#   Radboud University Nijmegen, The Netherlands
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

# Extensions (library, object, executable extensions)
LE=.a
OE=.o
EE=

# Libraries
LIBS=-ldai

# We use the BOOST Program Options library
BOOSTLIBS=-lboost_program_options

# Compile using GNU C++ Compiler
CC=g++

# Flags for the C++ compiler
CCFLAGS=-Wno-deprecated -Wall -W -Wextra -fpic -I./include -Llib -O3 #-pg #-static #-DVERBOSE
ifdef DEBUG
CCFLAGS:=$(CCFLAGS) -g -DDAI_DEBUG
endif

ifdef WINDOWS
CCFLAGS=$(CCFLAGS) -DWINDOWS
endif

OBJECTS:=exactinf$(OE)
ifdef WITH_BP
CCFLAGS:=$(CCFLAGS) -DWITH_BP
OBJECTS:=$(OBJECTS) bp$(OE)
endif
ifdef WITH_MF
CCFLAGS:=$(CCFLAGS) -DWITH_MF
OBJECTS:=$(OBJECTS) mf$(OE)
endif
ifdef WITH_HAK
CCFLAGS:=$(CCFLAGS) -DWITH_HAK
OBJECTS:=$(OBJECTS) hak$(OE)
endif
ifdef WITH_LC
CCFLAGS:=$(CCFLAGS) -DWITH_LC
OBJECTS:=$(OBJECTS) lc$(OE)
endif
ifdef WITH_TREEEP
CCFLAGS:=$(CCFLAGS) -DWITH_TREEEP
OBJECTS:=$(OBJECTS) treeep$(OE)
endif
ifdef WITH_JTREE
CCFLAGS:=$(CCFLAGS) -DWITH_JTREE
OBJECTS:=$(OBJECTS) jtree$(OE)
endif
ifdef WITH_MR
CCFLAGS:=$(CCFLAGS) -DWITH_MR
OBJECTS:=$(OBJECTS) mr$(OE)
endif

ifdef WITH_MATLAB
# Replace the following by the directory where Matlab has been installed
MATLABDIR=/opt/matlab/bin
# Replace the following with the extension of compiled MEX files on this platform, e.g. .mexglx for x86
ME=.mexglx
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

HEADERS=$(INC)/bipgraph.h $(INC)/diffs.h $(INC)/index.h $(INC)/var.h $(INC)/factor.h $(INC)/varset.h $(INC)/prob.h $(INC)/daialg.h $(INC)/properties.h $(INC)/alldai.h $(INC)/enum.h $(INC)/x2x.h $(INC)/exceptions.h

TARGETS=tests utils $(LIB)/libdai$(LE) example$(EE) testregression doc
ifdef WITH_MATLAB
TARGETS:=$(TARGETS) matlabs
endif
all : $(TARGETS)
	echo -e "\a"

matlabs : matlab/dai.$(ME) matlab/dai_readfg.$(ME) matlab/dai_writefg.$(ME) matlab/dai_potstrength.$(ME)

$(LIB)/libdai$(LE) : bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) x2x$(OE) exceptions$(OE) $(OBJECTS)
	ar rcus $(LIB)/libdai$(LE) bipgraph$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) x2x$(OE) exceptions$(OE) $(OBJECTS)

tests : tests/test$(EE)

utils : utils/createfg$(EE) utils/fg2dot$(EE) utils/fginfo$(EE)

testregression : tests/test
	echo Testing...this can take a while...
	cd tests; time ./testregression; cd ..

doc : $(INC)/*.h $(SRC)/*.cpp doxygen.conf
	doxygen doxygen.conf

clean :
	rm *$(OE) example$(EE) matlab/*.$(ME) matlab/*$(OE) tests/test$(EE) utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE) $(LIB)/libdai$(LE); echo
	rm -R doc; echo

bipgraph$(OE) : $(SRC)/bipgraph.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/bipgraph.cpp

daialg$(OE) : $(SRC)/daialg.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/daialg.cpp

exactinf$(OE) : $(SRC)/exactinf.cpp $(INC)/exactinf.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/exactinf.cpp

bp$(OE) : $(SRC)/bp.cpp $(INC)/bp.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/bp.cpp

lc$(OE) : $(SRC)/lc.cpp $(INC)/lc.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/lc.cpp

mf$(OE) : $(SRC)/mf.cpp $(INC)/mf.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/mf.cpp

factorgraph$(OE) : $(SRC)/factorgraph.cpp $(INC)/factorgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/factorgraph.cpp

util$(OE) : $(SRC)/util.cpp $(INC)/util.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/util.cpp

regiongraph$(OE) : $(SRC)/regiongraph.cpp $(INC)/regiongraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/regiongraph.cpp

hak$(OE) : $(SRC)/hak.cpp $(INC)/hak.h $(HEADERS) $(INC)/regiongraph.h
	$(CC) $(CCFLAGS) -c $(SRC)/hak.cpp

clustergraph$(OE) : $(SRC)/clustergraph.cpp $(INC)/clustergraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/clustergraph.cpp

jtree$(OE) : $(SRC)/jtree.cpp $(INC)/jtree.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h
	$(CC) $(CCFLAGS) -c $(SRC)/jtree.cpp

treeep$(OE) : $(SRC)/treeep.cpp $(INC)/treeep.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h $(INC)/jtree.h
	$(CC) $(CCFLAGS) -c $(SRC)/treeep.cpp

weightedgraph$(OE) : $(SRC)/weightedgraph.cpp $(INC)/weightedgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/weightedgraph.cpp

mr$(OE) : $(SRC)/mr.cpp $(INC)/mr.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/mr.cpp

properties$(OE) : $(SRC)/properties.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/properties.cpp

exceptions$(OE) : $(SRC)/exceptions.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/exceptions.cpp

alldai$(OE) : $(SRC)/alldai.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/alldai.cpp

x2x$(OE) : $(SRC)/x2x.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/x2x.cpp


# EXAMPLE
##########

example$(EE) : example.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCFLAGS) -o example$(EE) example.cpp $(LIBS)

# TESTS
########

tests/test$(EE) : tests/test.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCFLAGS) -o tests/test$(EE) tests/test.cpp $(LIBS) $(BOOSTLIBS)


# MATLAB INTERFACE
###################

matlab/dai.$(ME) : matlab/dai.cpp $(HEADERS) matlab/matlab$(OE) $(LIB)/libdai$(LE)
	$(MEX) $(MEXFLAGS) -o matlab/dai matlab/dai.cpp matlab/matlab$(OE) $(LIB)/libdai$(LE)

matlab/dai_readfg.$(ME) : matlab/dai_readfg.cpp $(HEADERS) factorgraph$(OE) matlab/matlab$(OE) exceptions$(OE)
	$(MEX) $(MEXFLAGS) -o matlab/dai_readfg matlab/dai_readfg.cpp factorgraph$(OE) matlab/matlab$(OE) exceptions$(OE)

matlab/dai_writefg.$(ME) : matlab/dai_writefg.cpp $(HEADERS) factorgraph$(OE) matlab/matlab$(OE) exceptions$(OE)
	$(MEX) $(MEXFLAGS) -o matlab/dai_writefg matlab/dai_writefg.cpp factorgraph$(OE) matlab/matlab$(OE) exceptions$(OE)

matlab/dai_potstrength.$(ME) : matlab/dai_potstrength.cpp $(HEADERS) matlab/matlab$(OE) exceptions$(OE)
	$(MEX) $(MEXFLAGS) -o matlab/dai_potstrength matlab/dai_potstrength.cpp matlab/matlab$(OE) exceptions$(OE)

matlab/matlab$(OE) : matlab/matlab.cpp matlab/matlab.h $(HEADERS)
	$(MEX) $(MEXFLAGS) -outdir matlab -c matlab/matlab.cpp


# UTILS
########

utils/createfg$(EE) : utils/createfg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCFLAGS) -o utils/createfg utils/createfg.cpp $(LIBS) $(BOOSTLIBS)

utils/fg2dot$(EE) : utils/fg2dot.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCFLAGS) -o utils/fg2dot utils/fg2dot.cpp $(LIBS)

utils/fginfo$(EE) : utils/fginfo.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCFLAGS) -o utils/fginfo utils/fginfo.cpp $(LIBS)
