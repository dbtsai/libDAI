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
WITH_BP = true
WITH_MF = true
WITH_HAK = true
WITH_LC = true
WITH_TREEEP = true
WITH_JTREE = true
WITH_MR = true
DEBUG = true
NEW_MATLAB = true
WITH_MATLAB =

# Directories
INC = include/dai
SRC = src
LIB = lib

# We use the BOOST Program Options library
BOOSTFLAGS = -lboost_program_options

# Compile using GNU C++ Compiler
CC = g++

# Flags for the C++ compiler
CCFLAGS = -Wall -W -Wextra -fpic -I./include -Llib -O3 #-pg #-static -DVERBOSE
ifdef DEBUG
CCFLAGS := $(CCFLAGS) -g -DDAI_DEBUG
else
CCFLAGS := $(CCFLAGS)
endif

ifdef WITH_BP
CCFLAGS := $(CCFLAGS) -DWITH_BP
OBJECTS := $(OBJECTS) bp.o
endif
ifdef WITH_MF
CCFLAGS := $(CCFLAGS) -DWITH_MF
OBJECTS := $(OBJECTS) mf.o
endif
ifdef WITH_HAK
CCFLAGS := $(CCFLAGS) -DWITH_HAK
OBJECTS := $(OBJECTS) hak.o
endif
ifdef WITH_LC
CCFLAGS := $(CCFLAGS) -DWITH_LC
OBJECTS := $(OBJECTS) lc.o
endif
ifdef WITH_TREEEP
CCFLAGS := $(CCFLAGS) -DWITH_TREEEP
OBJECTS := $(OBJECTS) treeep.o
endif
ifdef WITH_JTREE
CCFLAGS := $(CCFLAGS) -DWITH_JTREE
OBJECTS := $(OBJECTS) jtree.o
endif
ifdef WITH_MR
CCFLAGS := $(CCFLAGS) -DWITH_MR
OBJECTS := $(OBJECTS) mr.o
endif

ifdef WITH_MATLAB
# Replace the following by the directory where Matlab has been installed
MATLABDIR = /opt/matlab/bin
# Replace the following with the extension of compiled MEX files on this platform, e.g. .mexglx for x86
MEXEXT = .mexglx
MEX = $(MATLABDIR)/mex
MEXFLAGS = -I.
ifdef DEBUG
MEXFLAGS := $(MEXFLAGS) -g -DDAI_DEBUG
endif
ifdef NEW_MATLAB
MEXFLAGS := $(MEXFLAGS) -largeArrayDims
else
MEXFLAGS := $(MEXFLAGS) -DSMALLMEM
endif
endif

HEADERS = $(INC)/bipgraph.h $(INC)/diffs.h $(INC)/index.h $(INC)/var.h $(INC)/factor.h $(INC)/varset.h $(INC)/prob.h $(INC)/daialg.h $(INC)/properties.h $(INC)/alldai.h $(INC)/enum.h $(INC)/x2x.h

TARGETS = tests utils $(LIB)/libdai.a example testregression doc
ifdef WITH_MATLAB
TARGETS := $(TARGETS) matlabs
endif
all : $(TARGETS)
	echo -e "\a"

matlabs : matlab/dai.$(MEXEXT) matlab/dai_readfg.$(MEXEXT) matlab/dai_writefg.$(MEXEXT) matlab/dai_removeshortloops.$(MEXEXT) matlab/dai_potstrength.$(MEXEXT)

$(LIB)/libdai.a : daialg.o alldai.o clustergraph.o factorgraph.o properties.o regiongraph.o util.o weightedgraph.o x2x.o $(OBJECTS)
	ar rcs $(LIB)/libdai.a daialg.o alldai.o clustergraph.o factorgraph.o properties.o regiongraph.o util.o weightedgraph.o x2x.o $(OBJECTS)

tests : tests/test

utils : utils/createfg utils/fg2dot utils/remove_short_loops utils/fginfo

testregression : tests/test
	echo Testing...this can take a while...
	cd tests; time ./testregression; cd ..

doc : $(INC)/*.h $(SRC)/*.cpp doxygen.conf
	doxygen doxygen.conf

clean :
	rm *.o example matlab/*.$(MEXEXT) matlab/*.o tests/test utils/fg2dot utils/createfg utils/remove_short_loops utils/fginfo $(LIB)/libdai.a; echo
	rm -R doc; echo


daialg.o : $(SRC)/daialg.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/daialg.cpp

bp.o : $(SRC)/bp.cpp $(INC)/bp.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/bp.cpp

lc.o : $(SRC)/lc.cpp $(INC)/lc.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/lc.cpp

mf.o : $(SRC)/mf.cpp $(INC)/mf.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/mf.cpp

factorgraph.o : $(SRC)/factorgraph.cpp $(INC)/factorgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/factorgraph.cpp

util.o : $(SRC)/util.cpp $(INC)/util.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/util.cpp

regiongraph.o : $(SRC)/regiongraph.cpp $(INC)/regiongraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/regiongraph.cpp

hak.o : $(SRC)/hak.cpp $(INC)/hak.h $(HEADERS) $(INC)/regiongraph.h
	$(CC) $(CCFLAGS) -c $(SRC)/hak.cpp

clustergraph.o : $(SRC)/clustergraph.cpp $(INC)/clustergraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/clustergraph.cpp

jtree.o : $(SRC)/jtree.cpp $(INC)/jtree.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h
	$(CC) $(CCFLAGS) -c $(SRC)/jtree.cpp

treeep.o : $(SRC)/treeep.cpp $(INC)/treeep.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h $(INC)/jtree.h
	$(CC) $(CCFLAGS) -c $(SRC)/treeep.cpp

weightedgraph.o : $(SRC)/weightedgraph.cpp $(INC)/weightedgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/weightedgraph.cpp

mr.o : $(SRC)/mr.cpp $(INC)/mr.h $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/mr.cpp

properties.o : $(SRC)/properties.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/properties.cpp

alldai.o : $(SRC)/alldai.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/alldai.cpp

x2x.o : $(SRC)/x2x.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c $(SRC)/x2x.cpp

# EXAMPLE
##########

example : example.cpp $(HEADERS) $(LIB)/libdai.a
	$(CC) $(CCFLAGS) -o example example.cpp -ldai

# TESTS
########

tests/test : tests/test.cpp $(HEADERS) $(LIB)/libdai.a
	$(CC) $(CCFLAGS) -o tests/test tests/test.cpp -ldai $(BOOSTFLAGS)


# MATLAB INTERFACE
###################

matlab/dai.$(MEXEXT) : matlab/dai.cpp $(HEADERS) matlab/matlab.o $(LIB)/libdai.a
	$(MEX) $(MEXFLAGS) -o matlab/dai matlab/dai.cpp matlab/matlab.o $(LIB)/libdai.a

matlab/dai_readfg.$(MEXEXT) : matlab/dai_readfg.cpp $(HEADERS) factorgraph.o matlab/matlab.o
	$(MEX) $(MEXFLAGS) -o matlab/dai_readfg matlab/dai_readfg.cpp factorgraph.o matlab/matlab.o

matlab/dai_writefg.$(MEXEXT) : matlab/dai_writefg.cpp $(HEADERS) factorgraph.o matlab/matlab.o
	$(MEX) $(MEXFLAGS) -o matlab/dai_writefg matlab/dai_writefg.cpp factorgraph.o matlab/matlab.o

matlab/dai_removeshortloops.$(MEXEXT) : matlab/dai_removeshortloops.cpp $(HEADERS) factorgraph.o matlab/matlab.o
	$(MEX) $(MEXFLAGS) -o matlab/dai_removeshortloops matlab/dai_removeshortloops.cpp factorgraph.o matlab/matlab.o

matlab/dai_potstrength.$(MEXEXT) : matlab/dai_potstrength.cpp $(HEADERS) matlab/matlab.o
	$(MEX) $(MEXFLAGS) -o matlab/dai_potstrength matlab/dai_potstrength.cpp matlab/matlab.o

matlab/matlab.o : matlab/matlab.cpp matlab/matlab.h $(HEADERS)
	$(MEX) $(MEXFLAGS) -outdir matlab -c matlab/matlab.cpp


# UTILS
########

utils/createfg : utils/createfg.cpp $(HEADERS) factorgraph.o weightedgraph.o util.o
	$(CC) $(CCFLAGS) -o utils/createfg utils/createfg.cpp factorgraph.o weightedgraph.o util.o $(BOOSTFLAGS)

utils/fg2dot : utils/fg2dot.cpp $(HEADERS) factorgraph.o
	$(CC) $(CCFLAGS) -o utils/fg2dot utils/fg2dot.cpp factorgraph.o

utils/remove_short_loops : utils/remove_short_loops.cpp $(HEADERS) factorgraph.o
	$(CC) $(CCFLAGS) -o utils/remove_short_loops utils/remove_short_loops.cpp factorgraph.o

utils/fginfo : utils/fginfo.cpp $(HEADERS) factorgraph.o
	$(CC) $(CCFLAGS) -o utils/fginfo utils/fginfo.cpp factorgraph.o
