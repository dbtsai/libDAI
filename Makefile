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


# We use the BOOST Program Options library
BOOSTFLAGS = -lboost_program_options

# Compile using GNU C++ Compiler
CC = g++

# Flags for the C++ compiler
CCFLAGS = -Wall -fpic -g -DDEBUG -I. -I$(HOME)/include -O3 #-static #-pg #-DVERBOSE

# To enable the Matlab interface, define WITH_MATLAB = yes
WITH_MATLAB = 
ifdef WITH_MATLAB
# Replace the following by the directory where Matlab has been installed
MATLABDIR = /opt/matlab/bin
MEX = $(MATLABDIR)/mex
MEXFLAGS = -g -I. -DDEBUG -largeArrayDims #-g means debugging
endif

# Replace the following with the extension of compiled MEX files on this platform, e.g. .mexglx for x86
MEXEXT = .mexglx

HEADERS = bipgraph.h diffs.h index.h var.h factor.h varset.h prob.h daialg.h properties.h alldai.h enum.h x2x.h

# target matlabs is disabled by default since it only compiles with a very recent MatLab version
TARGETS = tests utils libdai.a example testregression
ifdef WITH_MATLAB
TARGETS := $(TARGETS) matlabs
endif
all : $(TARGETS)
	echo -e "\a"

matlabs : matlab/dai.$(MEXEXT) matlab/dai_readfg.$(MEXEXT) matlab/dai_writefg.$(MEXEXT) matlab/dai_removeshortloops.$(MEXEXT) matlab/dai_potstrength.$(MEXEXT)

libdai.a : daialg.o alldai.o bp.o clustergraph.o factorgraph.o hak.o jtree.o lc.o mf.o mr.o properties.o regiongraph.o util.o treeep.o weightedgraph.o x2x.o
	ar rcs libdai.a daialg.o alldai.o bp.o clustergraph.o factorgraph.o hak.o jtree.o lc.o mf.o mr.o properties.o regiongraph.o util.o treeep.o weightedgraph.o x2x.o

tests : tests/test

utils : utils/createfg utils/fg2dot utils/remove_short_loops utils/fginfo

testregression : tests/test
	echo Testing...this can take a while...
	cd tests; ./testregression; cd ..

doc : *.h *.cpp doxygen.conf
	doxygen doxygen.conf

clean :
	rm *.o *.$(MEXEXT) example matlab/*.$(MEXEXT) matlab/*.o tests/test utils/fg2dot utils/createfg utils/remove_short_loops utils/fginfo libdai.a; echo
	rm -R doc; echo


daialg.o : daialg.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c daialg.cpp

bp.o : bp.cpp bp.h $(HEADERS)
	$(CC) $(CCFLAGS) -c bp.cpp

lc.o : lc.cpp lc.h $(HEADERS)
	$(CC) $(CCFLAGS) -c lc.cpp

mf.o : mf.cpp mf.h $(HEADERS)
	$(CC) $(CCFLAGS) -c mf.cpp

factorgraph.o : factorgraph.cpp factorgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c factorgraph.cpp

util.o : util.cpp util.h $(HEADERS)
	$(CC) $(CCFLAGS) -c util.cpp

regiongraph.o : regiongraph.cpp regiongraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c regiongraph.cpp

hak.o : hak.cpp hak.h $(HEADERS) regiongraph.h
	$(CC) $(CCFLAGS) -c hak.cpp

clustergraph.o : clustergraph.cpp clustergraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c clustergraph.cpp

jtree.o : jtree.cpp jtree.h $(HEADERS) weightedgraph.h clustergraph.h regiongraph.h
	$(CC) $(CCFLAGS) -c jtree.cpp

treeep.o : treeep.cpp treeep.h $(HEADERS) weightedgraph.h clustergraph.h regiongraph.h jtree.h
	$(CC) $(CCFLAGS) -c treeep.cpp

weightedgraph.o : weightedgraph.cpp weightedgraph.h $(HEADERS)
	$(CC) $(CCFLAGS) -c weightedgraph.cpp

mr.o : mr.cpp mr.h $(HEADERS)
	$(CC) $(CCFLAGS) -c mr.cpp

properties.o : properties.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c properties.cpp

alldai.o : alldai.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c alldai.cpp

x2x.o : x2x.cpp $(HEADERS)
	$(CC) $(CCFLAGS) -c x2x.cpp

# EXAMPLE
##########

example : example.cpp $(HEADERS) libdai.a
	$(CC) $(CCFLAGS) -o example example.cpp -L. libdai.a

# TESTS
########

tests/test : tests/test.cpp $(HEADERS) libdai.a
	$(CC) $(CCFLAGS) -o tests/test tests/test.cpp -L. -ldai $(BOOSTFLAGS)


# MATLAB INTERFACE
###################

matlab/dai.$(MEXEXT) : matlab/dai.cpp $(HEADERS) matlab/matlab.o daialg.o alldai.o bp.o clustergraph.o factorgraph.o hak.o jtree.o lc.o mf.o mr.o properties.o regiongraph.o util.o treeep.o weightedgraph.o x2x.o
	$(MEX) $(MEXFLAGS) -o matlab/dai matlab/dai.cpp matlab/matlab.o daialg.o alldai.o bp.o clustergraph.o factorgraph.o hak.o jtree.o lc.o mf.o mr.o properties.o regiongraph.o util.o treeep.o weightedgraph.o x2x.o

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
