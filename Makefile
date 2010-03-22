# This file is part of libDAI - http://www.libdai.org/
#
# libDAI is licensed under the terms of the GNU General Public License version
# 2, or (at your option) any later version. libDAI is distributed without any
# warranty. See the file COPYING for more details.
#
# Copyright (C) 2006-2010  Joris Mooij  [joris dot mooij at libdai dot org]
# Copyright (C) 2006-2007  Radboud University Nijmegen, The Netherlands


# Load the platform independent build configuration file
include Makefile.ALL

# Load the local configuration from Makefile.conf
include Makefile.conf

# Set version and date
DAI_VERSION="git HEAD"
DAI_DATE="February 11, 2010 - or later"

# Directories of libDAI sources
# Location libDAI headers
INC=include/dai
# Location of libDAI source files
SRC=src
# Destination directory of libDAI library
LIB=lib

# Set final compiler flags
ifdef DEBUG
  CCFLAGS:=$(CCFLAGS) $(CCDEBUGFLAGS)
else
  CCFLAGS:=$(CCFLAGS) $(CCNODEBUGFLAGS)
endif

# Define build targets
TARGETS=tests utils lib examples unittests testregression testem
ifdef WITH_DOC
  TARGETS:=$(TARGETS) doc
endif
ifdef WITH_MATLAB
  TARGETS:=$(TARGETS) matlabs
  # Specify the same C++ compiler and flags to mex
  ifneq ($(OS),WINDOWS)
    MEXFLAGS=CXX\#$(CC) CXXFLAGS\#'$(CCFLAGS)'
  else
    MEXFLAGS=CXX\#$(CC) CXXFLAGS\#"$(CCFLAGS)"
  endif
  ifdef NEW_MATLAB
    MEXFLAGS:=$(MEXFLAGS) -largeArrayDims
  else
    MEXFLAGS:=$(MEXFLAGS) -DSMALLMEM
  endif
endif

# Define conditional build targets
OBJECTS:=exactinf$(OE) evidence$(OE) emalg$(OE)
ifdef WITH_BP
  CCFLAGS:=$(CCFLAGS) -DDAI_WITH_BP
  OBJECTS:=$(OBJECTS) bp$(OE)
endif
ifdef WITH_FBP
  CCFLAGS:=$(CCFLAGS) -DDAI_WITH_FBP
  OBJECTS:=$(OBJECTS) fbp$(OE)
endif
ifdef WITH_TRWBP
  CCFLAGS:=$(CCFLAGS) -DDAI_WITH_TRWBP
  OBJECTS:=$(OBJECTS) trwbp$(OE)
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
ifdef WITH_CBP
  CCFLAGS:=$(CCFLAGS) -DDAI_WITH_CBP
  OBJECTS:=$(OBJECTS) bbp$(OE) cbp$(OE) bp_dual$(OE)
endif

# Define standard libDAI header dependencies
HEADERS=$(INC)/bipgraph.h $(INC)/graph.h $(INC)/index.h $(INC)/var.h $(INC)/factor.h $(INC)/varset.h $(INC)/smallset.h $(INC)/prob.h $(INC)/daialg.h $(INC)/properties.h $(INC)/alldai.h $(INC)/enum.h $(INC)/exceptions.h $(INC)/util.h

# Setup final command for C++ compiler and MEX
ifneq ($(OS),WINDOWS)
  CC:=$(CC) $(CCINC) $(CCFLAGS) $(CCLIB)
else
  CC:=$(CC) $(CCINC) $(CCFLAGS)
  LIBS:=$(LIBS) $(CCLIB)
endif
MEX:=$(MEX) $(CCLIB) $(CCINC) $(MEXFLAGS)


# META TARGETS
###############

all : $(TARGETS)

examples : examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE) examples/example_permute$(EE) examples/example_sprinkler$(EE) examples/example_sprinkler_gibbs$(EE) examples/example_sprinkler_em$(EE)

matlabs : matlab/dai$(ME) matlab/dai_readfg$(ME) matlab/dai_writefg$(ME) matlab/dai_potstrength$(ME)

unittests : tests/unit/var$(EE) tests/unit/smallset$(EE) tests/unit/varset$(EE) tests/unit/graph$(EE) tests/unit/bipgraph$(EE)
	echo Running unit tests...
	tests/unit/var
	tests/unit/smallset
	tests/unit/varset
	tests/unit/graph
	tests/unit/bipgraph

tests : tests/testdai$(EE) tests/testem/testem$(EE) tests/testbbp$(EE) $(unittests)

utils : utils/createfg$(EE) utils/fg2dot$(EE) utils/fginfo$(EE)

lib: $(LIB)/libdai$(LE)


# OBJECTS
##########

%$(OE) : $(SRC)/%.cpp $(INC)/%.h $(HEADERS)
	$(CC) -c $< -o $@

bbp$(OE) : $(SRC)/bbp.cpp $(INC)/bbp.h $(INC)/bp_dual.h $(HEADERS)
	$(CC) -c $<

cbp$(OE) : $(SRC)/cbp.cpp $(INC)/cbp.h $(INC)/bbp.h $(INC)/bp_dual.h $(HEADERS)
	$(CC) -c $<

hak$(OE) : $(SRC)/hak.cpp $(INC)/hak.h $(HEADERS) $(INC)/regiongraph.h
	$(CC) -c $<

jtree$(OE) : $(SRC)/jtree.cpp $(INC)/jtree.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h
	$(CC) -c $<

treeep$(OE) : $(SRC)/treeep.cpp $(INC)/treeep.h $(HEADERS) $(INC)/weightedgraph.h $(INC)/clustergraph.h $(INC)/regiongraph.h $(INC)/jtree.h
	$(CC) -c $<

emalg$(OE) : $(SRC)/emalg.cpp $(INC)/emalg.h $(INC)/evidence.h $(HEADERS)
	$(CC) -c $<


# EXAMPLES
###########

examples/example$(EE) : examples/example.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_bipgraph$(EE) : examples/example_bipgraph.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_varset$(EE) : examples/example_varset.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_permute$(EE) : examples/example_permute.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_sprinkler$(EE) : examples/example_sprinkler.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_sprinkler_gibbs$(EE) : examples/example_sprinkler_gibbs.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

examples/example_sprinkler_em$(EE) : examples/example_sprinkler_em.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)


# UNIT TESTS
#############

tests/unit/var$(EE) : tests/unit/var.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)
tests/unit/smallset$(EE) : tests/unit/smallset.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)
tests/unit/varset$(EE) : tests/unit/varset.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)
tests/unit/graph$(EE) : tests/unit/graph.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)
tests/unit/bipgraph$(EE) : tests/unit/bipgraph.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_UTF)


# TESTS
########

tests/testdai$(EE) : tests/testdai.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)
tests/testem/testem$(EE) : tests/testem/testem.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)

tests/testbbp$(EE) : tests/testbbp.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)


# MATLAB INTERFACE
###################

matlab/dai$(ME) : $(SRC)/matlab/dai.cpp $(HEADERS) matlab$(OE) $(LIB)/libdai$(LE)
	$(MEX) -o$@ $< matlab$(OE) $(LIB)/libdai$(LE)

matlab/dai_readfg$(ME) : $(SRC)/matlab/dai_readfg.cpp $(HEADERS) factorgraph$(OE) matlab$(OE) exceptions$(OE) bipgraph$(OE)
	$(MEX) -o$@ $< factorgraph$(OE) matlab$(OE) exceptions$(OE) bipgraph$(OE)

matlab/dai_writefg$(ME) : $(SRC)/matlab/dai_writefg.cpp $(HEADERS) factorgraph$(OE) matlab$(OE) exceptions$(OE) bipgraph$(OE)
	$(MEX) -o$@ $< factorgraph$(OE) matlab$(OE) exceptions$(OE) bipgraph$(OE)

matlab/dai_potstrength$(ME) : $(SRC)/matlab/dai_potstrength.cpp $(HEADERS) matlab$(OE) exceptions$(OE)
	$(MEX) -o$@ $< matlab$(OE) exceptions$(OE)

matlab$(OE) : $(SRC)/matlab/matlab.cpp $(INC)/matlab/matlab.h $(HEADERS)
	$(MEX) -c $<


# UTILS
########

utils/createfg$(EE) : utils/createfg.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS) $(BOOSTLIBS_PO)

utils/fg2dot$(EE) : utils/fg2dot.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)

utils/fginfo$(EE) : utils/fginfo.cpp $(HEADERS) $(LIB)/libdai$(LE)
	$(CC) $(CCO)$@ $< $(LIBS)


# LIBRARY
##########

OBJECTS:=bipgraph$(OE) graph$(OE) varset$(OE) daialg$(OE) alldai$(OE) clustergraph$(OE) factor$(OE) factorgraph$(OE) properties$(OE) regiongraph$(OE) util$(OE) weightedgraph$(OE) exceptions$(OE) $(OBJECTS) 

ifneq ($(OS),WINDOWS)
$(LIB)/libdai$(LE) : $(OBJECTS)
	-mkdir -p lib
	ar rcus $(LIB)/libdai$(LE) $(OBJECTS)
else
$(LIB)/libdai$(LE) : $(OBJECTS)
	-mkdir lib
	lib /out:$(LIB)/libdai$(LE) $(OBJECTS)
endif


# REGRESSION TESTS
###################

testregression : tests/testdai$(EE)
	@echo Starting regression test...this can take a minute or so!
ifneq ($(OS),WINDOWS)
	cd tests && ./testregression && cd ..
else
	cd tests && testregression.bat && cd ..
endif

testem : tests/testem/testem$(EE)
	@echo Starting EM tests
ifneq ($(OS),WINDOWS)
	cd tests/testem && ./runtests && cd ../..
else
	cd tests\testem && runtests && cd ..\..
endif


# DOCUMENTATION
################

doc : $(INC)/*.h $(SRC)/*.cpp examples/*.cpp doxygen.conf
	doxygen doxygen.conf

README : doc scripts/makeREADME
	DAI_VERSION=$(DAI_VERSION) DAI_DATE=$(DAI_DATE) scripts/makeREADME

TAGS :
	etags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp
	ctags src/*.cpp include/dai/*.h tests/*.cpp utils/*.cpp


# CLEAN
########

ifneq ($(OS),WINDOWS)
.PHONY : clean
clean :
	-rm *$(OE)
	-rm matlab/*$(ME)
	-rm examples/example$(EE) examples/example_bipgraph$(EE) examples/example_varset$(EE) examples/example_permute$(EE) examples/example_sprinkler$(EE) examples/example_sprinkler_gibbs$(EE) examples/example_sprinkler_em$(EE)
	-rm tests/testdai$(EE) tests/testem/testem$(EE) tests/testbbp$(EE)
	-rm tests/unit/var$(EE) tests/unit/smallset$(EE) tests/unit/varset$(EE) tests/unit/graph$(EE) tests/unit/bipgraph$(EE)
	-rm utils/fg2dot$(EE) utils/createfg$(EE) utils/fginfo$(EE)
	-rm -R doc
	-rm -R lib
else
.PHONY : clean
clean :
	-del *$(OE)
	-del *.ilk
	-del *.pdb
	-del *$(EE)
	-del matlab\*$(ME)
	-del examples\*$(EE)
	-del examples\*.ilk
	-del examples\*.pdb
	-del tests\testdai$(EE)
	-del tests\testbbp$(EE)
	-del tests\testem\testem$(EE)
	-del tests\*.pdb
	-del tests\*.ilk
	-del tests\testem\*.pdb
	-del tests\testem\*.ilk
	-del utils\*$(EE)
	-del utils\*.pdb
	-del utils\*.ilk
	-del tests\unit\*.pdk
	-del tests\unit\*.ilk
	-del tests\unit\var$(EE)
	-del tests\unit\smallset$(EE)
	-del tests\unit\varset$(EE)
	-del tests\unit\graph$(EE)
	-del tests\unit\bipgraph$(EE)
	-del $(LIB)\libdai$(LE)
	-rmdir lib
endif
