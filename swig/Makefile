# This file is part of libDAI - http://www.libdai.org/
#
# libDAI is licensed under the terms of the GNU General Public License version
# 2, or (at your option) any later version. libDAI is distributed without any
# warranty. See the file COPYING for more details.
#
# Copyright (C) 2009  Patrick Pletscher  [pletscher at inf dot ethz dot ch]


include ../Makefile.conf

.PHONY:	all

all:			_dai.so dai.oct

_dai.so:		../lib/libdai.a dai.i
			$(SWIG) -python -classic -c++ dai.i
			g++ -Wall -c -O3 -g -fPIC dai_wrap.cxx -I$(INCLUDE_PYTHON) -I$(INCLUDE_BOOST) -I../include
			g++ -shared dai_wrap.o -o _dai.so ../lib/libdai.a

dai.oct:		../lib/libdai.a dai.i
			$(SWIG) -octave -c++ -o dai_wrap.cpp dai.i
			mkoctfile -I$(INCLUDE_BOOST) -I../include -o dai.oct dai_wrap.cpp ../lib/libdai.a

.PHONY: clean

clean:		
			rm -f dai_wrap.cxx dai_wrap.cpp dai_wrap.o _dai.so dai.oct dai.py dai.pyc sprinkler.fg
