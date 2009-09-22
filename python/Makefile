include ../Makefile.conf

.PHONY:	all

all:			_dai.so

_dai.so:		../lib/libdai.a dai.i
			$(SWIG) -python -c++ dai.i
			g++ -Wall -c -O3 -g -fPIC dai_wrap.cxx -I$(INCLUDE_PYTHON) -I$(INCLUDE_BOOST) -I../include
			g++ -shared dai_wrap.o -o _dai.so ../lib/libdai.a

.PHONY: clean

clean:		
			rm -f dai_wrap.cxx dai_wrap.o _dai.so dai.py dai.pyc
