## Insertion devices
## ========
## Author:      Accelerator Physics Group - LNLS
## contact:     xresende@gmail.com
## affiliation: Laboratorio Nacional de Luz Sincrotron
##
## The MIT License (MIT)
##
## Copyright (c) <year> <copyright holders>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.

#### READS LIB VERSION ####

FILE=VERSION
VERSION=$(shell cat ${FILE})

#### COMPILATION OPTIONS ####
CC          = gcc
CXX         = g++
AR          = ar
MACHINE     = -m64
OPT_FLAG    = -O3 -std=c++11 -fPIC
DBG_FLAG    = -O0 -g3 -std=c++11 -fPIC
ARFLAGS     = rcs
DFLAGS      = -DVERSION=$(VERSION)
PYTHON      = python3

INCDIR  = include
SRCDIR  = src
OBJDIR  = build
TGTDIR  = release
SWIGDIR = swig
PKGDIR  = python_package
PKGIDDIR = insertion_devices
BINDEST_DIR = /usr/local/bin
LIBDEST_DIR = /usr/local/lib
INCDEST_DIR = /usr/local/include

LIBSOURCES_CPP =  ap.cpp \
									alglibinternal.cpp \
									alglibmisc.cpp \
									linalg.cpp \
									integration.cpp \
									interpolation.cpp \
									optimization.cpp \
									solvers.cpp \
									specialfunctions.cpp \
									fieldmap.cpp \
									fieldmapcontainer.cpp \
									rungekutta.cpp \
									grid.cpp \
									mask.cpp \
									kickmap.cpp \

BINSOURCES_CPP =	generate_kickmap.cpp

BINSOURCES_CPP2 =	tests.cpp

AUXFILES  = VERSION

LIBS = -lm
INC  = -I./$(INCDIR) -I/usr/include/python3.4

$(shell touch $(SRCDIR)/generate_kickmap.cpp) # this is so that last compilation time always goes into executable
$(shell touch $(SRCDIR)/tests.cpp)

ifeq ($(MAKECMDGOALS),insertion_devices-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif


LIBOBJECTS  = $(addprefix $(OBJDIR)/$(TGTDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/$(TGTDIR)/, $(BINSOURCES_CPP:.cpp=.o))
BINOBJECTS2 = $(addprefix $(OBJDIR)/$(TGTDIR)/, $(BINSOURCES_CPP2:.cpp=.o))
LDFLAGS    = $(MACHINE)


.PHONY: all alllibs clean cleanall

#### TARGETS ####

all: libids lnls-generate-kickmap python_package

python_package: $(PKGDIR)/$(PKGIDDIR)/insertion_devices.py $(PKGDIR)/$(PKGIDDIR)/_insertion_devices.so

test: libids run_test

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)/, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP2)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP)) | sed 's/.*\.o/$(OBJDIR)\/$(TGTDIR)\/&/' > .depend)
-include .depend

lnls-generate-kickmap: $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap

run_test: $(OBJDIR)/$(TGTDIR)/run_test

libids: $(OBJDIR)/$(TGTDIR)/libids.a

$(OBJDIR)/$(TGTDIR)/libids.a: $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap: libids $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)/$(TGTDIR)/libids.a $(LIBS) -o $@

$(OBJDIR)/$(TGTDIR)/run_test: libids $(BINOBJECTS2)
	$(CXX) $(LDFLAGS) $(BINOBJECTS2) $(OBJDIR)/$(TGTDIR)/libids.a $(LIBS) -o $@

$(PKGDIR)/$(PKGIDDIR)/insertion_devices.py: $(PKGDIR)/$(SWIGDIR)/insertion_devices.py | $(PKGDIR)/$(PKGIDDIR)
	cp $(PKGDIR)/$(SWIGDIR)/insertion_devices.py $(PKGDIR)/$(PKGIDDIR)

$(PKGDIR)/$(PKGIDDIR)/_insertion_devices.so: libids $(OBJDIR)/$(TGTDIR)/libinsertion_devices.so | $(PKGDIR)/$(PKGIDDIR)
	cp $(OBJDIR)/$(TGTDIR)/libinsertion_devices.so $(PKGDIR)/$(PKGIDDIR)/_insertion_devices.so

$(OBJDIR)/$(TGTDIR)/libinsertion_devices.so: $(OBJDIR)/$(TGTDIR)/insertion_devices_wrap.o $(LIBOBJECTS) | $(OBJDIR)/$(TGTDIR)
	$(CXX) -shared $(LDFLAGS) $(LIBOBJECTS)  $(OBJDIR)/$(TGTDIR)/insertion_devices_wrap.o $(LIBS) -o $@

$(OBJDIR)/$(TGTDIR)/insertion_devices_wrap.o: $(PKGDIR)/$(SWIGDIR)/insertion_devices_wrap.cxx | $(OBJDIR)/$(TGTDIR)
	$(CXX) -c $(CFLAGS) $(INC) $< -o $@

$(PKGDIR)/$(SWIGDIR)/insertion_devices.py $(PKGDIR)/$(SWIGDIR)/insertion_devices_wrap.cxx: $(PKGDIR)/$(SWIGDIR)/insertion_devices.i $(LIBOBJECTS)
	swig -c++ -python $(INC) $(PKGDIR)/$(SWIGDIR)/insertion_devices.i

$(LIBOBJECTS): | $(OBJDIR)/$(TGTDIR)

$(BINOBJECTS): | $(OBJDIR)/$(TGTDIR)

$(BINOBJECTS2): | $(OBJDIR)/$(TGTDIR)

$(PKGDIR)/$(PKGIDDIR):
	mkdir -p $(PKGDIR)/$(PKGIDDIR)

$(OBJDIR)/$(TGTDIR):
	mkdir -p $(OBJDIR)/$(TGTDIR)

$(BINDEST_DIR):
	mkdir $(BINDEST_DIR)

$(LIBDEST_DIR):
	mkdir $(LIBDEST_DIR)

$(INCDEST_DIR):
	mkdir $(INCDEST_DIR)

install: uninstall all
	cp $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap $(BINDEST_DIR)
	cp $(OBJDIR)/$(TGTDIR)/libids.a $(LIBDEST_DIR)
	$(MAKE) install -C $(PKGDIR)

develop: uninstall all
	ln -srf $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap $(BINDEST_DIR)
	ln -srf $(OBJDIR)/$(TGTDIR)/libids.a $(LIBDEST_DIR)
	$(MAKE) develop -C $(PKGDIR)

clean:
	-rm -rf $(OBJDIR) .depend *.out *.dat *~ *.o *.a *.txt $(PKGDIR)/$(PKGIDDIR)/insertion_devices.py $(PKGDIR)/$(PKGIDDIR)/insertion_devices_wrap.cxx
	$(MAKE) clean -C $(PKGDIR)

uninstall:
	-rm -rf $(BINDEST_DIR)/lnls-generate-kickmap
	-rm -rf $(LIBDEST_DIR)/libids.a

cleanall: clean

#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/$(TGTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
