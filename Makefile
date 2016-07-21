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
ALGLIBDIR = alglib
SRCDIR  = src
OBJDIR  = build
TGTDIR  = release
SWIGDIR = swig
PKGDIR  = python_package
PKGIDDIR = idcpp
BINDEST_DIR = /usr/local/bin
LIBDEST_DIR = /usr/local/lib
INCDEST_DIR = /usr/local/include

ALGLIBSRC_CPP  = 	ap.cpp \
									alglibinternal.cpp \
									alglibmisc.cpp \
									linalg.cpp \
									integration.cpp \
									interpolation.cpp \
									optimization.cpp \
									solvers.cpp \
									specialfunctions.cpp \

LIBSOURCES_CPP = 	magnet.cpp \
									block.cpp \
									blockcontainer.cpp \
									cassette.cpp \
									cassettecontainer.cpp \
									fieldmap.cpp \
									fieldmap3d.cpp \
									fieldmapcontainer.cpp \
									functions.cpp \
									grid.cpp \
									mask.cpp \
									kickmap.cpp \
									idmodel.cpp \

BINSOURCES_CPP =	generatekickmap.cpp

AUXFILES  = VERSION

LIBS = -lm
INC  = -I./$(INCDIR) -I./$(INCDIR)/$(ALGLIBDIR) -I/usr/include/python3.4

$(shell touch $(SRCDIR)/generatekickmap.cpp) # this is so that last compilation time always goes into executable

ifeq ($(MAKECMDGOALS),idcpp-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif

ALGLIBOBJS  = $(addprefix $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/, $(ALGLIBSRC_CPP:.cpp=.o))
LIBOBJECTS  = $(addprefix $(OBJDIR)/$(TGTDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/$(TGTDIR)/, $(BINSOURCES_CPP:.cpp=.o))
LDFLAGS    = $(MACHINE)

.PHONY: all alllibs clean cleanall

#### TARGETS ####

all: alglib libidcpp lnls-generate-kickmap python_package

python_package: $(PKGDIR)/$(PKGIDDIR)/idcpp.py $(PKGDIR)/$(PKGIDDIR)/_idcpp.so

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)/, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)/$(ALGLIBDIR)/, $(ALGLIBSRC_CPP)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP)) | sed 's/.*\.o/$(OBJDIR)\/$(TGTDIR)\/&/' > .depend)
-include .depend

lnls-generate-kickmap: $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap

alglib: $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a

libidcpp: $(OBJDIR)/$(TGTDIR)/libidcpp.a

$(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a: $(ALGLIBOBJS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/$(TGTDIR)/libidcpp.a: $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap: $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a libidcpp $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a $(OBJDIR)/$(TGTDIR)/libidcpp.a $(LIBS) -o $@

$(PKGDIR)/$(PKGIDDIR)/idcpp.py: $(PKGDIR)/$(SWIGDIR)/idcpp.py | $(PKGDIR)/$(PKGIDDIR)
	cp $(PKGDIR)/$(SWIGDIR)/idcpp.py $(PKGDIR)/$(PKGIDDIR)

$(PKGDIR)/$(PKGIDDIR)/_idcpp.so: $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/alglib.a libidcpp $(OBJDIR)/$(TGTDIR)/libidcpp.so | $(PKGDIR)/$(PKGIDDIR)
	cp $(OBJDIR)/$(TGTDIR)/libidcpp.so $(PKGDIR)/$(PKGIDDIR)/_idcpp.so

$(OBJDIR)/$(TGTDIR)/libidcpp.so: $(OBJDIR)/$(TGTDIR)/idcpp_wrap.o $(LIBOBJECTS) | $(OBJDIR)/$(TGTDIR)
	$(CXX) -shared $(LDFLAGS) $(LIBOBJECTS) $(ALGLIBOBJS) $(OBJDIR)/$(TGTDIR)/idcpp_wrap.o $(LIBS) -o $@

$(OBJDIR)/$(TGTDIR)/idcpp_wrap.o: $(PKGDIR)/$(SWIGDIR)/idcpp_wrap.cxx | $(OBJDIR)/$(TGTDIR)
	$(CXX) -c $(CFLAGS) $(INC) $< -o $@

$(PKGDIR)/$(SWIGDIR)/idcpp.py $(PKGDIR)/$(SWIGDIR)/idcpp_wrap.cxx: $(PKGDIR)/$(SWIGDIR)/idcpp.i $(LIBOBJECTS)
	swig -c++ -python $(INC) $(PKGDIR)/$(SWIGDIR)/idcpp.i

$(ALGLIBOBJS): | $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)

$(LIBOBJECTS): | $(OBJDIR)/$(TGTDIR)

$(BINOBJECTS): | $(OBJDIR)/$(TGTDIR)

$(PKGDIR)/$(PKGIDDIR):
	mkdir -p $(PKGDIR)/$(PKGIDDIR)

$(OBJDIR)/$(TGTDIR):
	mkdir -p $(OBJDIR)/$(TGTDIR)

$(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR):
	mkdir -p $(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)

$(BINDEST_DIR):
	mkdir $(BINDEST_DIR)

$(LIBDEST_DIR):
	mkdir $(LIBDEST_DIR)

$(INCDEST_DIR):
	mkdir $(INCDEST_DIR)

install: uninstall all
	cp $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap $(BINDEST_DIR)
	cp $(OBJDIR)/$(TGTDIR)/libidcpp.a $(LIBDEST_DIR)
	$(MAKE) install -C $(PKGDIR)

develop: uninstall all
	ln -srf $(OBJDIR)/$(TGTDIR)/lnls-generate-kickmap $(BINDEST_DIR)
	ln -srf $(OBJDIR)/$(TGTDIR)/libidcpp.a $(LIBDEST_DIR)
	$(MAKE) develop -C $(PKGDIR)

clean:
	-rm -rf $(OBJDIR) .depend *.out *.dat *~ *.o *.a *.txt $(PKGDIR)/$(PKGIDDIR)/idcpp.py $(PKGDIR)/$(PKGIDDIR)/idcpp_wrap.cxx
	$(MAKE) clean -C $(PKGDIR)

uninstall:
	-rm -rf $(BINDEST_DIR)/lnls-generate-kickmap
	-rm -rf $(LIBDEST_DIR)/libidcpp.a

cleanall: clean

#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/$(TGTDIR)/$(ALGLIBDIR)/%.o: $(SRCDIR)/$(ALGLIBDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR)/$(ALGLIBDIR) $< -o $@;

$(OBJDIR)/$(TGTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
