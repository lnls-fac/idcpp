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

BINSOURCES_CPP =	generate_kickmap.cpp

BINSOURCES_CPP2 =	tests.cpp

AUXFILES  = VERSION

LIBS = -lm
INC  = -I./include

OBJDIR = build
SRCDIR = src
BINDEST_DIR = input_files

$(shell touch $(SRCDIR)/generate_kickmap.cpp) # this is so that last compilation time always goes into executable
$(shell touch $(SRCDIR)/tests.cpp) # this is so that last compilation time always goes into executable


ifeq ($(MAKECMDGOALS),insertion_devices-debug)
  CFLAGS    = $(MACHINE) $(DBG_FLAG) $(DFLAGS) -pthread
else
  CFLAGS    = $(MACHINE) $(OPT_FLAG) $(DFLAGS) -pthread
endif

LIBOBJECTS  = $(addprefix $(OBJDIR)/, $(LIBSOURCES_CPP:.cpp=.o))
BINOBJECTS  = $(addprefix $(OBJDIR)/, $(BINSOURCES_CPP:.cpp=.o))
BINOBJECTS2 = $(addprefix $(OBJDIR)/, $(BINSOURCES_CPP2:.cpp=.o))
LDFLAGS    = $(MACHINE)


.PHONY: all alllibs clean cleanall

#### TARGETS ####

all: libids generate_kickmap

test: libids run_test

#### GENERATES DEPENDENCY FILE ####
$(shell $(CXX) -MM $(CFLAGS) $(addprefix $(SRCDIR)/, $(LIBSOURCES_CPP)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP2)) $(addprefix $(SRCDIR)/, $(BINSOURCES_CPP)) | sed 's/.*\.o/$(OBJDIR)\/&/' > .depend)
-include .depend

generate_kickmap: $(OBJDIR)/generate_kickmap

run_test: $(OBJDIR)/run_test

libids: $(OBJDIR)/libids.a

$(OBJDIR)/libids.a: $(LIBOBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(OBJDIR)/generate_kickmap: libids $(BINOBJECTS)
	$(CXX) $(LDFLAGS) $(BINOBJECTS) $(OBJDIR)/libids.a $(LIBS) -o $@
	-rm -rf $(BINDEST_DIR)/generate_kickmap
	ln -srf $(OBJDIR)/generate_kickmap $(BINDEST_DIR)

$(OBJDIR)/run_test: libids $(BINOBJECTS2)
	$(CXX) $(LDFLAGS) $(BINOBJECTS2) $(OBJDIR)/libids.a $(LIBS) -o $@
	-rm -rf $(BINDEST_DIR)/run_test
	ln -srf $(OBJDIR)/run_test $(BINDEST_DIR)

$(LIBOBJECTS): | $(OBJDIR)

$(BINOBJECTS): | $(OBJDIR)

$(BINOBJECTS2): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

clean:
	-rm -rf $(OBJDIR) run_test generate_kickmap .depend *.out *.dat *~ *.o *.a *.txt
	-rm -rf $(BINDEST_DIR)/generate_kickmap
	-rm -rf $(BINDEST_DIR)/run_test

cleanall: clean


#### RULES ####

*.cpp: VERSION
	touch $(SRCDIR)/*.cpp
*.cc: VERSION
	touch $(SRCDIR)/*.cc
*.c: VERSION
	touch $(SRCDIR)/*.c

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CFLAGS) $(INC) -I./$(SRCDIR) $< -o $@;
