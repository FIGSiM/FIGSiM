####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

#Makefile

#Used to compile the Robinson MC code with the GNU C++ compiler for testing/deployment on any system with a relatively modern version of gcc.

#LEJ, 01/12/10
#AT, Jan 17,2011
# - works with GNU make
# - modularized a bit

PREFIX=../

include Makefile.choices
include $(PREFIX)Makefile.check

#Program name

PROG_NAME=MCfig

#Define object code files
ifeq ($(TEST_OPENCL), no)
MINICL_OBJS = btAlignedAllocator.o btThreadSupportInterface.o Win32ThreadSupport.o PosixThreadSupport.o SequentialThreadSupport.o MiniCLTaskScheduler.o MiniCLTask.o MiniCL.o
else
MINICL_OBJS = $()
endif
OBJS = ScalarMat.o VecMat.o Latticeoids3d.o MC_Elements.o ConfigReader.o MC_Config.o setup.o

# Find out where we're compiling on
UNAME_CMD := $(shell uname -s)
OBJLINK = -Wl,--format=binary -Wl,Source.tar.gz -Wl,--format=default
ifeq ($(UNAME_CMD),Darwin)
	OBJLINK = -sectcreate __DATA __Source_tar_gz Source.tar.gz 
endif

#Compile executable and link
all: $(OBJS) $(MINICL_OBJS)
	touch Source.tar.gz
	$(CPP) $(CPPFLAGS) -o $(PROG_NAME) $(PREFIX)MCmain.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	rm Source.tar.gz

#Compile object code files without linking
%.o:$(PREFIX)%.cpp
	$(CPP) $(CPPFLAGS) $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/MiniCL/%.cpp
	$(CPP) $(CPPFLAGS) $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/MiniCL/MiniCLTask/%.cpp
	$(CPP) $(CPPFLAGS) $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/BulletMultiThreaded/%.cpp
	$(CPP) $(CPPFLAGS) $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/LinearMath/%.cpp
	$(CPP) $(CPPFLAGS) $(ALL_INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o $(PROG_NAME) *.traj *_summary_*.dat fit.log configs/*.scale
