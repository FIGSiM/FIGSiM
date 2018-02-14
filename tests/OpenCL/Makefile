####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Makefile

# compile tests

PREFIX=../../

include $(PREFIX)Makefile.choices
include $(PREFIX)Makefile.check

# Program name
PROG_NAME=OpenCL

# Define object code files
ifeq ($(TEST_OPENCL), no)
MINICL_OBJS = btAlignedAllocator.o btThreadSupportInterface.o Win32ThreadSupport.o PosixThreadSupport.o SequentialThreadSupport.o MiniCLTaskScheduler.o MiniCLTask.o MiniCL.o
else
MINICL_OBJS = $()
endif
OBJS = ScalarMat.o VecMat.o 

# Compile executable and link
all: $(OBJS) $(MINICL_OBJS)
	$(CPP) $(CPPFLAGS) -o $(PROG_NAME) $(PROG_NAME).cpp $(LIB_OPENCL) $(OBJS) $(MINICL_OBJS)

# Compile object code files without linking
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

clean:
	rm -f *.o $(PROG_NAME)

test: all
	./$(PROG_NAME)
