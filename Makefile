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

PREFIX=./

include Makefile.choices
include Makefile.check

#Program name
PROG_NAME=MCfig
BENCH_TOOL=MCbench
X3D_TOOL=traj2x3d
MOL2_TOOL=traj2mol2
CORR_TOOL=traj2dat
STAT_TOOL=traj2stat
LOD_TOOL=fit2lod
LODPLOT_TOOL=lod2plot
PARAMS_TOOL=MCparams

# Define object code files
ifeq ($(TEST_OPENCL), no)
	MINICL_OBJS = btAlignedAllocator.o btThreadSupportInterface.o Win32ThreadSupport.o PosixThreadSupport.o SequentialThreadSupport.o MiniCLTaskScheduler.o MiniCLTask.o MiniCL.o
else
	MINICL_OBJS = $()
endif
OBJS = ScalarMat.o VecMat.o Latticeoids3d.o MC_Elements.o ConfigReader.o MC_Config.o setup.o
OBJS_X3DTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o setup.o
OBJS_MOL2TOOL = ScalarMat.o VecMat.o MC_Elements.o ConfigReader.o MC_Config.o setup.o
OBJS_CORRTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o setup.o
OBJS_STATTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o setup.o
OBJS_LODTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o setup.o
OBJS_LODPLOTTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o setup.o
OBJS_PARAMSTOOL = ScalarMat.o VecMat.o ConfigReader.o MC_Config.o

# Find out where we're compiling on
UNAME_CMD := $(shell uname -s)
OBJLINK = -Wl,--format=binary -Wl,Source.tar.gz -Wl,--format=default
ifeq ($(UNAME_CMD),Darwin)
	OBJLINK = -sectcreate __DATA __Source_tar_gz Source.tar.gz 
endif

ifeq ($(OPTIMIZE), y)
	PROFILE=-fprofile-generate
ifneq ("$(wildcard optimized.make)","")
	PROFILE=-fprofile-use
endif
else
	PROFILE=$()
endif

#Compile executable and link
all:
ifeq ($(PROFILE), -fprofile-generate)
	make optimize
	make $(PROG_NAME)
	make $(BENCH_TOOL)
	make -C $(PREFIX)tests $@ $(MAKECMDGOALS)
else
	make $(OBJS) $(MINICL_OBJS)
	make -C $(PREFIX)tests $@ $(MAKECMDGOALS)
	make repository
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(PROG_NAME) MCmain.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(BENCH_TOOL) MCbench.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	rm -f Source.tar.gz
endif
	$(CPP) $(CPPFLAGS) -o $(X3D_TOOL) traj2x3d.cpp $(OBJS_X3DTOOL) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(MOL2_TOOL) traj2mol2.cpp $(OBJS_MOL2TOOL) $(MINICL_OBJS) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(CORR_TOOL) traj2dat.cpp $(OBJS_CORRTOOL) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(STAT_TOOL) traj2stat.cpp $(OBJS_STATTOOL) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(LOD_TOOL) fit2lod.cpp $(OBJS_LODTOOL) $(MINICL_OBJS) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(LODPLOT_TOOL) lod2plot.cpp $(OBJS_LODPLOTTOOL) $(MINICL_OBJS) $(LIB_OPENCL)
	$(CPP) $(CPPFLAGS) -o $(PARAMS_TOOL) MCparams.cpp $(OBJS_PARAMSTOOL) $(LIB_OPENCL)

# Compile object code files without linking
%.o:$(PREFIX)%.cpp
	$(CPP) $(CPPFLAGS) $(PROFILE) $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/MiniCL/%.cpp
	$(CPP) $(CPPFLAGS) $(PROFILE)  $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/MiniCL/MiniCLTask/%.cpp
	$(CPP) $(CPPFLAGS) $(PROFILE)  $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/BulletMultiThreaded/%.cpp
	$(CPP) $(CPPFLAGS) $(PROFILE)  $(ALL_INCLUDES) -c $< -o $@

%.o:$(PREFIX)/bulletphysics/LinearMath/%.cpp
	$(CPP) $(CPPFLAGS) $(PROFILE)  $(ALL_INCLUDES) -c $< -o $@

# folder docs exists, hence is a "phony" target ...
.PHONY: $(PROG_NAME) $(BENCH_TOOL) $(X3D_TOOL) $(MOL2_TOOL) $(CORR_TOOL) $(STAT_TOOL) $(LOD_TOOL) $(LODPLOT_TOOL) $(PARAMS_TOOL) docs tests clean repository

clean:
	rm -f *.o *.gcda $(PROG_NAME) $(BENCH_TOOL) $(X3D_TOOL) $(MOL2_TOOL) $(CORR_TOOL) $(STAT_TOOL) $(LOD_TOOL) $(LODPLOT_TOOL) $(PARAMS_TOOL) optimized.make
	make -C $(PREFIX)docs $@ $(MAKECMDGOALS)
	make -C $(PREFIX)tests $@ $(MAKECMDGOALS)
	rm -f $(PREFIX)tests/test_opencl
	rm -f Source.tar.gz
	echo "" > $(PREFIX)tests/use_opencl.h

$(PROG_NAME):
ifeq ($(PROFILE), -fprofile-generate)
	make optimize
	make $(PROG_NAME)
else
	make $(OBJS) $(MINICL_OBJS)
	make repository
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(PROG_NAME) MCmain.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	rm -f Source.tar.gz
endif

$(BENCH_TOOL):
ifeq ($(PROFILE), -fprofile-generate)
	make optimize
	make $(BENCH_TOOL)
else
	make $(OBJS) $(MINICL_OBJS)
	make repository
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(BENCH_TOOL) MCbench.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	rm -f Source.tar.gz
endif

$(X3D_TOOL): $(OBJS_X3DTOOL)
	$(CPP) $(CPPFLAGS) -o $(X3D_TOOL) traj2x3d.cpp $(OBJS_X3DTOOL) $(LIB_OPENCL)

$(MOL2_TOOL): $(OBJS_MOL2TOOL) $(MINICL_OBJS)
	$(CPP) $(CPPFLAGS) -o $(MOL2_TOOL) traj2mol2.cpp $(OBJS_MOL2TOOL) $(MINICL_OBJS) $(LIB_OPENCL)

$(CORR_TOOL): $(OBJS_CORRTOOL)
	$(CPP) $(CPPFLAGS) -o $(CORR_TOOL) traj2dat.cpp $(OBJS_CORRTOOL) $(LIB_OPENCL)

$(STAT_TOOL): $(OBJS_STATTOOL)
	$(CPP) $(CPPFLAGS) -o $(STAT_TOOL) traj2stat.cpp $(OBJS_STATTOOL) $(LIB_OPENCL)

$(LOD_TOOL): $(OBJS_LODTOOL) $(MINICL_OBJS)
	$(CPP) $(CPPFLAGS) -o $(LOD_TOOL) fit2lod.cpp $(OBJS_LODTOOL) $(MINICL_OBJS) $(LIB_OPENCL)

$(LODPLOT_TOOL): $(OBJS_LODPLOTTOOL) $(MINICL_OBJS)
	$(CPP) $(CPPFLAGS) -o $(LODPLOT_TOOL) lod2plot.cpp $(OBJS_LODPLOTTOOL) $(MINICL_OBJS) $(LIB_OPENCL)

$(PARAMS_TOOL): $(OBJS_PARAMSTOOL)
	$(CPP) $(CPPFLAGS) -o $(PARAMS_TOOL) MCparams.cpp $(OBJS_PARAMSTOOL) $(LIB_OPENCL)

tests:
	make -C $(PREFIX)tests

docs:
	make -C $(PREFIX)docs $@ $(MAKECMDGOALS)

optimize:
	@echo "***********************************************************"
	@echo "*                                                         *"
	@echo "* Compiling an optimized build, this may take a while ... *"
	@echo "*                                                         *"
	@echo "***********************************************************"
	rm -f *.gcda *.o $(PROG_NAME) Source.tar.gz optimized.make
	touch Source.tar.gz
	make $(OBJS) $(MINICL_OBJS)
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(PROG_NAME) MCmain.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	$(CPP) $(CPPFLAGS) $(PROFILE) -o $(BENCH_TOOL) MCbench.cpp $(OBJS) $(MINICL_OBJS) $(LIB_OPENCL) $(OBJLINK)
	@echo "***********************************************************"
	@echo "*                                                         *"
	echo "*                  Running benchmarks ...                 *"
	@echo "*                                                         *"
	@echo "***********************************************************"
	./$(BENCH_TOOL) benchmark/configs/test0
	./$(PROG_NAME) benchmark/configs/test0
	./$(PROG_NAME) benchmark/configs/touch_scaling
	./$(PROG_NAME) benchmark/configs/test4
	rm -f *.o $(PROG_NAME) $(BENCH_TOOL) *_0.* Source.tar.gz
	touch optimized.make
	@echo "***********************************************************"
	@echo "*                                                         *"
	@echo "*                  Recompiling things ...                 *"
	@echo "*                                                         *"
	@echo "***********************************************************"

repository:
	rm -rf $(PREFIX)Source
	mkdir $(PREFIX)Source
	cp $(PREFIX)LICENSE $(PREFIX)Source
	cp $(PREFIX)README.md $(PREFIX)Source
	cp $(PREFIX)CHANGELOG $(PREFIX)Source
	cp $(PREFIX)MCstat $(PREFIX)Source
	cp $(PREFIX)MCsubAT $(PREFIX)Source
	cp $(PREFIX)MCqueueInfo $(PREFIX)Source
	cp $(PREFIX)MCrun.sh $(PREFIX)Source
	cp $(PREFIX)Makefile* $(PREFIX)Source
	cp $(PREFIX)*.plt $(PREFIX)Source
	cp $(PREFIX)*.png $(PREFIX)Source
	cp $(PREFIX)*.cpp $(PREFIX)Source
	cp $(PREFIX)*.cl $(PREFIX)Source
	cp $(PREFIX)*.h $(PREFIX)Source
	mkdir $(PREFIX)Source/docs
	cp $(PREFIX)docs/Makefile $(PREFIX)Source/docs
	cp $(PREFIX)docs/Doxyfile $(PREFIX)Source/docs
	cp $(PREFIX)docs/*.pdf $(PREFIX)Source/docs
	cp $(PREFIX)docs/*.eps $(PREFIX)Source/docs
	mkdir $(PREFIX)Source/tests
	sed "s/TESTS =/\#TESTS =/g" $(PREFIX)tests/Makefile | sed "s/Makefile/Makefile\\nTESTS = /g" > $(PREFIX)Source/tests/Makefile
	cp $(PREFIX)tests/*.sh $(PREFIX)Source/tests
	cp $(PREFIX)tests/*.cpp $(PREFIX)Source/tests
	cp $(PREFIX)tests/*.h $(PREFIX)Source/tests
	cp -r $(PREFIX)bulletphysics $(PREFIX)Source/bulletphysics
	tar -cf Source.tar Source
	rm -rf $(PREFIX)Source
	gzip -f Source.tar


