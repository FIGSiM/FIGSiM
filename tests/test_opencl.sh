#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Test if OpenCL can be used for compiling

current_dir=`pwd`
script_dir=`dirname $0`
cd "$script_dir"
echo "#define OpenCL" > use_opencl.h
rm test_opencl &> /dev/null
$1 $2 -o test_opencl test_opencl.cpp &> /dev/null
test -e test_opencl && echo yes || echo no

if test -e test_opencl ; then
	echo "#define OpenCL" > use_opencl.h
else
	echo "" > use_opencl.h
fi
cd "$current_dir"

