#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Test if PThread can be used for compiling

current_dir=`pwd`
script_dir=`dirname $0`
cd "$script_dir"
rm test_pthread &> /dev/null
$1 -pthread -o test_pthread test_pthread.cpp &> /dev/null
test -e test_pthread && echo yes || echo no
cd "$current_dir"

