####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Makefile
TESTS = OpenCL invsqrt random testcomplex testRF testmatrix

SUBDIRS = 
SUBDIRS += $(TESTS)

export TOPDIR = $(CURDIR)

all: $(SUBDIRS) 

.PHONY: subdirs $(SUBDIRS) clean test quicktest

subdirs: $(SUBDIRS)

$(SUBDIRS):
	make -C $@ $(MAKECMDGOALS)

clean: subdirs

test: subdirs

quicktest: subdirs
