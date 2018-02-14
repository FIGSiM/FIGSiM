####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

# Makefile

all: docs

.PHONY: docs clean

clean:
	rm -rf html
	rm -f *.aux *.log *.nav *.out *Model.pdf *.snm *.toc

docs:
	doxygen
	pdflatex *.tex
	pdflatex *.tex

