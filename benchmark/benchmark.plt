#!/usr/bin/gnuplot
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

set term postscript color enhanced
set output "benchmark_results.eps"
set style fill solid 1.0 border -1
set boxwidth 0.8

set ylabel "Runtime in seconds"

plot 'benchmark_results.dat' u 0:5:xticlabels("") t "Overall" w boxes lc rgb "#4AA66F",'' u 0:5:5 w labels font "Helvetica, 8" center offset 0,-0.5 notitle,'' u 0:7 t "Randomization" w boxes lc rgb "#2A64A6",'' u 0:7:7 w labels font "Helvetica, 8" tc rgb "#FFFFFF" center offset 0,-0.5 notitle,'' u 0:(0.0):3 w labels font "Helvetica, 12" center offset 0,-0.8 tc rgb "#000000" notitle,'' u 0:(0.0):4 w labels font "Helvetica, 12" center offset 0,-2 tc rgb "#000000" notitle,'' u 0:($5/2):2 w labels font "Helvetica, 12" rotate by 90 center tc rgb "#000000" notitle

