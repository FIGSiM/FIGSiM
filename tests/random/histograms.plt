#!/usr/bin/gnuplot
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

filename="gauss"

set term png enhanced font Arial 14 size 1600,1200
set samples 1000

A=2.75E6
sigma=1.0/3.4
x0=0.5
gaussian(x,x0,A,sigma)=A*exp(-0.5*((x-x0)/sigma)**2)

set style histogram clustered gap 1
set style fill solid border -1

set xlabel "x"
set ylabel "counts"

fit gaussian(x,x0,A,sigma) filename.".dat" u 1:2 via A,sigma

set output filename.".png"
plot filename.".dat" u 1:2 w boxes t 'normal distribution', gaussian(x,x0,A,sigma) t 'Fit' w l ls 1 lw 4 lc rgb "#202080"

filename="gauss_notail"
fit gaussian(x,x0,A,sigma) filename.".dat" u 1:2 via A,sigma

set output filename.".png"
plot filename.".dat" u 1:2 w boxes t 'normal distribution', gaussian(x,x0,A,sigma) t 'Fit' w l ls 1 lw 4 lc rgb "#202080"

filename="ran2"
a=1E6
fit a filename.".dat" u 1:2 via a

set output filename.".png"
plot filename.".dat" u 1:2 w boxes t 'normal distribution', a t 'Fit' w l ls 1 lw 4 lc rgb "#202080"

filename="CMWC4096"
a=1E6
fit a filename.".dat" u 1:2 via a

set output filename.".png"
plot [0:] [0:] filename.".dat" u 1:2 w boxes t 'normal distribution', a t 'Fit' w l ls 1 lw 4 lc rgb "#202080"

