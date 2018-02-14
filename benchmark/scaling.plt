#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

for datafile in `ls *scaling_results*.dat`; do
	datafile=${datafile:0:${#datafile}-4}
	export FILENAME="${datafile}"
	echo "Creating plot for ${datafile} ...";
	gnuplot << EOF
filename="`echo $FILENAME`"

steps = system("awk -F'\\t' '(\$1==\"#steps:\"){print \$2};(NR>16){exit}' ".filename.".dat")
if(steps eq "") steps=2000
randsteps = system("awk -F'\\t' '(\$1==\"#randsteps:\"){print \$2};(NR>16){exit}' ".filename.".dat")
if(randsteps eq "") randsteps=1000
fullsteps=steps-randsteps
Title = system("awk -F'\\t' '(\$1!=\"#\"){first=\$1};(NR>16){exit}END{print first}' ".filename.".dat")
Title = Title." scaling results"

set term png enhanced font Arial 16 size 1280,1024
set samples 1000

set style line 1 lt 1 lw 4 lc rgb "#820808"
set style line 2 lt 2 lw 4 lc rgb "#088208"
set style line 3 lt 2 lw 4 lc rgb "#080882"
set style line 4 lt 2 lw 4 lc rgb "#820882"

set title Title
set xlabel "N"
set ylabel "normalized runtime [µs/cycle]"

set output filename.".png"

set fit errorvariables
a2=0.5
b2=2.0
b2_err=0.0
fit a2*x**b2 filename.'.dat' u 2:(\$3/randsteps*1E6) via a2
af2=1.0
bf2=2.0
bf2_err=0.0
fit af2*x**bf2 filename.'.dat' u 2:((\$4-\$3)/fullsteps*1E6) via af2
plot filename.'.dat' u 2:(\$3/randsteps*1E6) t "randomization phase" w p pt 7 ps 2, a2*x**b2 t sprintf("a*N^n exponent fit (a=%.3f+/-%.3f µs/cycle, n=%.3f+/-%.3f)",a2,a2_err,b2,b2_err) w l ls 1, filename.'.dat' u 2:((\$4-\$3)/fullsteps*1E6) t "full energy calculations" w p pt 7 ps 2, af2*x**bf2 t sprintf("a*N^n exponent fit (a=%.3f+/-%.3f µs/cycle, n=%.3f+/-%.3f)",af2,af2_err,bf2,bf2_err) w l ls 2
EOF
	echo "-> Done."
done
