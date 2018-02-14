#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

for datafile in `ls *.dat`; do
	datafile=${datafile:0:${#datafile}-4}
	export FILENAME="${datafile}"
	echo "Creating plot for ${datafile} ...";
	gnuplot << EOF
filename="`echo $FILENAME`"

set term png enhanced font Arial 14 size 1600,1200
set samples 1000

set style line 1 lt 1 lw 4 lc rgb "#820808"
set style line 2 lt 2 lw 4 lc rgb "#088208"

set xlabel "r"
set ylabel "Lennard-Jones potential"

set output filename.".png"
plot filename.'.dat' u 3:6 t "fully atomistic" w l ls 1, '' u 3:7 t "LOD" w l ls 2
EOF
	echo "-> Done."
done
