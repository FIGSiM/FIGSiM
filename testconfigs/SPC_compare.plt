#!/usr/bin/gnuplot

#######################################
# Correlation function GnuPlot script #
#  created Sep 2011, Andreas Tillack  #
#######################################

SPC='SPC-LOD_0_analysis_1'
LOD='SPC-LOD_1_analysis_1'
EXP='testconfigs/Water_experimental'

set term pngcairo dashed enhanced font "Arial,36" size 1600,1200
set samples 2000
set border lw 4

set style histogram clustered gap 1
set style fill solid border -1

set style line 1 lt 1 lw 4 lc rgb "#820808"
set style line 2 lt 2 lw 4 lc rgb "#088208"

set xlabel "r (\305)"
set ylabel "g(r)"

kT=298.0*1.38E-4

set output "SPC H2O O-O radial distribution function.png"
plot [0:9] [0:] SPC.'.dat' u 2:3 t "SPC" smooth csplines w l lw 4 lc rgb "#820808", LOD.'.dat' u 2:3 t "LOD" smooth csplines w l lw 4 lc rgb "#080882", EXP.'.dat' u 1:2:3 t "Experimental" smooth csplines w l lw 4 lc rgb "#8208FF"

