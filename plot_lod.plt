#!/usr/bin/gnuplot
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

filename="Acetonitrile.LOD1_1"

set datafile separator "\t"

set term pngcairo dashed enhanced font "Arial,36" size 1600,1200
set samples 1000
set border lw 4

r_diff_x = system("awk -F'\\t' '(\$1==\"#r_diff:\"){print \$2};(NR>16){exit}' ".filename."_x.dat")
r_diff_y = system("awk -F'\\t' '(\$1==\"#r_diff:\"){print \$2};(NR>16){exit}' ".filename."_y.dat")
r_diff_z = system("awk -F'\\t' '(\$1==\"#r_diff:\"){print \$2};(NR>16){exit}' ".filename."_z.dat")
avg_width = system("awk -F'\\t' '(\$1==\"#avg_width:\"){print \$2};(NR>16){exit}' ".filename."_x.dat")

set output filename."_interactions.png"
set xlabel "Distance [Angström]"
set ylabel "LJ potential energy [perg]"
plot [avg_width:2*8*avg_width] [-0.5:0.3] filename.'_interactions_x.dat' u 1:2 t "fully atomistic" w l lt 2 lw 4 lc rgb "#820808", '' u 1:4 t "best fit epsilon" w l lt 1 lw 4 lc rgb "#080882", '' u 1:6 t "textured epsilon" w l lt 1 lw 4 lc rgb "#088208",\
                       filename.'_interactions_y.dat' u 1:2 not w l lt 2 lw 4 lc rgb "#820808", '' u 1:4 not w l lt 1 lw 4 lc rgb "#080882", '' u 1:6 t "y-y" w l lt 1 lw 4 lc rgb "#088208",\
                       filename.'_interactions_z.dat' u 1:2 not w l lt 2 lw 4 lc rgb "#820808", '' u 1:4 not w l lt 1 lw 4 lc rgb "#080882", '' u 1:6 t "z-z" w l lt 1 lw 4 lc rgb "#088208"

set output filename."_potentials.png"
set xlabel "Distance [Angström]"
set ylabel "LJ potential energy [perg]"
plot [avg_width:8*avg_width] [-0.6:0.2] filename.'_x.dat' u 1:2 t "fully atomistic" w l lt 2 lw 4 lc rgb "#820808", '' u ($1+r_diff_x):3 t "best fit epsilon" w l lt 1 lw 4 lc rgb "#080882", '' u ($1+r_diff_x):4 t "textured epsilon" w l lt 1 lw 4 lc rgb "#088208",\
                       filename.'_y.dat' u 1:2 not w l lt 2 lw 4 lc rgb "#820808", '' u ($1+r_diff_y):3 not w l lt 1 lw 4 lc rgb "#080882", '' u ($1+r_diff_y):4 not w l lt 1 lw 4 lc rgb "#088208",\
                       filename.'_z.dat' u 1:2 not w l lt 2 lw 4 lc rgb "#820808", '' u ($1+r_diff_z):3 not w l lt 1 lw 4 lc rgb "#080882", '' u ($1+r_diff_z):4 not w l lt 1 lw 4 lc rgb "#088208"

# obtain sqrt(eps) fit coefficients
sigma_LOD=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$2};(NR>2){exit}' ".filename.".epsilon.dat")))
series_6=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$3};(NR>2){exit}' ".filename.".epsilon.dat")))
series_5=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$4};(NR>2){exit}' ".filename.".epsilon.dat")))
series_4=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$5};(NR>2){exit}' ".filename.".epsilon.dat")))
series_3=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$6};(NR>2){exit}' ".filename.".epsilon.dat")))
series_2=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$7};(NR>2){exit}' ".filename.".epsilon.dat")))
series_1=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$8};(NR>2){exit}' ".filename.".epsilon.dat")))
series_0=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$9};(NR>2){exit}' ".filename.".epsilon.dat")))

a(x)=series_0*x**6+series_1*x**5+series_2*x**4+series_3*x**3+series_4*x**2+series_5*x+series_6
model(x) = a(x)/(x+sigma_LOD)**6

FIT_LIMIT=1E-12
#fit model(x) filename.'.epsilon.dat' u 1:2 via sigma_LOD,series_6,series_5,series_4,series_3,series_2,series_1

set output filename."_epsilon.png"
set lmargin 10
set rmargin 10
set key right bot

set multiplot
set xlabel "{r_{test} [Angström]"
set y2label "relative error %"
set ylabel ""
set y2tics
set noytics
set size 1,0.4
set origin 0.0,0.0
set bmargin 3
set tmargin 0
plot filename.'.epsilon.dat' u 1:((model($1)-$2)/$2*100) t "MC fit relative error %" w l lw 4 lc rgb "#080882" axis x1y2 # , '' u 1:($3/$2*100) t "epsilon std.dev. %" w l lw 2 lc rgb "#088208" axis x1y2, '' u 1:(-$3/$2*100) not w l lw 2 lc rgb "#088208" axis x1y2
set size 1,0.6
set origin 0.0,0.4
set ytics
set format x ""
set xlabel ""
set ylabel "<{/Symbol \326e}>"
set y2label ""
set ytics
set noy2tics
set bmargin 0
set tmargin 1
plot filename.'.epsilon.dat' u 1:2 t "{/Symbol \326e}" w l lw 4 lc rgb "#820808", model(x) t "MC fit" w l lw 4 lc rgb "#088208"
set nomultiplot

reset
set term pngcairo dashed enhanced font "Arial,48" size 1600,1200
set samples 1000
set border lw 4
set output filename."_epsilon_zoom.png"
set lmargin 10
set rmargin 10
set key right bot

set multiplot
set xlabel ""
unset y2label
unset ylabel
set y2tics
set noytics
set size 1,0.4
set origin 0.0,0.0
set bmargin 3
set tmargin 0
plot [0:8] filename.'.epsilon.dat' u 1:((model($1)-$2)/$2*100) not w l lw 4 lc rgb "#080882" axis x1y2 # , '' u 1:($3/$2*100) t "epsilon std.dev. %" w l lw 2 lc rgb "#088208" axis x1y2, '' u 1:(-$3/$2*100) not w l lw 2 lc rgb "#088208" axis x1y2
set size 1,0.6
set origin 0.0,0.4
set ytics
set format x ""
set xlabel ""
set ylabel ""
set y2label ""
set ytics
set noy2tics
set bmargin 0
set tmargin 1
plot [0:8] filename.'.epsilon.dat' u 1:2 not w l lw 4 lc rgb "#820808", model(x) not w l lw 4 lc rgb "#088208"
set nomultiplot
