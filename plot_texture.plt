#!/usr/bin/gnuplot
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

filename="CLD.LOD2_1"

set datafile separator "\t"

set term pngcairo dashed enhanced font "Arial,36" size 1600,1200
set samples 1000
set border lw 4

k=1.3806488E-4 # perg/K

set output filename.".texture.png"
set xlabel "cos({/Symbol q})"
set ylabel "fraction"

step=1.0
ystep=1.0

xmin=0.685
ymin=0.6

xmax=1.45
ymax=1.51

#a_h=0.12
#b_h=(ymax-ystep-a_h*(xmax**2-step**2))/(xmax-step)
#c_h=ystep-a_h*step**2-b_h*step
#a_l=(ymin-ystep-2*a_h*step*(xmin-step)-b_h*(xmin-step))/(xmin-step)**2

step=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$2};(NR>2){exit}' ".filename.".IA_vs_texture.dat")))
c_h=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$3};(NR>2){exit}' ".filename.".IA_vs_texture.dat")))
b_h=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$4};(NR>2){exit}' ".filename.".IA_vs_texture.dat")))
a_h=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$5};(NR>2){exit}' ".filename.".IA_vs_texture.dat")))
a_l=real(system(sprintf("awk -F'\\t' '(\$1==\"#coefficients:\"){print \$6};(NR>2){exit}' ".filename.".IA_vs_texture.dat")))

adj_square(x) = (x>step) ? a_h*x**2+b_h*x+c_h : a_l*x**2+(2*step*(a_h-a_l)+b_h)*x+c_h-step**2*(a_h-a_l)
#fit adj_square(x) filename.'.IA_vs_texture.dat' u 3:4 via a_h,b_h,c_h,a_l,step

plot filename.'.IA_vs_texture.dat' every 256::0 u 1:4 t "texture ({/Symbol f}=0)" w l lt 1 lw 4 lc rgb "#820808", '' every 256::64 u 1:4 t "texture ({/Symbol f}={/Symbol p}/2)" w l lt 1 lw 4 lc rgb "#820882", '' every 256::128 u 1:4 t "texture ({/Symbol f}={/Symbol p})" w l lt 1 lw 4 lc rgb "#088282", '' every 256::0 u 1:3 t "IA ({/Symbol f}=0)" w l lt 2 lw 4 lc rgb "#088208", '' every 256::0 u 1:(adj_square($3)) t "piecewise square ({/Symbol f}=0)" w l lt 2 lw 4 lc rgb "#080808", '' every 256::64 u 1:3 t "IA ({/Symbol f}={/Symbol p}/2)" w l lt 2 lw 4 lc rgb "#088208", '' every 256::64 u 1:(adj_square($3)) t "piecewise square ({/Symbol f}={/Symbol p}/2)" w l lt 2 lw 4 lc rgb "#080808"

set output filename.".IA_vs_texture.png"
set xlabel "IA"
set ylabel "texture"
plot filename.'.IA_vs_texture.dat' u 3:4 t "texture ({/Symbol f}=0)" w p pt 7 ps 2 lc rgb "#088208", x t "IA" w l lw 4 lc rgb "#280882", adj_square(x) t sprintf("fit (x_0=%.2f, a_h=%.2f, b_h=%.2f, c_h=%.2f, a_l=%.2f)",step,a_h,b_h,c_h,a_l) w l lt 2 lw 4 lc rgb "#080808"
