#!/usr/bin/gnuplot
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

#######################################
# Correlation function GnuPlot script #
#  created Sep 2011, Andreas Tillack  #
#######################################

filename='acetonitrile-one_0'

set term png enhanced font Arial 14 size 1600,1200
set samples 500

boltzmann(x) = (x>d_b) ? a_b*sqrt(2.0/(pi*(b_b**6)))*(x-d_b)**2*exp(-(x-d_b)**2/(2*b_b**2)) : 0.0
gauss(x) = a_g*exp(-(x-d_g)**2/(2*b_g**2))
decay(x) = a_g*exp(-(x-d_g)/b_g)


set output filename.".statistics.volume.png"
plot filename.'.stat.dat' u 1:2 t "Volume" w l lw 2

set output filename.".statistics.volume-logscale.png"
set logscale y
plot filename.'.stat.dat' u 1:2 t "Volume" w l lw 2
unset logscale

set output filename.".statistics.M_x.png"
plot filename.'.stat.dat' u 1:4 t "<M_x>" w l lw 2

set output filename.".statistics.M_y.png"
plot filename.'.stat.dat' u 1:5 t "<M_y>" w l lw 2

set output filename.".statistics.M_z.png"
plot filename.'.stat.dat' u 1:6 t "<M_z>" w l lw 2

set output filename.".statistics.M2.png"
plot filename.'.stat.dat' u 1:10 t "<M^2>" w l lw 2

set output filename.".statistics.varM.png"
plot filename.'.stat.dat' u 1:12 t "var(M)" w l lw 2

set output filename.".statistics.epsilon_comparison.png"
plot filename.'.stat.dat' u 1:14 t "Epsilon" w l lw 2, '' u 1:16 t "Epsilon_{RF}" w l lw 2

set output filename.".statistics.heatcapacity.png"
plot filename.'.stat.dat' u 1:18 t "Heat capacity" w l lw 2

set output filename.".statistics.potentials_ES_LJ.png"
set y2tics
plot filename.'.stat.dat' u 1:24 t "U_{ES}" w l lw 2, '' u 1:26 t "U_{LJ}" w l lw 2 axis x1y2

set output filename.".statistics.potentials_muE_G.png"
plot filename.'.stat.dat' u 1:22 t "U_{muE}" w l lw 2, '' u 1:28 t "U_{G}" w l lw 2 axis x1y2

set output filename.".statistics.potentials_total.png"
unset y2tics
plot filename.'.stat.dat' u 1:20 t "U_{total}" w l lw 2

a_g=1.0
d_g=0.0
b_g=100

set fit errorvariables
fit decay(x) filename.'.stat.autocorr.dat' via b_g

set output filename.".statistics.autocorr.png"
plot filename.'.stat.autocorr.dat' w l lw 2, decay(x) t sprintf("Fit (decay constant=%.3f +/- %.3f)",b_g,b_g_err) w l lw 2

binwidth=0.2
bin(x,width)=width*floor(x/width)+width/2.0
set boxwidth binwidth

set output filename.".statistics.epsilon.png"

set x2tics
set x2label "Step"
set y2label "Epsilon"
set y2tics

set style histogram clustered gap 1
set style fill solid border -1

set xlabel "Epsilon"
set ylabel "Counts"

set table 'histogram.dat'
plot filename.'.stat.dat' u (bin($14,binwidth)):(1.0) smooth freq
unset table

d_g=10
a_g=1
b_g=2

set fit errorvariables
fit gauss(x) 'histogram.dat' via a_g,b_g,d_g

d_b=d_g
a_b=a_g/d_g**2
b_b=b_g
fit boltzmann(x) 'histogram.dat' via a_b,b_b,d_b

print ""
print "Average epsilon from Maxwell-Boltzmann distribution fit:"
print sprintf("%.3f +/- %.3f",sqrt(8.0/pi*b_b**2)+d_b,sqrt(b_b_err**2+d_b_err**2))
print ""
print "Average epsilon from Gaussian fit:"
print sprintf("%.3f +/- %.3f",d_g,d_g_err)
print ""

plot 'histogram.dat' t "Epsilon histogram" w boxes lc rgb "#820808", boltzmann(x) t sprintf("Maxwell-Boltzmann distribution fit (<eps> = %.3f +/- %.3f)",sqrt(8.0/pi*b_b**2)+d_b,sqrt(b_b_err**2+d_b_err**2)) w l lw 2 lc rgb "#080882", gauss(x) t sprintf("Gaussian fit (<eps> = %.3f +/- %.3f)",d_g,d_g_err) w l lw 2 lc rgb "#088208", filename.'.stat.dat' u 1:14 t "Epsilon" w l lw 2 lc rgb "#488292" axis x2y2

set output filename.".statistics.epsilon_s.png"

set y2label "Epsilon_{RF}"

set table 'histogram.dat'
plot filename.'.stat.dat' u (bin($16,binwidth)):(1.0) smooth freq
unset table

set fit errorvariables
fit boltzmann(x) 'histogram.dat' via a_b,b_b,d_b
fit gauss(x) 'histogram.dat' via a_g,b_g,d_g

print ""
print "Average epsilon from Maxwell-Boltzmann distribution fit:"
print sprintf("%.3f +/- %.3f",sqrt(8.0/pi*b_b**2)+d_b,sqrt(b_b_err**2+d_b_err**2))
print ""
print "Average epsilon from Gaussian fit:"
print sprintf("%.3f +/- %.3f",d_g,d_g_err)
print ""

plot 'histogram.dat' t "Epsilon_{RF} histogram" w boxes lc rgb "#820808", boltzmann(x) t sprintf("Maxwell-Boltzmann distribution fit (<eps> = %.3f +/- %.3f)",sqrt(8.0/pi*b_b**2)+d_b,sqrt(b_b_err**2+d_b_err**2)) w l lw 2 lc rgb "#080882", gauss(x) t sprintf("Gaussian fit (<eps> = %.3f +/- %.3f)",d_g,d_g_err) w l lw 2 lc rgb "#088208", filename.'.stat.dat' u 1:14 t "Epsilon" w l lw 2 lc rgb "#488292" axis x2y2
