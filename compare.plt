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

filename='H2O_0_analysis_1'
plot_nr=1
ga_jump=0

filename2='H2O_1_analysis_1'
plot_nr2=1
ga_jump2=0

auto_axis=0
divide=0
gammaalpha=0
Qtensor=0
if(divide>0) auto_axis=0

set term png enhanced font Arial 14 size 1600,1200
set samples 1000

# Get variables from data file (limit ourselves to first 8 lines)
# Title of plot (with original filename, have to substitute "_" with "\_" for gnuplot not to turn it into subscripts)
PlotTitle = system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename.".dat")." (".system("echo \"".filename."\" | sed s/_/\\\\\\\\_/g").") "
PlotTitle2 = system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename2.".dat")." (".system("echo \"".filename2."\" | sed s/_/\\\\\\\\_/g").") "
if(divide>0) PlotTitle=system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename.".dat")."\n divided by ".system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Results per frame
NpF = system("awk -F'\\t' '(\$1==\"#N per frame:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Number of plots in data file
nr_plots = system("awk -F'\\t' '(\$1==\"#N correlations:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Total number of data points
Ntotal = system("awk -F'\\t' '(\$1==\"#N total:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Get axis labels
x_label = system(sprintf("awk -F'\\t' '(\$1==\"#x_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
if(x_label eq "cos(alpha)") gammaalpha=gammaalpha+1
if(x_label eq "lambda_1") Qtensor=Qtensor+38
if(gammaalpha>0) print "Gammaalpha plot comparison not supported."; exit
x_label2 = system(sprintf("awk -F'\\t' '(\$1==\"#x_label:\"){print \$%i};(NR>16){exit}' ".filename2.".dat",plot_nr2+ga_jump2+1))
if(x_label2 eq "cos(alpha)") gammaalpha=gammaalpha+1
if(gammaalpha>0) print "Gammaalpha plot comparison not supported."; exit
y_label = system(sprintf("awk -F'\\t' '(\$1==\"#y_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
y_label2 = system(sprintf("awk -F'\\t' '(\$1==\"#y_label:\"){print \$%i};(NR>16){exit}' ".filename2.".dat",plot_nr2+ga_jump2+1))
# Minimum x value
x_min = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_min:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)))
# Maximum x value
x_max = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_max:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)))
# Minimum x value
x_min2 = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_min:\"){print \$%i};(NR>16){exit}' ".filename2.".dat",plot_nr2+ga_jump2+1)))
# Maximum x value
x_max2 = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_max:\"){print \$%i};(NR>16){exit}' ".filename2.".dat",plot_nr2+ga_jump2+1)))

print "-> Plotting ".y_label." vs ".y_label2." ..."

set angle degrees

set style line 1 lt 1 lw 4 lc rgb "#820808"
set style line 2 lt 2 lw 4 lc rgb "#088208"

set xlabel x_label
if(Qtensor>0) set xlabel "Steps"; else set xrange [x_min:x_max]
if(divide>0) y_label=y_label." / ".y_label2
set ylabel y_label
set macros
axis_to_use="x1y1"
if((y_label ne y_label2) && (auto_axis>0)) axis_to_use="x1y2"; set y2tics; set y2label y_label2; set ylabel textcolor lt 1; set y2label textcolor lt 2; if(x_label ne x_label2) axis_to_use="x2y2"; set x2label x_label2; set y2label y_label2; set x2range [x_min2:x_max2]; set y2tics; set x2tics; set ylabel textcolor lt 1; set y2label textcolor lt 2; set xlabel textcolor lt 1; set x2label textcolor lt 2;
if((x_label ne x_label2) && (auto_axis>0)) if(y_label eq y_label2) axis_to_use="x2y1"; set x2label x_label2; set x2range [x_min2:x_max2]; set x2tics; set xlabel textcolor lt 1; set x2label textcolor lt 2;

if(divide>0) outfile=filename."_plot".plot_nr."_over_".filename2."_plot".plot_nr2; else outfile=filename."_plot".plot_nr."_vs_".filename2."_plot".plot_nr2
outfile=system(sprintf("echo \"%s\" | sed s/\\\\//\\_/g | sed s/\\\\.//g",outfile))
if(Qtensor==0) set output outfile.".png"
if((divide==0) && (Qtensor==0)) plot filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t PlotTitle smooth csplines w l ls 1, filename2.'.dat' u (column(2*plot_nr2+ga_jump2)):(column(2*plot_nr2+ga_jump2+1)) t PlotTitle2 smooth csplines w l ls 2 axis @axis_to_use; else if(Qtensor==0) plot filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)/column(2*plot_nr2+ga_jump2+1)) t PlotTitle w l ls 1

if(Qtensor>0) set output outfile."_lambda.png"
if(Qtensor>0) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t PlotTitle."{/Symbol l}_1" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1))  notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+2)):(column(2*plot_nr+ga_jump+3)) t PlotTitle."{/Symbol l}_2" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+2)):(column(2*plot_nr+ga_jump+3)) notitle smooth csplines lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+4)):(column(2*plot_nr+ga_jump+5)) t PlotTitle."{/Symbol l}_3" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#080882", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+4)):(column(2*plot_nr+ga_jump+5)) notitle smooth csplines lw 2 lc rgb "#080882", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2)):(column(2*plot_nr2+ga_jump2+1)) t PlotTitle2."{/Symbol l}_1" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#FF8200", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2)):(column(2*plot_nr2+ga_jump2+1))  notitle smooth csplines lw 2 lc rgb "#FF8200", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+2)):(column(2*plot_nr2+ga_jump2+3)) t PlotTitle2."{/Symbol l}_2" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#82FF00", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+2)):(column(2*plot_nr2+ga_jump2+3)) notitle smooth csplines lw 2 lc rgb "#82FF00", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+4)):(column(2*plot_nr2+ga_jump2+5)) t PlotTitle2."{/Symbol l}_3" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#8200FF", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+4)):(column(2*plot_nr2+ga_jump2+5)) notitle smooth csplines lw 2 lc rgb "#8200FF"
if(Qtensor>0) set output outfile."_P1.png"
if(Qtensor>0) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+24)):(column(2*plot_nr+ga_jump+25)) t PlotTitle."P_1(v_1)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+24)):(column(2*plot_nr+ga_jump+25)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+26)):(column(2*plot_nr+ga_jump+27)) t PlotTitle."P_1(v_2)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+26)):(column(2*plot_nr+ga_jump+27)) notitle smooth csplines lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+28)):(column(2*plot_nr+ga_jump+29)) t PlotTitle."P_1(v_3)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#080882", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+28)):(column(2*plot_nr+ga_jump+29)) notitle smooth csplines lw 2 lc rgb "#080882", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+24)):(column(2*plot_nr2+ga_jump2+25)) t PlotTitle2."P_1(v_1)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#FF8200", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+24)):(column(2*plot_nr2+ga_jump2+25)) notitle smooth csplines lw 2 lc rgb "#FF8200", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+26)):(column(2*plot_nr2+ga_jump2+27)) t PlotTitle2."P_1(v_2)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#82FF00", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+26)):(column(2*plot_nr2+ga_jump2+27)) notitle smooth csplines lw 2 lc rgb "#82FF00", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+28)):(column(2*plot_nr2+ga_jump2+29)) t PlotTitle2."P_1(v_3)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#8200FF", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+28)):(column(2*plot_nr2+ga_jump2+29)) notitle smooth csplines lw 2 lc rgb "#8200FF"
if(Qtensor>0) set output outfile."_SR.png"
if(Qtensor>0) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+30)):(column(2*plot_nr+ga_jump+31)) t PlotTitle."S" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+30)):(column(2*plot_nr+ga_jump+31)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+32)):(column(2*plot_nr+ga_jump+33)) t PlotTitle."R" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+32)):(column(2*plot_nr+ga_jump+33)) notitle smooth csplines lw 2 lc rgb "#088208", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+30)):(column(2*plot_nr2+ga_jump2+31)) t PlotTitle2."S" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#FF8200", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+30)):(column(2*plot_nr2+ga_jump2+31)) notitle smooth csplines lw 2 lc rgb "#FF8200", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+32)):(column(2*plot_nr2+ga_jump2+33)) t PlotTitle2."R" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#82FF00", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+32)):(column(2*plot_nr2+ga_jump2+33)) notitle smooth csplines lw 2 lc rgb "#82FF00"
if(Qtensor>0) set output outfile."_Pmax.png"
if(Qtensor>0) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+34)):(column(2*plot_nr+ga_jump+35)) t PlotTitle."P_{max}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+34)):(column(2*plot_nr+ga_jump+35)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+36)):(column(2*plot_nr+ga_jump+37)) t PlotTitle."P_{dipole}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+36)):(column(2*plot_nr+ga_jump+37)) notitle smooth csplines lw 2 lc rgb "#088208", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+34)):(column(2*plot_nr2+ga_jump2+35)) t PlotTitle2."P_{max}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#FF8200", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+34)):(column(2*plot_nr2+ga_jump2+35)) notitle smooth csplines lw 2 lc rgb "#FF8200", filename2.'.dat' u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+36)):(column(2*plot_nr2+ga_jump2+37)) t PlotTitle2."P_{dipole}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#82FF00", "" u (column(2*plot_nr2+ga_jump2-1)):(column(2*plot_nr2+ga_jump2+36)):(column(2*plot_nr2+ga_jump2+37)) notitle smooth csplines lw 2 lc rgb "#82FF00"

set angle degrees

set xlabel x_label
if(Qtensor>0) set xlabel "Steps"
set ylabel y_label
if(gammaalpha>0) set zlabel z_label; set dgrid3d nr_bins,nr_bins; set zlabel rotate by -270; set hidden3d; set contour
set style line 8 lt 1 lw 3 lc rgb "#020208"

if(Qtensor==0) set output filename."-correlation_".plot_nr.".png"

if(Qtensor>0) set view equal xyz
if(Qtensor>0) set parametric
if(Qtensor>0) set isosample 13
if(Qtensor>0) set view 60, 136, 1.22, 1.26
if(Qtensor>0) set urange [-90.0:90.0] noreverse nowriteback
if(Qtensor>0) set vrange [0:360.0] noreverse nowriteback
if(Qtensor>0) set output outfile."_vs.png"
if(Qtensor>0) set ticslevel 0.0
if(Qtensor>0) set xlabel "x"
if(Qtensor>0) set ylabel "y"
if(Qtensor>0) set zlabel "z"

set style line 8 lt 1 lw 3 lc rgb "#020208"
set xlabel x_label
if(Qtensor>0) set xlabel "Steps"
set ylabel y_label
if(Qtensor>0) splot cos(u)*cos(v),cos(u)*sin(v),sin(u) notitle w l lt 3 lw 1 lc rgb "#2882ED", filename.'.dat' u (column(2*plot_nr+ga_jump+6)):(column(2*plot_nr+ga_jump+8)):(column(2*plot_nr+ga_jump+10)) t "Major director axis ".PlotTitle w p lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", filename2.'.dat' u (column(2*plot_nr+ga_jump+6)):(column(2*plot_nr+ga_jump+8)):(column(2*plot_nr+ga_jump+10)) t "Major director axis ".PlotTitle2 w p lt 1 pt 6 ps 2 lw 2 lc rgb "#088208"


print "<- Done."
