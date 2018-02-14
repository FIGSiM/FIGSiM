#!/usr/bin/gnuplot

#######################################
# Correlation function GnuPlot script #
#  created Sep 2011, Andreas Tillack  #
#######################################

filename='H2O_2_analysis_1'
fit_potential_of_mean_force=0
T=400

set term png enhanced font Arial 14 size 1600,1200
set samples 500

if(!exists("plot_nr")) plot_nr=1; gammaalpha=0; Qtensor=0; costheta=0;
ga_jump=0
if(gammaalpha>1) ga_jump=gammaalpha+Qtensor-1

# Get variables from data file (limit ourselves to first 8 lines)
# Title of plot (with original filename, have to substitute "_" with "\_" for gnuplot not to turn it into subscripts)
PlotTitle = system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename.".dat")."\n(".system("echo \"".filename."\" | sed s/_/\\\\\\\\_/g").")"
# Results per frame
NpF = system("awk -F'\\t' '(\$1==\"#N per frame:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Number of plots in data file
nr_plots = system("awk -F'\\t' '(\$1==\"#N correlations:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Total number of data points
Ntotal = system("awk -F'\\t' '(\$1==\"#N total:\"){print \$2};(NR>16){exit}' ".filename.".dat")
# Get axis labels
x_label = system(sprintf("awk -F'\\t' '(\$1==\"#x_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
if(x_label eq "cos(alpha)") gammaalpha=gammaalpha+1
if(x_label eq "cos(theta)") costheta=costheta+6
if(x_label eq "lambda_1") Qtensor=Qtensor+38
if(gammaalpha>0) y_label = system(sprintf("awk -F'\\t' '(\$1==\"#x_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+2)); else y_label = system(sprintf("awk -F'\\t' '(\$1==\"#y_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
z_label=y_label;
if(gammaalpha>0) z_label = system(sprintf("awk -F'\\t' '(\$1==\"#y_label:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
if(y_label eq "Normalized counts") average=system(sprintf("awk -F'\\t' '(\$1==\"#averages:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)); stddev=system(sprintf("awk -F'\\t' '(\$1==\"#stddevs:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1))
# Minimum x value
x_min = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_min:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)))
# Maximum x value
x_max = real(system(sprintf("awk -F'\\t' '(\$1==\"#x_max:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)))
# Histogram bin width
binwidth = real(system(sprintf("awk -F'\\t' '(\$1==\"#binwidth:\"){print \$%i};(NR>16){exit}' ".filename.".dat",plot_nr+ga_jump+1)))
frames = Ntotal/NpF
nr_bins=1
if(binwidth>0) nr_bins=int((x_max-x_min)/binwidth+0.5)

print "-> Plotting ".z_label." (".nr_bins." bins) ..."

set title PlotTitle

set boxwidth binwidth

set style histogram clustered gap 1
set style fill solid border -1

set angle degrees

set xlabel x_label
if(Qtensor) set xlabel "Steps"
set ylabel y_label
if(gammaalpha>0) set zlabel z_label; set dgrid3d nr_bins,nr_bins; set zlabel rotate by -270; set hidden3d; set contour
set style line 8 lt 1 lw 3 lc rgb "#020208"

if(Qtensor==0) set output filename."-correlation_".plot_nr.".png"

if(Qtensor) set output filename."-correlation_".plot_nr."_lambda.png"
if(Qtensor) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t "{/Symbol l}_1" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+2)):(column(2*plot_nr+ga_jump+3)) t "{/Symbol l}_2" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+2)):(column(2*plot_nr+ga_jump+3)) notitle smooth csplines lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+4)):(column(2*plot_nr+ga_jump+5)) t "{/Symbol l}_3" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#080882", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+4)):(column(2*plot_nr+ga_jump+5)) notitle smooth csplines lw 2 lc rgb "#080882"
if(Qtensor) set output filename."-correlation_".plot_nr."_P1.png"
if(Qtensor) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+24)):(column(2*plot_nr+ga_jump+25)) t "P_1(v_1)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+24)):(column(2*plot_nr+ga_jump+25)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+26)):(column(2*plot_nr+ga_jump+27)) t "P_1(v_2)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+26)):(column(2*plot_nr+ga_jump+27)) notitle smooth csplines lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+28)):(column(2*plot_nr+ga_jump+29)) t "P_1(v_3)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#080882", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+28)):(column(2*plot_nr+ga_jump+29)) notitle smooth csplines lw 2 lc rgb "#080882"
if(Qtensor) set output filename."-correlation_".plot_nr."_SR.png"
if(Qtensor) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+30)):(column(2*plot_nr+ga_jump+31)) t "S" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+30)):(column(2*plot_nr+ga_jump+31)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+32)):(column(2*plot_nr+ga_jump+33)) t "R" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+32)):(column(2*plot_nr+ga_jump+33)) notitle smooth csplines lw 2 lc rgb "#088208"
if(Qtensor) set output filename."-correlation_".plot_nr."_Pmax.png"
if(Qtensor) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+34)):(column(2*plot_nr+ga_jump+35)) t "P_{max}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+34)):(column(2*plot_nr+ga_jump+35)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+36)):(column(2*plot_nr+ga_jump+37)) t "P_{dipole}" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+36)):(column(2*plot_nr+ga_jump+37)) notitle smooth csplines lw 2 lc rgb "#088208"
if(Qtensor) set output filename."-correlation_".plot_nr."_cos3.png"
if(Qtensor) plot filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+38)):(column(2*plot_nr+ga_jump+39)) t "cos^3(v_1)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+38)):(column(2*plot_nr+ga_jump+39)) notitle smooth csplines lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+40)):(column(2*plot_nr+ga_jump+41)) t "cos^3(v_2)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+40)):(column(2*plot_nr+ga_jump+41)) notitle smooth csplines lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+42)):(column(2*plot_nr+ga_jump+43)) t "cos^3(v_3)" w yerr lt 1 pt 6 ps 2 lw 2 lc rgb "#080882", "" u (column(2*plot_nr+ga_jump-1)):(column(2*plot_nr+ga_jump+42)):(column(2*plot_nr+ga_jump+43)) notitle smooth csplines lw 2 lc rgb "#080882"
if(Qtensor+costheta) set view equal xyz
if(Qtensor) set parametric;
if(Qtensor) set isosample 13
if(Qtensor) set view 60, 136, 1.22, 1.26
if(Qtensor) set urange [-90.0:90.0] noreverse nowriteback
if(Qtensor) set vrange [0:360.0] noreverse nowriteback
if(Qtensor) set output filename."-correlation_".plot_nr."_vs.png"
if(Qtensor+costheta) set ticslevel 0.0
if(Qtensor+costheta) set xlabel "x"
if(Qtensor+costheta) set ylabel "y"
if(Qtensor+costheta) set zlabel "z"
if(costheta) set angle radians
if(Qtensor) splot cos(u)*cos(v),cos(u)*sin(v),sin(u) notitle w l lt 3 lw 1 lc rgb "#2882ED", filename.'.dat' u (column(2*plot_nr+ga_jump+6)):(column(2*plot_nr+ga_jump+8)):(column(2*plot_nr+ga_jump+10)) t "{v}_1" w p lt 1 pt 6 ps 2 lw 2 lc rgb "#820808", filename.'.dat' u (column(2*plot_nr+ga_jump+12)):(column(2*plot_nr+ga_jump+14)):(column(2*plot_nr+ga_jump+16)) t "{v}_2" w p lt 1 pt 6 ps 2 lw 2 lc rgb "#088208", filename.'.dat' u (column(2*plot_nr+ga_jump+18)):(column(2*plot_nr+ga_jump+20)):(column(2*plot_nr+ga_jump+22)) t "{v}_3" w p lt 1 pt 6 ps 2 lw 2 lc rgb "#080882"
if(gammaalpha) splot [-1.0:1.0] [x_min:x_max] filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)):(column(2*plot_nr+ga_jump+2)) notitle w pm3d; else 
if(costheta) splot filename.'.dat' u (column(2*plot_nr+ga_jump+2)*sqrt(1.0-column(2*plot_nr+ga_jump)**2)*cos(column(2*plot_nr+ga_jump+1))):(column(2*plot_nr+ga_jump+2)*sqrt(1.0-column(2*plot_nr+ga_jump)**2)*sin(column(2*plot_nr+ga_jump+1))):(column(2*plot_nr+ga_jump+2)*column(2*plot_nr+ga_jump)) t 'sdf(r)'; else 
if(y_label eq "Normalized counts") plot [x_min:x_max] [0:] filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t sprintf('binwidth: %.3f',binwidth) w boxes lc rgb "#820808"; else plot [x_min:x_max] filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t sprintf('binwidth: %.3f',binwidth) w l lt 1 lw 4 lc rgb "#820808"

a=1.0
b=1.0
c=0.1
d=3.5
xa=3.0
xm=18.0
k=1.3806488E-4 # Boltzmann constant in perg

model(x,a,b,c,d,xa,xm)=(exp((xa-x)/a)+exp((x-xm)/b))+c*x+d

set fit errorvariables
if((y_label eq "g(r)") && (fit_potential_of_mean_force>0)) fit model(x,a,b,c,d,xa,xm) filename.'.dat' u (column(2*plot_nr+ga_jump)):(-k*T*log(column(2*plot_nr+ga_jump+1))) via a,b,c,d,xa,xm

if((y_label eq "g(r)") && (fit_potential_of_mean_force>0)) set title system("awk -F'\\t' '(\$1==\"#Title:\"){print \$2};(NR>16){exit}' ".filename.".dat")."\n".sprintf('model fit: a=%.2f\261%.2f {\305}, b=%.2f\261%.2f {\305}, m=%.4f\261%.4f perg/{\305}, V_0=%.2f\261%.2f perg, l_a=%.1f\261%.1f {\305}, l_b=%.1f\261%.1f {\305}',a,a_err,b,b_err,c,c_err,d,d_err,xa,xa_err,xm,xm_err)
if((y_label eq "g(r)") && (fit_potential_of_mean_force>0)) set ylabel "-kT*ln(g(r)) [perg]"; set y2label "g(r)"; set y2tics; set output filename."-correlation_".plot_nr."-potential.png"; plot [x_min:xm+xa] [0:(model(xa+a+binwidth,a,b,c,d,xa,xm)+model(xm-b+binwidth,a,b,c,d,xa,xm))/2] filename.'.dat' u (column(2*plot_nr+ga_jump)):(-k*T*log(column(2*plot_nr+ga_jump+1))) t "potential of mean force" w p pt 7 ps 4 lc rgb "#820808" axis x1y1, model(x,a,b,c,d,xa,xm) t "model fit"  w l lt 1 lw 4 lc rgb "#088208" axis x1y1, filename.'.dat' u (column(2*plot_nr+ga_jump)):(column(2*plot_nr+ga_jump+1)) t "g(r)" w l lt 1 lw 4 lc rgb "#080882" axis x1y2

print "<- Done."

plot_nr=plot_nr+1
if(plot_nr<=nr_plots) reread

