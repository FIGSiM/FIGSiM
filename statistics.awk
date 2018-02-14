####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

function round(d)
{
	return int(d+0.5)
}
function abs(val)
{
	if(val<0.0) return -val
	return val
}
function safe_sqrt(val)
{
	if(abs(val)<1.2E-7) return 0.0
	return sqrt(val)
}
function format_value(val, err)
{
	if(abs(err)>1.2E-7){
		factor=10
		while(int(err*factor+0.5)<10){
			factor=factor*10
		}
		factor=factor/10
		return round(factor*val)/factor" +/- "round(factor*err)/factor
	}
	return val" +/- 0"
}
BEGIN{
	racemic=0
	firstU=0
}
(!/^#/){
	if($4<0) Efield[$1]=-$4; else Efield[$1]=$4
	nr_datapoints[$1]=nr_datapoints[$1]+1
	if(firstU==0){ # for numerical precision (e^+/-100 can be a large number) shift energies by first value we encounter
		firstU=1
		dU=$7*1000
	}
	bw=exp(-($7*1000-dU)/(8.3144621*$2)) # e^(-U/RT) if U is in J/mol
	
	boltzmann_weight[$1]=boltzmann_weight[$1]+bw
	boltzmann_weight2[$1]=boltzmann_weight2[$1]+bw*bw
	
	density[$1]=density[$1]+bw*$5 ; density2[$1]=density2[$1]+bw*$5*$5
	
	N[$1]=N[$1]+bw*$6/1E20 ; N2[$1]=N2[$1]+bw*$6*$6/1E40 # number density in 10^20 molecules/cm3
	
	Utot[$1]=Utot[$1]+bw*$7 ; Utot2[$1]=Utot2[$1]+bw*$7*$7
	
	Cvp[$1]=Cvp[$1]+bw*$8 ; Cvp2[$1]=Cvp2[$1]+bw*$8*$8
	
	epsilon[$1]=epsilon[$1]+bw*$9 ; epsilon2[$1]=epsilon2[$1]+bw*$9*$9
	
	cs[$1]=cs[$1]+bw*$10 ; cs2[$1]=cs2[$1]+bw*$10*$10
	
	cs_2[$1]=cs_2[$1]+bw*$12 ; cs_22[$1]=cs_22[$1]+bw*$12*$12
	
	cs_3[$1]=cs_3[$1]+bw*$14 ; cs_32[$1]=cs_32[$1]+bw*$14*$14
}
END{
	for(dd in cs){
		var_div=1.0
		nrsims=nr_datapoints[dd]
		if(nrsims>1){
			print dd "\t" nrsims "\t" Efield[dd] "\t" format_value(density[dd]/boltzmann_weight[dd],safe_sqrt((var_div*density2[dd]*boltzmann_weight[dd]-density[dd]*density[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(N[dd]/boltzmann_weight[dd],safe_sqrt((var_div*N2[dd]*boltzmann_weight[dd]-N[dd]*N[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(Utot[dd]/boltzmann_weight[dd],safe_sqrt((var_div*Utot2[dd]*boltzmann_weight[dd]-Utot[dd]*Utot[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(Cvp[dd]/boltzmann_weight[dd],safe_sqrt((var_div*Cvp2[dd]*boltzmann_weight[dd]-Cvp[dd]*Cvp[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(epsilon[dd]/boltzmann_weight[dd],safe_sqrt((var_div*epsilon2[dd]*boltzmann_weight[dd]-epsilon[dd]*epsilon[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(cs[dd]/boltzmann_weight[dd],safe_sqrt((var_div*cs2[dd]*boltzmann_weight[dd]-cs[dd]*cs[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(cs_2[dd]/boltzmann_weight[dd],safe_sqrt((var_div*cs_22[dd]*boltzmann_weight[dd]-cs_2[dd]*cs_2[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(cs_3[dd]/boltzmann_weight[dd],safe_sqrt((var_div*cs_32[dd]*boltzmann_weight[dd]-cs_3[dd]*cs_3[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))) "\t" format_value(cs_3[dd]/boltzmann_weight[dd]*N[dd]/boltzmann_weight[dd],safe_sqrt((cs_3[dd]/boltzmann_weight[dd])**2*((var_div*N2[dd]*boltzmann_weight[dd]-N[dd]*N[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))+(N[dd]/boltzmann_weight[dd])**2*((var_div*cs_32[dd]*boltzmann_weight[dd]-cs_3[dd]*cs_3[dd])/(boltzmann_weight[dd]*boltzmann_weight[dd]-boltzmann_weight2[dd]))))
		} else{
			print dd "\t1\t" Efield[dd] "\t" density[dd]/boltzmann_weight[dd] "\t-\t" N[dd]/boltzmann_weight[dd] "\t-\t" Utot[dd]/boltzmann_weight[dd] "\t-\t" Cvp[dd]/boltzmann_weight[dd] "\t-\t" epsilon[dd]/boltzmann_weight[dd] "\t-\t" cs[dd]/boltzmann_weight[dd] "\t-\t" cs_2[dd]/boltzmann_weight[dd] "\t-\t" cs_3[dd]/boltzmann_weight[dd] "\t-"
		}
	}
}

