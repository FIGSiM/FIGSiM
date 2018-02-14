#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

echo "# name	T [K]	p [atm]	E [V/µm]	density [g/cm³]	N [molecules/cm³]	U_total [kJ/mol]	heat capacity [J/mol-K]	dielectric	cos	cos²	cos³" > results_poled.dat
echo "# name	T [K]	p [atm]	E [V/µm]	density [g/cm³]	N [molecules/cm³]	U_total [kJ/mol]	heat capacity [J/mol-K]	dielectric	cos	cos²	cos³" > results_unpoled.dat
for datafile in `ls *.stat`; do
	datafile=${datafile:0:${#datafile}-5}
	conffile=${datafile:0:${#datafile}-2}.conf
	if [ ! -e ${conffile} ]; then
		conffile=${datafile:0:${#datafile}-3}.conf
		if [ ! -e ${conffile} ]; then
			conffile=${datafile:0:${#datafile}-4}.conf # go up to hundred subruns
		fi
	fi
	if [ ! -e ${conffile} ]; then
		conffile=$datafile.traj
		if [ ! -e ${conffile} ]; then
			conffile=$datafile.traj.bz2
			if [ -e ${conffile} ]; then
				head -c 512000 ${conffile} | bunzip2 | head -n 10 > $datafile.conf.traj
				conffile=$datafile.conf.traj
			fi
		fi
	fi
	if [ ! -e ${conffile} ]; then
		echo "Cannot find configuration file."
		exit
	fi
	name=$(MCparams ${conffile} "Simulation Parameters" "FileOut")
	Efield=$(MCparams ${conffile} "Simulation Parameters" "Efield")
	Efield_zero=$(echo "${Efield}" | awk '{printf("%d\n",$1+0.5)}')
	outname=results_unpoled.dat
	if [ ${Efield_zero} -ne 0 ]; then
		outname=results_poled.dat
	fi
	p="NVT"
	if [ $(MCparams ${conffile} "Simulation Parameters" "NpT") -gt 0 ]; then
		p=$(MCparams ${conffile} "Simulation Parameters" "pext")
	fi
	density=`awk -F'\t' '(($1=="Average Density:")){printf $2}' ${datafile}.stat`
	N=`awk -F'\t' '(($1=="Number Density:")){printf $2}' ${datafile}.stat`
	Utot_string="U_total"
	Utot=`awk -F'\t' '(($1=="<U>")){ Utot=$2 }END{ print Utot }' ${datafile}.stat`
	if [ "${Utot}" == "" ]; then
		Utot_string="H_total"
		Utot=`awk -F'\t' '(($1=="<H>")){ Utot=$2 }END{ print Utot }' ${datafile}.stat`
	fi
	Cvp=`awk -F'\t' '(($1=="<Cv>") || ($1=="<Cp>")){printf $2}' ${datafile}.stat`
	
	density_nr=`echo "${density}*10/1" | bc`
	if [ "${p}" == "NVT" ]; then # pass through user specified volumes
		density_nr=8
	fi
	if [ $density_nr -ge 5 ]; then
		dielectric=`awk -F'\t' '(($1=="<epsilon>")){printf $2}' ${datafile}.stat`
		n2=$(MCparams ${conffile} "Simulation Parameters" "n2")
		T=$(MCparams ${conffile} "Simulation Parameters" "T")
		good_dielectric=$(echo "scale=3; (${dielectric}>${n2})" | bc)
		if [ $good_dielectric -gt 0 ]; then
			echo -e -n "${name}\t" >> ${outname}
			echo -e -n "${T}\t" >> ${outname}
			echo -e -n "${p}\t" >> ${outname}
			echo -e -n "${Efield}\t" >> ${outname}
			echo -e -n "${density}\t" >> ${outname}
			echo -e -n "${N}\t" >> ${outname}
			echo -e -n "${Utot}\t" >> ${outname}
			echo -e -n "${Cvp}\t" >> ${outname}
			echo -e -n "${dielectric}\t" >> ${outname}
			fn=$(ls -1 ${datafile}_Analysis_*.dat 2>/dev/null)
			if [ ! -e "${fn}" ]; then
				fn="${datafile}_analysis_1.dat"
			fi
			if [ -e "${fn}" ]; then
				awk -F'\t' '(($1=="#averages:")){printf $2}' ${fn} >> ${outname}
				echo -e -n "\t" >> ${outname}
				awk -F'\t' '(($1=="#stddevs:")){printf $2}' ${fn} >> ${outname}
				echo -e -n "\t" >> ${outname}
				awk -F'\t' '(($1=="#averages:")){printf $3}' ${fn} >> ${outname}
				echo -e -n "\t" >> ${outname}
				awk -F'\t' '(($1=="#stddevs:")){printf $3}' ${fn} >> ${outname}
				echo -e -n "\t" >> ${outname}
				awk -F'\t' '(($1=="#averages:")){printf $4}' ${fn} >> ${outname}
				echo -e -n "\t" >> ${outname}
				awk -F'\t' '(($1=="#stddevs:")){printf $4}' ${fn} >> ${outname}
			else
				echo "No analysis for ${datafile}."
			fi
			echo >> ${outname}
		else
			echo "Ignoring ${datafile} - dielectric is below n^2 (simulation got stuck)."
		fi
	else
		echo "Ignoring ${datafile} - density is below 0.5 g/cm^3."
	fi
done

echo "# name	# sims	E [V/µm]	density [g/cm³]	N [10^20 molecules/cm³]	${Utot_string} [kJ/mol]	heat capacity [J/mol-K]	dielectric	<cos>	<cos²>	<cos³>	N<cos³>" > tmp.dat
gawk -f statistics.awk results_unpoled.dat | sort -n >> tmp.dat
mv tmp.dat results_unpoled_overall.dat
echo "# name	# sims	E [V/µm]	density [g/cm³]	N [10^20 molecules/cm³]	${Utot_string} [kJ/mol]	heat capacity [J/mol-K]	dielectric	<cos>	<cos²>	<cos³>	N<cos³>" > tmp.dat
gawk -f statistics.awk results_poled.dat | sort -n >> tmp.dat
mv tmp.dat results_poled_overall.dat

