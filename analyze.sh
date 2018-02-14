#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

MAX_NR_JOBS=1

# find out how many processors exist on this node (if available)
if [ -e "/proc/cpuinfo" ]; then
	nr_cpus=$(cat /proc/cpuinfo | grep "processor" | tail -n 1 | awk "(\$1==\"processor\"){print \$3+1}")
	if [ $nr_cpus != "" ]; then
		MAX_NR_JOBS=$nr_cpus
	fi
fi

OVERALL_COUNT=0
JOB_NR=0
JOB_PIDS=""
for datafile in `ls *.traj.bz2 *.traj`; do
	bztest=${datafile: -3}
	if [ "${bztest}" = "bz2" ]; then
		data=$(head -c 400000 "${datafile}" | bunzip2 -c 2>/dev/null)
		datafile=${datafile:0:${#datafile}-9}
	else
		data=$(head -c 10000 "${datafile}")
		datafile=${datafile:0:${#datafile}-5}
	fi
	conf=$(echo "${data}" | grep configuration | sed "s/configuration = *//g")
	steps=$(MCparams "${conf}" "Simulation Parameters" steps)
	last_step=$(echo "${data}" | grep last_step | sed "s/last_step = *//g")
	if [ ! $last_step -lt $steps ]; then # only analyze if last step in trajectory is last step of simulation (aka we're looking at the last trajectory file and the simulation finished)
		if [ ! -f ${datafile}.stat ]; then
			if [ "${bztest}" = "bz2" ]; then
				echo "-> Uncompressing ${datafile}"
				pbunzip2 ${datafile}.traj.bz2
			fi
			if [ $JOB_NR -lt $MAX_NR_JOBS ]; then
				(traj2stat ${datafile}.traj > ${datafile}.traj2stat && traj2dat ${datafile}.traj > ${datafile}.traj2dat && bzip2 ${datafile}.traj) &
				JOB_PIDS=$(ps -o pgrp,pid -C traj2stat,traj2dat | grep $$ | sed s/$$//g | sed s/' '//g)
				OVERALL_COUNT=$(($OVERALL_COUNT+1))
				JOB_NR=$(echo -n "${JOB_PIDS}" | wc -l)
				echo "Running analysis (and compression) on ${datafile} ($OVERALL_COUNT - $JOB_NR/$MAX_NR_JOBS)"
			fi
			while [ $JOB_NR -ge $MAX_NR_JOBS ]; do
				JOB_PIDS=$(ps -o pgrp,pid -C traj2stat,traj2dat | grep $$ | sed s/$$//g | sed s/' '//g)
				JOB_NR=$(echo -n "${JOB_PIDS}" | wc -l)
				sleep 1
			done
		else
			echo "Already analyzed ${datafile} ..."
			if [ "${bztest}" != "bz2" ]; then
				echo "-> Compressing ${datafile}"
				pbzip2 ${datafile}.traj
			fi
		fi
	fi
done
echo "Waiting for calculations to finish ..."
while [ $JOB_NR -gt 0 ]; do
	JOB_PIDS=$(ps -o pgrp,pid -C traj2stat,traj2dat | grep $$ | sed s/$$//g | sed s/' '//g)
	JOB_NR=$(echo -n "${JOB_PIDS}" | wc -l)
	sleep 1
done
