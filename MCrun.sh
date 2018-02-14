#!/bin/bash
####################################################
# This file is distributed under the               #
# University of Illinois/NCSA Open Source License. #
# See LICENSE file in top directory for details.   #
#                                                  #
# Copyright (c) 2016 FIGSiM developers             #
####################################################

######################################################################
# Run queue script for Robinson group Monte-Carlo simulation package #
######################################################################
#                                                                    #
# - please do not call directly, is called by MCsubAT                #
#                                                                    #
######################################################################
#                                                                    #
# - first version May 2012, Andreas Tillack                          #
# - complete rewrite July 2012 to use as scheduler across nodes      #
#   - when /tmp/MCqueue is shared (i.e. NFS share) information is    #
#     shared across all nodes (self-announced automatically)         #
#   - code will be expanded to run jobs on any node involved         #
#                                                                    #
######################################################################

MCfig="MCfig"
MAX_NR_JOBS=8

DIR=$PWD
TMP=/tmp/MCqueue
if [ ! -d ${TMP} ]; then
	mkdir ${TMP}
	chmod 777 ${TMP} 2>/dev/null
fi
IN_FSTAB=$(cat /etc/fstab | grep "$TMP")
MOUNTED=$(cat /etc/mtab | grep "$TMP")
if [ "$MOUNTED" == "" -a ! "$IN_FSTAB" == "" ]; then
	mount $TMP
fi
GLOBALCOUNT=${TMP}/MCcount
USER=$(whoami)
NODE=$(hostname -s)
got_count_lock=0

TMPLOCK=${TMP}/MClock.${NODE}
QUEUE=${TMP}/MCqueue.${NODE}
MCPING=${TMP}/MCping

GLOBALCOUNTLOCK_TMP=$GLOBALCOUNT.lock.$NODE.$$
TMPLOCK_TMP=${TMP}/MClock.${NODE}.$$

function GetCountLock {
# get lock on global counter file
	if [ $got_count_lock -eq 0 ]; then
		if [ ! -e "$GLOBALCOUNTLOCK_TMP" ]; then
			touch $GLOBALCOUNTLOCK_TMP 2>/dev/null
			chmod 664 $GLOBALCOUNTLOCK_TMP 2>/dev/null
		fi
		echo -e "$NODE\t$$" > $GLOBALCOUNTLOCK_TMP 2>/dev/null
		tNODE=""
		tPID=""
		while [ "$tPID" != "$$" ] && [ "$tNODE" != "$NODE" ]; do
			if ln "$GLOBALCOUNTLOCK_TMP" "$GLOBALCOUNT.lock" 2>&-; then # we got the lock
				IFS="	" read tNODE tPID 2>/dev/null < $GLOBALCOUNT.lock
			else # we don't
				IFS="	" read tNODE tPID 2>/dev/null < $GLOBALCOUNT.lock
				sleep 0.5 # wait a bit and try again later
			fi
		done
		rm -f $GLOBALCOUNTLOCK_TMP 2>/dev/null
		got_count_lock=1
		echo -e "[$(date)] $NODE got the lock ($$)" >> "$QUEUE.log" 2>/dev/null
	fi
}

function ReleaseCountLock {
	if [ $got_count_lock -eq 1 ]; then
		rm -f $GLOBALCOUNT.lock
		got_count_lock=0
		echo -e "[$(date)] $NODE released lock ($$)" >> "$QUEUE.log" 2>/dev/null
	fi
}

if [ "$2" != "" ]; then
	runfile="$1"
	IFS="	" read JOBNR JOBUSER JOBNODE JOBDIR CONF run_nr 2>/dev/null < $runfile
	
	cd ${TMP}/Q${JOBNR}
	
	outname=$(MCparams $CONF.conf "Simulation Parameters" "FileOut")
	outname=$outname"_"$run_nr
	configname=${CONF##*/}
	configname=$configname"_"$run_nr
	steps=$(MCparams $CONF.conf "Simulation Parameters" steps)
	randsteps=$(MCparams $CONF.conf "Simulation Parameters" randsteps)
	touch ${TMP}/Q${JOBNR}/$outname.traj
	echo steps=${steps} > ${TMP}/Q${JOBNR}/$configname.stats
	echo randsteps=${randsteps} >> ${TMP}/Q${JOBNR}/$configname.stats
	
	MCpid=0
	$MCfig ${CONF}.conf $run_nr 2>&1 | tee $outname.output | gawk "BEGIN{print \"start_time=\"systime(); fflush();}(\$1==${randsteps}){print \"randsteps_time=\"systime(); fflush();}(/^[1-9]/){print \$1\"=\"systime(); fflush();}" >> ${TMP}/Q${JOBNR}/$configname.stats &
	sleep 1
	MCpid=$(ps -o pid,ppid,comm | grep "MCfig" | gawk "(\$2==$$){print \$1}")
	echo -e "\t${MCpid}" >> ${runfile}
	mv ${runfile} ${TMP}/MCrun.$JOBNR
# wait until job is finished
	wait ${MCpid}
# let managing script on original node know that we are finished
	mv $TMP/MCrun.$JOBNR $TMP/_MCfinished.$JOBNR
	exit
fi

# find out how many processors exist on this node (if available)
if [ -e "/proc/cpuinfo" ]; then
	nr_cpus=$(cat /proc/cpuinfo | grep "processor" | tail -n 1 | awk "(\$1==\"processor\"){print \$3+1}")
	if [ $nr_cpus != "" ]; then
		MAX_NR_JOBS=$nr_cpus
	fi
fi
# find out how many scheduled jobs are running
NODE_RUNNING=0
for datafile in `ls -rt $TMP/MCrun.* 2>/dev/null`; do
	IFS="	" read NR USERNAME PROCESS_NODE HOMEDIR CONFIGFILE RUN_NR PID < "${datafile}"
	if [ "$PROCESS_NODE" == "$NODE" ]; then
		let NODE_RUNNING=NODE_RUNNING+1
	fi
done

# find out if we are supposed to be the manager
MANAGER_PID=0
if [ -e "${QUEUE}" ]; then
	IFS="	" read MANAGER_PID NR_RUNNING 2>/dev/null < "${QUEUE}"
	if [ "$MANAGER_PID" == "" ]; then
		echo -e "$$\t$NODE_RUNNING" > "${QUEUE}"
		touch "${QUEUE}" # to update changes over the net
		MANAGER_PID=$$
		echo -e "[$(date)] $NODE manager $$ starting ... (no manager pid in queue file)." >> $QUEUE.log 2>/dev/null
	fi
	RUNNING=`ps ${MANAGER_PID} | grep MCrun.sh | wc -l`
	if [ $RUNNING -lt 1 ]; then
		echo -e "$$\t$NODE_RUNNING" > "${QUEUE}"
		touch "${QUEUE}" # to update changes over the net
		MANAGER_PID=$$
		echo -e "[$(date)] $NODE manager $$ starting ... (no running manager)." >> $QUEUE.log 2>/dev/null
	fi
else
	echo -e "$$\t$NODE_RUNNING" > "${QUEUE}"
	chmod 664 ${QUEUE} 2>/dev/null
	touch "${QUEUE}" 2>/dev/null
	MANAGER_PID=$$
	echo -e "[$(date)] $NODE manager $$ starting ..." >> $QUEUE.log 2>/dev/null
	chmod 664 ${QUEUE}.log 2>/dev/null
fi

if [ "$MANAGER_PID" == "$$" ]; then # this script now acts as the queue manager of this node ...
# make sure there are no other managers around
	echo -e "[$(date)] Making absolutely sure no other managers are running (takes upto 8 secs) ..." >> $QUEUE.log 2>/dev/null
	if [ -e "${MCPING}" ]; then
		rm -f $MCPING
		sleep 4
	fi
	touch $MCPING
	chmod 664 ${MCPING} 2>/dev/null
	sleep 4
	MANAGER_NR=1
	while read ONODE OPID; do
		let MANAGER_NR=MANAGER_NR+1
		if [ "$ONODE" == "$NODE" ]; then
			echo -e "[$(date)] Stopping manager $OPID ..." >> $QUEUE.log 2>/dev/null
			kill $OPID 2>/dev/null
		fi
	done < $MCPING
	rm -f $MCPING
# clean up things which may or may not be there from previous managers
	rm -f $TMPLOCK* 2>/dev/null
	rm -f $GLOBALCOUNT.lock.$NODE* 2>/dev/null
	if [ $MANAGER_NR -eq 1 ]; then # there no other managers out there and it is safe to delete global counter lock file should it exist
		rm -f $GLOBALCOUNT.lock 2>/dev/null
	fi
# make sure we'll exit safely
	trap "rm -f $TMPLOCK 2>/dev/null; rm -f $GLOBALCOUNT.lock 2>/dev/null; rm -f $GLOBALCOUNTLOCK_TMP 2>/dev/null; echo -e \"[$(date)] $NODE manager $$ exiting.\" >> $QUEUE.log 2>/dev/null; exit" 1 2 3 9 13 15
# Since we are manager now we still need to schedule the job ;-)
	if [ "$1" != "" ]; then
		nohup bash MCrun.sh $1 > /dev/null 2>/dev/null < /dev/null &
	fi
	cd $TMP
	
	echo -e "[$(date)] $NODE manager $$ started." >> $QUEUE.log 2>/dev/null
	
	answered_ping=0
	while [ "$MANAGER_PID" == "$$" ]; do
# reply to existence request (Am I, and, if yes, how many?)
		if [ -e "${MCPING}" ]; then
			if [ $answered_ping -eq 0 ]; then
				echo -e "$NODE\t$$" >> $MCPING
				echo -e "[$(date)] Answered ping request ($$)." >> $QUEUE.log 2>/dev/null
				answered_ping=1
			fi
		else
			answered_ping=0
		fi
# take care of finished jobs
		finished=`ls MCfinished.* 2>/dev/null`
		for jobdone in ${finished}; do
			if [ -e "$jobdone" ]; then # safety net in case other node took care of things in the meantime
				IFS="	" read JOBNR JOBUSER JOBNODE JOBDIR JOBCONF NR_SIMILAR JOBPID 2>/dev/null < $jobdone
				if [ "$JOBNODE" == "$NODE" ]; then # only take care of this node's jobs
					GetCountLock
					IFS="	" read MYPID NODE_RUNNING 2>/dev/null < $QUEUE
					let NODE_RUNNING=NODE_RUNNING-1
					echo -e "$MANAGER_PID\t$NODE_RUNNING" > $QUEUE
					touch "${QUEUE}" 2>/dev/null
					chmod 664 ${QUEUE} 2>/dev/null
					IFS="	" read COUNTER NR_RUNNING 2>/dev/null < $GLOBALCOUNT
					let NR_RUNNING=NR_RUNNING-1
					echo -e "$COUNTER\t$NR_RUNNING" > $GLOBALCOUNT
					chmod 664 $GLOBALCOUNT 2>/dev/null
					rm -f $jobdone
					echo -e "[$(date)] -> Job $JOBNR on $NODE finished." >> $QUEUE.log 2>/dev/null
				fi
			fi
		done
# assign number to new jobs
		newcomers=`ls -rt _MCwait.* 2>/dev/null` # oldest first
		for waitfile in ${newcomers}; do
			GetCountLock
			if [ -e "$waitfile" ]; then # safety net in case other node took care of things in the meantime
				echo -e "[$(date)] -> Assigning queue number to $waitfile ($NODE)" >> "$QUEUE.log" 2>/dev/null
				COUNTER=1
				NR_RUNNING=0
				if [ -e $GLOBALCOUNT ]; then
					IFS="	" read COUNTER NR_RUNNING 2>/dev/null < $GLOBALCOUNT
					let COUNTER=COUNTER+1
				fi
				echo -e "$COUNTER\t$NR_RUNNING" > $GLOBALCOUNT
				echo -en "\t${COUNTER}" >> $waitfile
				mv "${waitfile}" "Q${COUNTER}${waitfile}"
				touch "Q${COUNTER}${waitfile}" 2>/dev/null
				chmod 664 $GLOBALCOUNT 2>/dev/null
			fi
		done
# decide if any can run
		waiting=`ls -rt Q*_MCwait.* 2>/dev/null` # oldest first
		for waitfile in ${waiting}; do # do one at a time only
			IFS="	" read MYPID NODE_RUNNING 2>/dev/null < $QUEUE
			if [ $NODE_RUNNING -lt $MAX_NR_JOBS ]; then
				GetCountLock
				if [ -e "$waitfile" ]; then # in case other script took care of it
					IFS="	" read JOBUSER JOBNODE JOBDIR JOBCONF JOBNR RUN 2>/dev/null < $waitfile
					if [ "$RUN" == "" ]; then # mark jobs to start which have not already been scheduled for that
						let NODE_RUNNING=NODE_RUNNING+1
						echo -e "$MYPID\t$NODE_RUNNING" > $QUEUE
						touch "${QUEUE}" 2>/dev/null
						chmod 664 ${QUEUE} 2>/dev/null
						IFS="	" read COUNTER NR_RUNNING 2>/dev/null < $GLOBALCOUNT
						let NR_RUNNING=NR_RUNNING+1
						echo -e "$COUNTER\t$NR_RUNNING" > $GLOBALCOUNT
						chmod 664 $GLOBALCOUNT 2>/dev/null
						if [ "$JOBNODE" == "$NODE" ]; then
							echo -e "\t1" >> $waitfile
						else
							echo -e "$JOBUSER\t$NODE\t$JOBDIR\t$JOBCONF\t$JOBNR\t1" > $waitfile
						fi
						echo -e "[$(date)] -> Job $JOBNR on $NODE marked to start." >> $QUEUE.log 2>/dev/null
					fi
				fi
			fi
		done
		ReleaseCountLock
		sleep 2
	done
else
	configfile=$1
	if [ "$1" == "" ]; then # this should not happen but let's be sure
		exit 1
	fi
	configname=${configfile##*/}
	waitfile=_MCwait.$$
	echo -en "$USER\t$NODE\t$DIR\t$configfile" > $TMP/$waitfile
	chmod 664 $TMP/$waitfile 2>/dev/null
# wait to get a queue number assigned
	while [ -e "$TMP/$waitfile" ]; do
		sleep 0.2
	done
	sleep 0.2
	newwait=$(ls $TMP/Q*${waitfile} 2>/dev/null)
	if [ -e "${newwait}" ]; then
# create run environment
		IFS="	" read JOBUSER JOBNODE JOBDIR JOBCONF JOBNR RUN 2>/dev/null < ${newwait}
		mkdir "${TMP}/Q${JOBNR}" 2>/dev/null
		MCparams "$configfile.conf" "Simulation Parameters" "files_used" > "${TMP}/Q${JOBNR}.files"
		while IFS="	" read run_dir run_file 2>/dev/null; do
			mkdir -p "${TMP}/Q${JOBNR}/${run_dir}" 2>/dev/null
			cp -f "${run_dir}${run_file}" "${TMP}/Q${JOBNR}/${run_dir}${run_file}"
		done < "${TMP}/Q${JOBNR}.files"
		rm -f "${TMP}/Q${JOBNR}.files"
# wait until scheduled to run
		while [ "$RUN" == "" ]; do
			IFS="	" read JOBUSER JOBNODE JOBDIR JOBCONF JOBNR RUN 2>/dev/null < $newwait
			if [ "$RUN" == "" ]; then # no need to wait if we're scheduled to run
				sleep 2
			fi
		done
# find similar (same output filename) jobs
		touch $TMPLOCK_TMP 2>/dev/null
		chmod 664 $TMPLOCK_TMP 2>/dev/null
		echo $$ > $TMPLOCK_TMP
		tPID=""
		while [ "$tPID" != "$$" ]; do
			if ln "$TMPLOCK_TMP" "$TMPLOCK" 2>&-; then # got the lock
				IFS="	" read tPID 2>/dev/null < $TMPLOCK
			else
				IFS="	" read tPID 2>/dev/null < $TMPLOCK
				sleep 0.1 # wait a bit and try again
			fi
		done
		rm -f $TMPLOCK_TMP
		nr_similar=0
		outname=$(MCparams $configfile.conf "Simulation Parameters" "FileOut")
		for similars in `ls $outname\_*.traj $outname\_*.traj.bz2 2>/dev/null`; do
			similar=$(echo "$similars" | sed -n "s/$outname\_[0-9]*./&/p")
			nr_ending=$(echo "$similars" | sed -n "s/$outname\_[0-9]*.//p")
			tnr=${similar:${#outname}+1:${#similar}-${#outname}-${#nr_ending}-2}
			if [ "$tnr" != "" ]; then
				if [ "$tnr" -ge "$nr_similar" ]; then
					nr_similar=$(($tnr+1))
				fi
			fi
		done
		outname=$outname"_"$nr_similar
# start job
		runfile=${TMP}/_MCrun.$JOBNR
		echo -en "$JOBNR\t$USER\t$JOBNODE\t$DIR\t$configfile\t$nr_similar" > $runfile
		chmod 664 $runfile 2>/dev/null
		
		run_nr=" "$nr_similar
		configname=$configname"_"$nr_similar
		
		steps=$(MCparams $configfile.conf "Simulation Parameters" steps)
		randsteps=$(MCparams $configfile.conf "Simulation Parameters" randsteps)
		
		if [ "${NODE}" == "${JOBNODE}" ]; then
			touch $outname.traj
			echo steps=${steps} > ${TMP}/Q${JOBNR}/$configname.stats
			echo randsteps=${randsteps} >> ${TMP}/Q${JOBNR}/$configname.stats
			
			MCpid=0
			$MCfig ${configfile}.conf$run_nr 2>&1 | tee $outname.output | gawk "BEGIN{print \"start_time=\"systime(); fflush();}(\$1==${randsteps}){print \"randsteps_time=\"systime(); fflush();}(/^[1-9]/){print \$1\"=\"systime(); fflush();}" >> ${TMP}/Q${JOBNR}/$configname.stats &
			sleep 1
			MCpid=$(ps -o pid,ppid,command | grep "$MCfig ${configfile}.conf$run_nr" | gawk "(\$2==$$){print \$1}")
			echo -e "\t${MCpid}" >> ${runfile}
		else
			ssh ${JOBNODE} "nohup bash MCrun.sh ${runfile} 1 > /dev/null 2>/dev/null < /dev/null &"
			ln -sf ${TMP}/Q${JOBNR}/$outname.output $outname.output
			ln -sf ${TMP}/Q${JOBNR}/$outname.traj $outname.traj
		fi
		
		rm -f $TMPLOCK # no need to lock anymore
		
		rm -f $newwait # keep queing information around for as long as possible
		
		if [ "${NODE}" == "${JOBNODE}" ]; then
			mv ${runfile} ${TMP}/MCrun.$JOBNR
# wait until job is finished
			wait ${MCpid}
# to let manager know that we are ...
			mv $TMP/MCrun.$JOBNR $TMP/MCfinished.$JOBNR
		else
# wait until job has been finished on other node
			while [ ! -e "$TMP/_MCfinished.$JOBNR" ]; do
				sleep 8
			done
			mv -f ${TMP}/Q${JOBNR}/$outname.output $outname.output
			mv -f ${TMP}/Q${JOBNR}/$outname.traj $outname.traj
			mv $TMP/_MCfinished.$JOBNR $TMP/MCfinished.$JOBNR
		fi
		rm -rf ${TMP}/Q${JOBNR}
	fi
fi

