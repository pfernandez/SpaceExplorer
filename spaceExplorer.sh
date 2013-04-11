#!/bin/bash
#
# This script is a "master loop," wherein multiple simulations can be run with a 
# single command. It is intentionally written in a general form, and is meant to be
# used in conjunction with a separate, simulation-specific script. Input parameters
# are set below, and a separate folder is generated for the output of each simulation.
#
#####################################################################################



##### START PARAMETERS ##############################################################
#
# INITPARAMs can be cycled through a series of nested FOR loops using multiple values
# separated by spaces, with the entire group enclosed in double quotes. TAGs are very
# short descriptions (preferably single letters) describing each INITPARAM.
# CONSTPARAMs will be carried into every simulation.
#
	tag1="n";	initParam1="1.5"
	tag2="q";	initParam2="1.5"
	tag3="m";	initParam3="4"	
	tag4="M";	initParam4="0.01"  #"0.0 0.01 0.1 1.0 5.0 10.0 25.0 50.0 100.0"  # star mass
	tag5="j";	initParam5="104"   #"28 53 104 155 206 257 308"  # jin
	
	constParam1="50000"		# maxSteps; set to 0 to get results only
	constParam2="512"		# jmax
#
#####################################################################################



##### MASTER RUN SCRIPT #############################################################
#
# Prepare a separate bash script to begin your simulation using the above parameters
# as input arguments, save it to the current directory, and enter the filename below.
# It may also contain follow-up commands to gather results.
# 
	runScript="starDisks"
	modelsDir="Models"
#					
#####################################################################################



	chmod +x $runScript
	scriptDir=`pwd`
	
	for param1 in $initParam1; do
		for param2 in $initParam2; do
			for param3 in $initParam3; do
				for param4 in $initParam4; do
					for param5 in $initParam5; do
	
	runTitle=$tag1$param1$tag2$param2$tag3$param3$tag4$param4$tag5$param5
	runDir=$scriptDir/$modelsDir/$runTitle
	mkdir -p $runDir
	cd $runDir
	echo "cd $runDir" > $runTitle
	echo "$scriptDir/./$runScript $param1 $param2 $param3 $param4 $param5 $constParam1 $constParam2" >> $runTitle



##### QSUB RUN COMMAND ##############################################################
#
# Starts the simulation on the computational processing nodes. The arguments to this
# command can be modified, and this line can be copied into a conditional statement
# at the end your script to run a simulation recursively. $(basename `pwd`) is the
# name of the working directory, which is also the name of the run command file.
#
	qsub -l nodes=1:ppn=1,walltime=168:00:00 -j oe -o output.txt ./$(basename `pwd`)
#					
#####################################################################################

	
	
	cd $scriptDir
					done
				done
			done
		done
	done
