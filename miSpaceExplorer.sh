#!/bin/bash
#
# This script is a "master loop," wherein multiple simulations can be run with a 
# single command. It is intentionally written in a general form, and is meant to be
# used in conjunction with a separate, simulation-specific script. Input parameters
# are set below, and a separate folder is generated for the output of each simulation.



##### START PARAMETERS ##############################################################
#
	runScript="migrate.sh"
	modelsDir="migration"

	tag1="n";	initParam1="1.5"
	tag2="q";	initParam2="1.75"
	tag3="m";	initParam3="2"	
	tag4="M";	initParam4="50.0"  	# star mass (as compared to disk mass)
	tag5="j";	initParam5="53" 	# jin
	tag6="pm";	initParam6="0.001"	# planet mass (percent of star mass)
	tag7="pr";	initParam7="0.1" # 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"	# planet location on disk (normalized from inner edge)
	
	constParam1="1000"		# maxSteps; set to 0 to get results only
	constParam2="512"		# jmax
#
#####################################################################################



chmod +x $runScript
scriptDir=`pwd`

for param1 in $initParam1; do
	for param2 in $initParam2; do
		for param3 in $initParam3; do
			for param4 in $initParam4; do
				for param5 in $initParam5; do
					for param6 in $initParam6; do
						for param7 in $initParam7; do
	
	runTitle=$tag1$param1$tag2$param2$tag3$param3$tag4$param4$tag5$param5$tag6$param6$tag7$param7
	runDir=$scriptDir/$modelsDir/$runTitle
	mkdir -p $runDir
	cd $runDir
	echo "cd $runDir" > $runTitle
	echo "$scriptDir/./$runScript $param1 $param2 $param3 $param4 $param5 $param6 $param7 $constParam1 $constParam2" >> $runTitle

		# execute on compute nodes
	qsub -l nodes=1:ppn=1,walltime=168:00:00 -j oe -o output.txt ./$(basename `pwd`)

	cd $scriptDir
						done
					done
				done
			done
		done
	done
done
