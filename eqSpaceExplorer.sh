#!/bin/bash



##### START PARAMETERS ##############################################################
#
# INITPARAMs can be cycled through a series of nested FOR loops using multiple values
# separated by spaces, with the entire group enclosed in double quotes. TAGs are very
# short descriptions (preferably single letters) describing each INITPARAM.
# CONSTPARAMs will be carried into every simulation.
#
	runScript="eqStarDisks"
	saveDirectory="eqModels"

	tag1="";	initParam1="eq_"
	tag2="n";	initParam2="1.5"
	tag3="q";	initParam3="1.5 1.75 2.0"
	tag4="M";	initParam4="0.000"  #"0.261"  #"0.079"	# star mass start value; increments logarithmically, so zero won't work without modifying eqStarDisks
	tag5="j";	initParam5="30"  #"222"  #"102"		# jin start value
	
	constParam1="512"	# jmax
	constParam2="1"		# number M models to run; values incremented in eqStarDisks
	constParam3="50"		# number jin models to run; values incremented in eqStarDisks
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
	runDir=$scriptDir/$saveDirectory/$runTitle
	mkdir -p $runDir
	cd $runDir
	echo "cd $runDir" > $runTitle
	echo "$scriptDir/./$runScript $param1 $param2 $param3 $param4 $param5 $constParam1 $constParam2 $constParam3" >> $runTitle

		# QSUB RUN COMMAND
	qsub -l nodes=1:ppn=1,walltime=168:00:00 -j oe -o output.txt ./$(basename `pwd`)

	cd $scriptDir
				done
			done
		done
	done
done
