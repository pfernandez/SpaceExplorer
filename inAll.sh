#!/bin/bash
#
# Perform actions in EVERY run folder.
#
# Search Models folder for all models, descend into each
# folder and execute desired script on compute nodes.
# Currently set to run "starDisks" with modified settings.
#
# *Excludes currently running models by default.


maxSteps=0
excludeRunningModels=true
jmax=512 # overridden if polyout is present
runScript="migrate.sh"
scriptDir=$(pwd)
modelsDir=$scriptDir/migration

chmod +x $runScript
for run in $(ls $modelsDir); do
	runDir=$modelsDir/$run
	if [ -d $runDir ] && [ -n "$(echo `ls $runDir` | grep fort.)" ]; then
		if [ "$excludeRunningModels" = "true" ] && [ $($scriptDir/./getValues $runDir isRunning) ]; then
			echo "Omitting active model $run."
		else
			cd $runDir
			if [ "$runScript" = "migrate.sh" ]; then
				pr=${run##*pr}
				pm=${run%%pr*}; pm=${pm##*pm}
				j=${run%%pm*}; j=${j##*j}
			else
				j=${run##*j}
			fi
			M=${run%%j*}; M=${M##*M}
			m=${run%%M*}; m=${m##*m}
			q=${run%%m*}; q=${q##*q}
			n=${run%%q*}; n=${n##*n}
			if [ -f $runDir/polyout ]; then
				jmax=$($scriptDir/./getValues $runDir jmax)
			fi
			
			#if [ $($scriptDir/./getValues $runDir notDone) ]; then
			echo $run
			echo "cd $runDir" > $run
			if [ "$runScript" = "migrate.sh" ]; then
				echo "$scriptDir/./$runScript $n $q $m $M $j $pm $pr $maxSteps $jmax" >> $run
			else
				echo "$scriptDir/./$runScript $n $q $m $M $j $maxSteps $jmax" >> $run
			fi
			qsub -l nodes=1:ppn=1,walltime=168:00:00 -j oe -o output.txt ./$run
			#fi
		fi
	fi
done
