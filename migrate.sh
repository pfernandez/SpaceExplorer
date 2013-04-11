#!/bin/bash
#
# This script will be executed when miSpaceExplorer.sh is run.
#
# Compiles and runs equilibrium and migration models and obtains results.
# To restart models automatically if Y2 has not converged when maxSteps
# has been reached, set "neverRestart" to "false" in the notDone
# function at the top of the getValues script.


n=$1
q=$2
m=$3
starMass=$4
jin=$5
planetMass=$6
planetRad=$7
maxSteps=$8
jmax=$9
scriptDir=$(dirname `dirname $(pwd)`)
val="$scriptDir/./getValues $(pwd)"
chmod +x $scriptDir/getValues



####################################################################################################
# RUN EQUILIBRIUM MODEL

echo "Compiling hscf.f..."
cp $scriptDir/hscf.f .
echo "
	  parameter (jmax = ${jmax},
     1           kmax = ${jmax}, 
     2           max = 25000)

c     32x32     wfw(452)
c     64x64     wfw(1030)
c     128x128   wfw(2315)
c     256x256   wfw(5132)
c     512x512   wfw(11279)
c     1024x1024 wfw(24593)
c     2048x2048 wfw(53264)
" > param.h

# ifort -O -mp -r8 -o hscf hscf.f
gfortran -fdefault-real-8 -O2 -o hscf hscf.f

echo "Running equilibrium model..."
echo "1                                      :  0 = white dwarf   1 = polytrope
0                                      :  0 = dead start    1 = jump start
${n} -${q} ${jmax} ${jmax} ${jmax} -${jin} 0.0 0.8 ${starMass}  :  n n' jold kold jout kout log(rho) del starm
" > hscf.in
./hscf < hscf.in
rm param.h hscf hscf.in hscf.f
	
	#plot eqContour
echo "Generating contour plot..."
r=`echo "scale=2; $($val rInOut)/1" | bc`
echo -e "
reset\n
set terminal png\n
set output \"eqContour.png\"\n
set title \"n${n}/q${q}/r-+${r}-j${jin}/M${starMass}/${jmax}\"\n
set contour\n
set cntrparam levels incremental 1.0e-30,$(echo "scale=7; $($val rhomax)/10" | bc),$($val rhomax)\n
set size square\n
set grid\n
set mxtics\n
set mytics\n
set nosurface\n
set view 0,0\n
set data style lines\n
set nokey\n
splot 'fort.47' ti \"fort.47\"\n
" | gnuplot

#
####################################################################################################
# RUN MIGRATION MODEL

if [ $maxSteps -gt 0 ]; then
	
echo "Compiling migrate.f..."
cp $scriptDir/migrate.f .
echo "
c      version: 9 May 2005

      parameter (jmax=${jmax},jmax1=jmax+1,jmax2=jmax+2,jmax3=jmax-1,
     1           kmax=${jmax},kmax1=kmax+1,kmax2=kmax+2,kmax3=kmax-1,
     2           lmax=16,lmax1=lmax/2+1,lmax2=lmax-2,
     3           max=25000,jkm1=jmax+kmax-1,neq=8)

c     32x32     wfw(452) 
c     64x64     wfw(1030)
c     128x128   wfw(2312)
c     256x256   wfw(5130)
c     512x512   wfw(11276)
c     1024x1024 wfw(24590)
" > param.h

# ifort -O -mp -r8 -o hscf hscf.f
gfortran -fdefault-real-8 -O2 -o migrate migrate.f
	
echo "Running migration model..."
rerun=0
Time=0
if [ -f fort.11 ] && [ ! -f fort.51 ]; then
	mv fort.11 fort.51
fi
	###### restart if fort.51 already present
if [ -f fort.51 ]; then
	rerun=1
	mv fort.51 fort.11
	if [ -f fort.22 ]; then
		numlines=$(cat fort.22 | wc -l)
		ourline=$(( $numlines - $numlines % 100 ))
		Time=$(sed -n ${ourline}p fort.22 | cut -c 1-15)
	fi
	if [ -f fort.23 ]; then 
		if [ -f fort.23a ]; then
			cp fort.23a temp
			cat fort.23 >> temp
			rm fort.23 fort.23a
			mv temp fort.23
		fi
		cp fort.23 fort.23a
	fi
fi
	
	##### create migrate.in file
echo "fort.2
${m}                          : m
1                          : model number
${n},000,.75               : n, 0 to use data file and 100xq otherwise, rin
${q}                        : q
${planetRad},${planetMass},${starMass},10.0e10  : see below
${maxSteps}                         : max steps  
${rerun},${Time}                        : 1 to restart,time; mv fort.51 fort.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%
Key:

! background (equilibrium) model location
! azimuthal mode wave number
! background model number in file
! polytropic index, Omega exponent x 100, disk edge / radius of density maxium
! semi-major axis of planet's orbit, planet's mass,  stellar mass, growth rate 
! 0=dead start, 1= continuation; time at start of run (the restart time) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2%
" > migrate.in

	##### execute with migrate.in
./migrate < migrate.in
fi
	# remove unnecessary files
rm param.h migrate migrate.in migrate.f
#rm fort.2 fort.15 fort.16 fort.17 fort.30 fort.31 fort.42 fort.52 fort.55 fort.56 fort.57
#rm fort.71 fort.93 anal.dat eta fact1 outfile pulsefile rho.dat wdout



#
####################################################################################################
# GET RESULTS


echo "Getting model results..."

	# Combine new fort.23 if restarting
if [ -f fort.23a ]; then
	cp fort.23a temp
	cat fort.23 >> temp
	rm fort.23 fort.23a
	mv temp fort.23
fi


	##### Generate and save images with gnuplot.
r=`echo "scale=2; $($val rInOut)/1" | bc`
title="m${m}/n${n}/q${q}/r-+${r}-j${jin}/M${starMass}/pr${planetRad}/pm${planetMass}/${jmax}"


	#plotit
plotCommand="
reset\n
set terminal png\n
set output \"plotit.png\"\n
set title \"${title}\"\n
set logscale y\n
set autoscale\n
plot 'fort.23' using 1:2 with lines ti \"fort.23\",\
	 'fort.23' using 1:4 with lines notitle,\
	 'fort.23' using 1:6 with lines notitle\n
"
echo -e $plotCommand | gnuplot


	#plotit2
echo -e "
reset\n
set terminal png\n
set output \"plotit2.png\"\n
set title \"${title}\"\n
set auto\n
plot 'fort.23' using 1:3 with lines ti \"fort.23\",\
	 'fort.23' using 1:5 with lines notitle,\
	 'fort.23' using 1:7 with lines notitle\n
" | gnuplot


	#torque1
echo -e "
reset\n
set terminal png\n
set output \"torque1.png\"\n
set title \"${title}\"\n
set auto\n
plot 'fort.27' using 1:2 with lines ti \"total_torque\"\n
" | gnuplot


	#torque2
echo -e "
reset\n
set terminal png\n
set output \"torque2.png\"\n
set title \"${title}\"\n
set auto\n
plot 'fort.27' using 1:3 with lines ti \"torkp\"\n
" | gnuplot


	#torque3
echo -e "
reset\n
set terminal png\n
set output \"torque3.png\"\n
set title \"${title}\"\n
set auto\n
plot 'fort.27' using 1:4 with lines ti \"torksum\"\n
" | gnuplot


#if [ ]; then #BEGIN COMMENTBLOCK
####################################################################################################
# RESTART IF REQUIRED

if [ $($val notDone) ] && [ $maxSteps -gt 0 ]; then
	echo -e "\nRestarting..."
	qsub -l nodes=1:ppn=1,walltime=168:00:00 -j oe -o output.txt ./$(basename `pwd`)
else
	echo -e "\nDone."
fi

#fi #END COMMENTBLOCK
####################################################################################################