#!/bin/bash

##########  Check this   ############

	### If adding new processes, check submit, some hard coded stuff at beginning  ###
#proc="HNLMonoPhoton.13TeV.UFO"
proc="CrossSection.HNLMonoPhoton.13TeV.UFO"

#proc="HNLPhotonLepton.8TeV.UFO"
#proc="CrossSection.HNLPhotonLepton.8TeV.UFO"

	#Choose number of events
nevents=30000

	#Choose lastFile
firstFile=1
lastFile=1

	#Choose masses
#masses=(0.03 0.3 3 10 30)
masses=(4000)

ncores=2
# Currently, Anal and Write both don't work
pyth="Mad" #Mad Anal Write

	#Choose code
code="HNL"



############## End Check This ##############



for id in `seq $firstFile 1  $lastFile`
do
	for mass in ${masses[*]}
	do
		echo Submitting: ID $id Mass $mass Proc "$proc" Code "$code"
		./submit.Sim.sh $id $mass $pyth $proc $ncores $code $nevents


	done
done
