#!/bin/bash

name=$1 
proc=$2

SCRATCH=Scratchy."$name"

mkdir $SCRATCH

cp HNL.py $SCRATCH
cp -r lib $SCRATCH

cd $SCRATCH



if [[ $proc == *"Mono"* ]]; then 
	python3 HNL.py ../../HNL_UFO_Repo/ HNLMonoPhoton.13TeV.UFO.
fi


if [[ $proc == *"Lepton"* ]]; then 
	python3 HNL.py ../../HNL_UFO_Repo/ HNLPhotonLepton.8TeV.UFO.
fi
