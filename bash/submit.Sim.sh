#!/bin/bash
ID=$1
mass=$2
pyth=$3
proc=$4
ncores=$5
code=$6
nevents=$7

rseed=$RANDOM
SCRATCH=Scratchy."$ID"."$mass"."$proc"."$rseed"
#nmass=$((${#mass}-1))

bwscal=10.0
bwmasscut=0.01
if (( $(echo "$mass < $bwmasscut" |bc -l) )); then
        bwscal=20.0
fi
bwvalue=$(python -c "print $mass/$bwscal")




########   CHECK THESE  #############


HOME=/Users/chaffy/Work/Repository_HNL/HNL_UFO_Repo
MADGRAPH=/Users/chaffy/Work/MG5_aMC_v2_4_0
BASH=/Users/chaffy/Work/Repository_HNL/bash

madspin="off"

if [[ $proc == *"HNL"* ]]; then
        name="$proc"."$mass"GeV."$ID"
        nameNoID="$proc"."$mass"GeV
        nameWrite="$code"."$proc"."$mass"
fi

#########   END CHECK    ############

if [ $ncores = 1 ]; then
	madMode=0
else
	madMode=2
fi

cd $BASH
mkdir $SCRATCH

############   RUNNING MADGRAPH  ###############
if [[ $pyth == *"Mad"* ]]; then

	cd $BASH
	cp -r $MADGRAPH/$proc $SCRATCH

	cd $SCRATCH/$proc/Cards/

	line=$rseed"   = iseed   ! rnd seed (0=assigned automatically=default))"
	gsed -i "s/.*iseed.*/  $line/" run_card.dat

	line=$nevents" = nevents ! Number of unweighted events requested"
	gsed -i "s/.*nevents ! Number of unweighted.*/  $line/" run_card.dat
	

	line=$bwvalue"  = bwcutoff"
	gsed -i "s/.*= bwcutoff.*/  $line/" run_card.dat

	
	if [[ $proc == *"HNL"* ]]; then
		gsed -i '/^  999999/c\  999999 '"$mass"' # n : mN' param_card.dat
		gsed -i '/Block frblock $/{n;{s/.*/    1 '"$mass"' # mN/}}' param_card.dat
#		gsed -i '/Block frblock $/{n;{s/.*/    1 '"${mass:0:1}"'.'"${mass:1:1}"''"${mass:2:2}"'0000e+0'"$nmass"' # mN/}}' param_card.dat
	fi

	../bin/generate_events $madMode $ncores $proc #used to be 2 8 $name for pythia interface..
	
	cd ../Events/

	if [ $madspin = off ]; then
		gunzip "$proc"/unweighted_events.lhe.gz
		cp "$proc"/unweighted_events.lhe $HOME/"$name"_unweighted_events.lhe

                gunzip 1/unweighted_events.lhe.gz
                cp 1/unweighted_events.lhe $HOME/"$name"_unweighted_events.lhe

		cd $HOME
		gsed -i "/<wgt id=/d" "$name"_unweighted_events.lhe

	elif [ $madspin = on ]; then
                gunzip "$proc"_decayed_1/unweighted_events.lhe.gz
                cp "$proc"_decayed_1/unweighted_events.lhe $HOME/"$name"_unweighted_events.lhe

                gunzip 1_decayed_1/unweighted_events.lhe.gz
                cp 1_decayed_1/unweighted_events.lhe $HOME/"$name"_unweighted_events.lhe
	fi

	if [[ $proc == *"CrossSection"* ]]; then
		cd $HOME
		gsed -i '/#  Integrated weight (pb)  :/!d' "$name"_unweighted_events.lhe
		gsed -i 's/#  Integrated weight (pb)  ://' "$name"_unweighted_events.lhe
		gsed -i 's/ //g' "$name"_unweighted_events.lhe
	fi

fi

################## Analyse Code  #####################

if [[ $pyth == *"Anal"* ]]; then

	cd $BASH
	cp Makefile myEvent.cc myEvent.hh fjcore.cc fjcore.hh "$code".cc $SCRATCH;
	cd $SCRATCH

	if [ $code = LHCO2Lepton ]; then
		g++ -o "$code".o "$code".cc fjcore.cc
		./"$code".o $HOME/"$nameNoID" $ID $ID $HOME/Accepted/accepted."$nameWrite".dat $rseed
	elif [[ $code == *"scalar_fourjet_resonance"* ]]; then
		make "$code"
	        ./"$code".o $HOME/hepmc-cairns."$name".hepmc $HOME/Accepted/accepted."$nameWrite"."$ID".dat
	fi
fi


########################   WRITE   ########################

if [[ $pyth == *"Write"* ]]; then
	cd $HOME/Accepted/

	bad="0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1"

	if [ $code = LHCO2Lepton ]; then
		echo $(grep -r -i "$bad" accepted."$nameWrite".dat | wc -l) >> nEvaccepted."$nameWrite".dat
		sed -i "/$bad/d" accepted."$nameWrite".dat
	elif [[ $code == *"scalar_fourjet_resonance"* ]]; then
		echo $(grep -r -i "$bad" accepted."$nameWrite"."$ID".dat | wc -l) $ID >> nEvaccepted."$nameWrite".dat
		sed -i "/$bad/d" accepted."$nameWrite"."$ID".dat
		cat accepted."$nameWrite"."$ID".dat >> accepted."$nameWrite".dat
		rm accepted."$nameWrite"."$ID".dat
	fi
fi

cd $BASH
rm -r $SCRATCH

if [[ $pyth == *"Pyth"* ]]; then
	rm $HOME/hepmc-cairns."$name".hepmc
fi
