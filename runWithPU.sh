#!/bin/sh

#for eta in 20 25 30 35
#do
#for set in 0 1 2 3 4 5 6 7 8 9
#do
    for type in Higgs #PedroPU
    do
	#for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/$type/eta$eta/*.dat`;
	#for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/version_4/$type/*.dat`;
	for justpileup in `ls /afs/cern.ch/work/p/pdauncey/public/annemarie/version_5/Signal/$type/$type*.dat`;
	do
	    echo "$justpileup" > tmp1
	    #sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/${type}/eta${eta}/pileup0\{0,3\}||" tmp1 > tmp2
	    #sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/version_4/${type}/Photon0\{0,3\}||" tmp1 > tmp2
	    sed "s|/afs/cern.ch/work/p/pdauncey/public/annemarie/version_5/Signal/${type}/${type}Eta||" tmp1 > tmp2
	    #subdir=eta$eta"_e"`sed "s|[A-Z]\.dat||" tmp2`
	    #subdir=set$set"_e"`sed "s|[A-Za-z]\{0,3\}\.dat||" tmp2`
	    subdir="eta"`sed "s|\.dat||" tmp2`

	    echo " -- processing file "$justpileup 
	    echo " -- subdir " $subdir
	
	    ./submitRunHEPMC.py -q 2nd -v 8 -m 2 -f $justpileup -d $type -s $subdir -n 1000 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/ -e /store/user/amagnan/HGCalEEGeant4

	done
    done
#done
#done
