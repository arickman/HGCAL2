#!/bin/sh

#grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh
#source eosenv.sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/
#INPATH=root://eoscms//eos/cms/store/user/msun/V12
GITVERSION=gitV00-02-09

rm check.log

for alpha in 0.361 #0.297 0.244 0.200 0.164 0.134 0.110
do
    for v in 12
    do
	for p in gamma
	do
	    #for r in 0
	    #do
	    for et in 20 30 40 50 60 70 80 90 100 125 150 175 200
	    do
		echo " - Processing point alpha=$alpha, et=$et : " | tee -a check.log
		eos ls /store/user/amagnan/HGCalEEGeant4/$GITVERSION/$p/HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root > /dev/null
		if (( "$?" != 0 )); then
		    echo " --- Sim File does not exist!" | tee -a check.log
		else 
		    eos ls /store/user/amagnan/HGCalEEGeant4/$GITVERSION/$p/Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root > /dev/null
		    if (( "$?" != 0 )); then
			echo " --- Rec File does not exist!" | tee -a check.log
		    else
			eos ls /store/user/amagnan/HGCalEEGeant4/$GITVERSION/$p/PuMix140_version${v}_model2_BOFF_et${et}_alpha${alpha}.root > /dev/null
			if (( "$?" != 0 )); then
			    echo " --- PU140 File does not exist!" | tee -a check.log
			else
			    ./bin/checkProduction $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root PuMix140_version${v}_model2_BOFF_et${et}_alpha${alpha}.root | tee -a check.log
			fi
		    fi
		fi
	    done
	done
    done
done
