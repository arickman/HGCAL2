#!/bin/sh

INPATH=root://eoscms//eos/cms/store/user/amagnan/HGCalEEGeant4/
#INPATH=root://eoscms//eos/cms/store/user/msun/V12
GITVERSION=gitV00-02-09

etaset=(17 19 21 23 25 27 29)

for doPU in 140
do
counter=0
    for alpha in 0.361 #0.297 0.244 0.200 0.164 0.134 0.110
    do
	eta=${etaset[$counter]}
	echo "counter:$counter, eta=$eta"
	let counter=$counter+1
	for v in 12
	do
	    for p in gamma
	    do
		t=2
		MYDIR=PLOTS/$GITVERSION/version$v/$p/${t}00um/
		mkdir -p $MYDIR
	    #for r in 0
	    #do
		for et in 30 50 70 90 100 150 200
		do
		    echo "Processing pu${doPU}_v${v}_p${p}_et${et}" #_run${r}
		#if (( "$et"!="80" )); then
		    mkdir -p $MYDIR/eta${eta}_et${et}_pu${doPU}
		    if (( "$doPU"=="0" ))
		    then
			./bin/egammaResolution 0 $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root $MYDIR/eta${eta}_et${et}_pu0.root ${t} | tee egreso_eta${eta}_et${et}_pu0.log
		    else 
			./bin/egammaResolution 0 $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root PuMix${doPU}_version${v}_model2_BOFF_et${et}_alpha${alpha}.root $MYDIR/eta${eta}_et${et}_pu${doPU}.root ${t} | tee egreso_eta${eta}_et${et}_pu${doPU}.log
		    fi
		#else
		#    ./bin/egammaResolution 0 $INPATH/$GITVERSION/$p/ HGcal_version${v}_model2_BOFF_et${et}_alpha${alpha}.root Digi_version${v}_model2_BOFF_et${et}_alpha${alpha}.root $MYDIR/eta${eta}_et${et}.root ${t} >& egreso_eta${eta}_et${et}.log
		#fi
		done
	    done
	done
    done
done