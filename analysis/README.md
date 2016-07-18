##################################
## Example analysis

#to compile

mkdir {lib,bin,obj,PLOTS}
make

#any .cc and .hh files will be compiled into obj files and into a library in lib/
#any .cpp file in test/ will be compiled as executable, in bin/

#./bin/studyTransverseProperties -i ../ -v 20 -e 5 <- is this needed
./bin/compareCaloStackPerformances

###################################
EOS file paths
###################################
root://eoscms//eos/cms/
/store/user/amagnan/
/store/cmst3/group/hgcal/Geant4


###################################
Info on existing Executables:
###################################

#MixPUSignal
#to submit all jobs:
for alpha in 0.361 0.297 0.244 0.200 0.164 0.134 0.110; do ./submitPuMixing.py -S -q 1nd -t V00-02-09 -g -v 12 -m 2 -d gamma -a $alpha -n 0 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/ -e /store/user/amagnan/HGCalEEGeant4 ; done

# Example submit script in:
./runAllFill.sh

# Example plotting macros in macros/plotE.C, plotXY.C, etc...
cd macros
root plotE.C++