# PFCalEE

Geant4 simulation of a Si-base sampling calorimeter

Check https://twiki.cern.ch/twiki/bin/view/CMS/HGCalPerformanceStudiesWithG4

Geometry implementation is instantiated by detector versions in an enum - cf. src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits - cf. src/EventAction.cc 

When changing the ttree content, adding homemade classes rather than
simple objects, the classes should be added to userlib. The dictionary
for root to understand the classes also need to be remade. Use "make
dictionary" before make, inside of userlib/. Follow instructions in
https://twiki.cern.ch/twiki/bin/view/Sandbox/AnnemarieMagnanSandbox
for what needs to be put in the homemade classes and makefile for root to
understand them (see also example in class userlib/include/HGCSSSimHit.hh).

## Setup the environment (SLC6)

./compile.sh

## Submit in parallel the runs submitProd.py (for lxplus)
## Submit in parallel the runs submitProdLPC.py (for fermilab)

## For loop is to generate several samples with same stat in parallel.
for i in `seq 0 5`; do python submitProd.py  -r ${i} -o /afs/cern.ch/work/o/ocolegro/test/ -n 1000 ; done

## use option -S to test the script locally
## to write to eos, one might replace -o /afs/cern.ch/work/o/ocolegro/test/ with -e /store/user/ocolegro/test.
