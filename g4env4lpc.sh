#!/bin/bash
export ARCH=x86_64-slc6-gcc46-opt
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.6.3/x86_64-slc6/setup.sh
export QTHOME=/cvmfs/sft.cern.ch/lcg/external/qt/4.8.4/${ARCH}/
export G4BASE=/cvmfs/geant4.cern.ch/geant4
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/${ARCH}/
export FASTJET_INSTALL=/cvmfs/sft.cern.ch/lcg/external/fastjet/3.0.3/${ARCH}/
source $G4BASE/9.6.p02/x86_64-slc6-gcc46-dbg/share/Geant4-9.6.2/geant4make/geant4make.sh
#for boost latest version
export BOOSTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms
export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}:$XERCESCROOT/lib:$HEPMC_DIR/lib:$FASTJET_INSTALL/lib:$BOOSTSYS/lib
source /cvmfs/sft.cern.ch/lcg/app/releases//ROOT/5.34.19/${ARCH}/root/bin/thisroot.sh
#set path = ($DAWNHOME/bin $path $FASTJET_INSTALL/bin)
