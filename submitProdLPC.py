#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
#parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
#parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,     type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=110,    type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=5,      type=int)
parser.add_option('-a', '--eta'         ,    dest='eta'                , help='incidence eta'                , default=0.0,    type=float)
parser.add_option('-p', '--phi'         ,    dest='phi'                , help='incidence phi angle in pi unit' , default=0.0,  type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0.0,    type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-N', '--njobs'       ,    dest='njobs'              , help='number of jobs'           , default=10,   type=int)
parser.add_option('-n', '--nevtsperjob' ,    dest='nevtsperjob'        , help='number of events per job' , default=500,  type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

print "nevents = ",opt.nevtsperjob
print "njobs = ",opt.njobs

enlist=[0]
if opt.dogun : 
    #enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
    #enlist=[2,5,10,20,40,60,80,100,150,200]
    #enlist=[3,5,10,30,50,70,100,200]
    #enlist=[2,5,10,50,100,300,500]
    #enlist=[500]
    #enlist=[8,16,32]
    enlist=[8]

#hgg seeds
#for seed in 1417791355 1417791400 1417791462 1417791488 1417791672 1417791741 1417791747 1417791766 1417791846
#command:
#run=1; for seed in 1417791355 1417791400 1417791462 1417791488 1417791672 1417791741 1417791747 1417791766 1417791846; do ./submitProd.py -s 1nw -q 1nw -t V00-02-14 -v 12 -m 2 -d Hgg -f /afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_${seed}.dat -r ${run} -n 1300 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEGeant4/ -e /store/cmst3/group/hgcal/HGCalEEGeant4; let run=$run+1; done
#vbfHgg seeds
#for seed in 1420833683 1420833689 1420833693 1420833695 1420833696 1420833717

wthick=''
pbthick=''
droplayers=''
label=''

##30
#wthick='1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2'
#pbthick='1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4'
#droplayers=''
#label=''
#label='v5_30'
##28
#wthick='1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2'
#pbthick='1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,5.6,5.6'
#droplayers='25,27'
#label='v5_28'
##24
#wthick='1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2'
#pbthick='2.2,2.2,1,1,2.2,2.1,2.1,3.3,2.1,2.1,4.4,4.4,5.6,5.6'
#droplayers='1,3,10,15,25,27'
#label='v5_24'
##18
#wthick='1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2'
#pbthick='2.2,2.2,2.2,2.2,2.2,2.2,2.1,3.3,3.3,3.3,4.4,5.6,5.6,5.6'
#droplayers='1,3,5,7,10,12,15,18,20,23,25,27'
#label='v5_18'

#granularity='0-20:4,21-63:6'
#noise='0-63:0.12'
#threshold='0-63:2'

#if (opt.version==8) :
#    granularity='0-20:4,21-30:6'
#    noise='0-30:0.14'
#    threshold='0-30:2'
#elif opt.version<20 :
#    granularity='0-19:4,20-29:4'
#    noise='0-29:0.14'
#    threshold='0-29:2'
#elif (opt.version==21 or opt.version==24):
#    granularity='0-23:6,24-33:8'
#    noise='0-33:0.14'
#    threshold='0-33:2'
#elif opt.version==22:
#    granularity='0-9:8'
#    noise='0-9:0.14'
#    threshold='0-9:2'
#elif opt.version==23:
#    granularity='0-53:12'
#    noise='0-53:0.12'
#    threshold='0-53:2'
#elif opt.version>24:
#    granularity='0-29:4,30-53:4,54-65:8'
#    noise='0-53:0.14,54-65:0.2'
#    threshold='0-53:2,54-65:4'

for en in enlist :

    nevents=opt.nevtsperjob
    #if en>150: nevents=nevents/2
    
#    myqueue=opt.lqueue
#    if en>0 and en<60 : myqueue=opt.squeue
    
    bval="BOFF"
    if opt.Bfield>0 : bval="BON" 
    
    outDir='%s/git_%s/version_%d/model_%d/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
    outDir='%s/%s'%(outDir,label) 
    if en>0 : outDir='%s/en_%d'%(outDir,en)
    #if len(opt.eos)>0: eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
    if opt.eta>0 : outDir='%s/eta_%3.3f/'%(outDir,opt.eta)
    if opt.phi!=0.5 : outDir='%s/phi_%3.3fpi/'%(outDir,opt.phi) 
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    os.system('mkdir -p %s'%outDir)

    outTag='%s_version%d_model%d_%s'%(label,opt.version,opt.model,bval)
    if en>0 : outTag='%s_%sen%d'%(outTag,opt.datatype,en)
    if opt.eta>0 : outTag='%s_eta%3.3f'%(outTag,opt.eta) 
    if opt.phi!=0.5 : outTag='%s_phi%3.3fpi'%(outTag,opt.phi) 
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)

    #wrapper
    scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('JOBNUM=$1\n')
    scriptFile.write('shift\n')                         # the shift is necessary to consume the argument prior to sourcing any other scripts!
    scriptFile.write('echo "job number = $JOBNUM"\n')
    scriptFile.write('ls * >> g4_%s_${JOBNUM}.log\n'%outTag)
    scriptFile.write('source ./g4env4lpc.sh\n')
    #scriptFile.write('cd %s\n'%(outDir))
    #scriptFile.write('cp %s/g4steer.mac .\n'%(outDir))
    scriptFile.write('./PFCalEE g4steer_${JOBNUM}.mac %d %d %f %s %s %s >g4_%s_${JOBNUM}.log\n'%(opt.version,opt.model,opt.eta,wthick,pbthick,droplayers,outTag))
    scriptFile.write('mv PFcal.root HGcal_%s_${JOBNUM}.root\n'%(outTag))
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('echo "--Local directory is " $localdir >> g4_%s_${JOBNUM}.log\n'%outTag)
    scriptFile.write('ls * >> g4_%s_${JOBNUM}.log\n'%outTag)

    # if len(opt.eos)>0:
    #     scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
    #     scriptFile.write('source eosenv.sh\n')
    #     scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
    #     scriptFile.write('cmsStage -f HGcal_%s.root %s/HGcal_%s.root\n'%(outTag,eosDir,outTag))
    #     scriptFile.write('if (( "$?" != "0" )); then\n')
    #     scriptFile.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log\n')
    #     scriptFile.write('else\n')
    #     scriptFile.write('eossize=`$myeos ls -l %s/HGcal_%s.root | awk \'{print $5}\'`\n'%(eosDir,outTag))
    #     scriptFile.write('localsize=`ls -l HGcal_%s.root | awk \'{print $5}\'`\n'%(outTag))
    #     scriptFile.write('if (( "$eossize" != "$localsize" )); then\n')
    #     scriptFile.write('echo " --- Copy of sim file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> g4.log\n')
    #     scriptFile.write('else\n')
    #     scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> g4.log\n')
    #     scriptFile.write('echo " --- File PFcal.root successfully copied to EOS: %s/HGcal_%s.root" >> g4.log\n'%(eosDir,outTag))
    #     scriptFile.write('rm HGcal_%s.root\n'%(outTag))
    #     scriptFile.write('fi\n')
    #     scriptFile.write('fi\n')

    scriptFile.write('echo "--deleting core files: too heavy!!"\n')
    scriptFile.write('rm -f core.*\n')
    #scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()
    
    #write geant 4 macro
    for j in xrange(0,opt.njobs):
        g4Macro = open('%s/g4steer_%d.mac'%(outDir,j), 'w')
        g4Macro.write('/control/verbose 0\n')
        g4Macro.write('/control/saveHistory\n')
        g4Macro.write('/run/verbose 0\n')
        g4Macro.write('/event/verbose 0\n')
        g4Macro.write('/tracking/verbose 0\n')
        g4Macro.write('/N03/det/setField %1.1f T\n'%opt.Bfield)
        g4Macro.write('/N03/det/setModel %d\n'%opt.model)
        g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
        if opt.dogun :
            g4Macro.write('/generator/select particleGun\n')
            g4Macro.write('/gun/particle %s\n'%(opt.datatype))
        #if opt.eta<5 : en=et*math.cosh(opt.eta)
        #else : en=et
            g4Macro.write('/gun/energy %f GeV\n'%(en))
            if opt.model==5 :
                g4Macro.write('/gun/direction 0 0 1\n')
            elif opt.model!=2 :
                alpha = 2*math.atan(math.exp(-1.*opt.eta));
                g4Macro.write('/gun/direction %f %f %f\n'%(math.cos(math.pi*opt.phi)*math.sin(alpha),math.sin(math.pi*opt.phi)*math.sin(alpha),math.cos(alpha)))
        #g4Macro.write('/gun/direction %f %f %f\n'%(math.cos(math.pi*opt.phi)*math.sin(opt.alpha),math.sin(math.pi*opt.phi)*math.sin(opt.alpha),math.cos(opt.alpha)))
        #g4Macro.write('/gun/direction %f %f %f\n'%(random.uniform(0,1000)/100.-5.,math.sin(opt.alpha),math.cos(opt.alpha)))
            #else gun direction hardcoded in PrimaryGeneratorAction
        else :
            g4Macro.write('/generator/select hepmcAscii\n')
            g4Macro.write('/generator/hepmcAscii/open %s\n'%(opt.datafile))
            g4Macro.write('/generator/hepmcAscii/verbose 0\n')
        g4Macro.write('/run/beamOn %d\n'%(nevents))
        g4Macro.close()

    # Condor jdl file
    jdlFile = open('%s/submitProd.jdl'%(outDir), 'w')
    jdlFile.write('Executable = %s/runJob.sh\n'%outDir)
    jdlFile.write('Universe = vanilla\n')
    jdlFile.write('Requirements = FileSystemDomain=="fnal.gov" && Arch=="X86_64"\n')
    jdlFile.write('Notification = ERROR\n')
    jdlFile.write('Should_Transfer_Files = YES\n')
    jdlFile.write('WhenToTransferOutput = ON_EXIT\n')
    for j in xrange(0,opt.njobs):
        jdlFile.write('Arguments = %d\n'%j)
        jdlFile.write('transfer_input_files = %s/g4env4lpc.sh,%s/g4steer_%d.mac,%s/bin/Linux-g++/PFCalEE,%s/tmp/Linux-g++/PFCalEE/libPFCalEE.so,%s/analysis/lib/libPFCalEEAnalysis.so,%s/userlib/lib/libPFCalEEuserlib.so\n'%(opt.out,outDir,j,os.environ['G4WORKDIR'],os.environ['G4WORKDIR'],opt.out,opt.out))
        jdlFile.write('Error = pfcalee_%s_%d.stderr\n'%(outTag,j))
        jdlFile.write('Output = pfcalee_%s_%d.stdout\n'%(outTag,j))
        jdlFile.write('Queue\n')
    jdlFile.close()
    
    #submit
    os.system('chmod u+rwx %s/runJob.sh'%outDir)
    if opt.nosubmit : os.system('echo condor_submit %s/submitProd.jdl'%outDir) 
    else: os.system('condor_submit %s/submitProd.jdl'%outDir)
