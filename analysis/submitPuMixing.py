#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

enlist=[0]
if opt.dogun : 
    #enlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200]
    #enlist=[3,5,7,10,20,30,40,50]
    #enlist=[60,70,175,200]
    enlist=[60,70,80,90,100,125,150,175,200]

#for alpha in 0.361 #0.297 0.244 0.200 0.164 0.134 0.110
INPATHPU="root://eoscms//eos/cms/store/user/msun/V12/MinBias/"
nPuVtx=140

for et in enlist :

    nevents=opt.nevts
    
    myqueue=opt.lqueue
    
    bval="BOFF"
    if opt.Bfield>0 : bval="BON" 
    
    outDir='%s/git_%s/version_%d/model_%d/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval)
    if et>0 : outDir='%s/et_%d'%(outDir,et)
    eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
    if opt.alpha>0 : outDir='%s/a_%3.3f/'%(outDir,opt.alpha) 
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    outlog='%s/pumixing%s.log'%(outDir,nPuVtx)
    g4log='pumixjob%s.log'%(nPuVtx)
    os.system('mkdir -p %s'%outDir)
    
    #wrapper
    scriptFile = open('%s/runPUMixJob%s.sh'%(outDir,nPuVtx), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
    #scriptFile.write('cd %s\n'%(outDir))
    outTag='version%d_model%d_%s'%(opt.version,opt.model,bval)
    if et>0 : outTag='%s_et%d'%(outTag,et)
    if opt.alpha>0 : outTag='%s_alpha%3.3f'%(outTag,opt.alpha) 
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('%s/bin/MixPUSignal %s %s \* root://eoscms//eos/cms%s/ Digi_%s.root HGcal_%s.root $localdir/ %s | tee %s\n'%(os.getcwd(),opt.nevts,INPATHPU,eosDir,outTag,outTag,nPuVtx,outlog))
    scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
    scriptFile.write('ls * >> %s\n'%(g4log))
    if len(opt.eos)>0:
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('cmsStage -f PuMix.root %s/PuMix%s_%s.root\n'%(eosDir,nPuVtx,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file PuMix.root to EOS. Keeping locally." >> %s\n'%(g4log))
        scriptFile.write('mv PuMix.root PuMix%s.root\n'%(nPuVtx))
        scriptFile.write('else\n')
        scriptFile.write('eossize=`$myeos ls -l %s/PuMix%s_%s.root | awk \'{print $5}\'`\n'%(eosDir,nPuVtx,outTag))
        scriptFile.write('localsize=`ls -l PuMix.root | awk \'{print $5}\'`\n')
        scriptFile.write('if (( "$eossize" != "$localsize" )); then\n')
        scriptFile.write('echo " --- Copy of pumix file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> %s\n'%(g4log))
        scriptFile.write('mv PuMix.root PuMix%s.root\n'%(nPuVtx))
        scriptFile.write('else\n')
        scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> %s\n'%(g4log))
        scriptFile.write('echo " --- File PuMix.root successfully copied to EOS: %s/PuMix%s_%s.root" >> %s\n'%(eosDir,nPuVtx,outTag,g4log))
        scriptFile.write('rm PuMix.root\n')
        scriptFile.write('fi\n')
        scriptFile.write('fi\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()
    
    #submit
    os.system('chmod u+rwx %s/runPUMixJob%s.sh'%(outDir,nPuVtx))
    if opt.nosubmit : os.system('echo bsub -q %s %s/runPUMixJob%s.sh'%(myqueue,outDir,nPuVtx)) 
    else: os.system("bsub -q %s \'%s/runPUMixJob%s.sh\'"%(myqueue,outDir,nPuVtx))

