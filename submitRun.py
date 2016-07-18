#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-f', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=0,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=0,      type=int)
parser.add_option('-a', '--alpha'       ,    dest='alpha'              , help='incidence angle in rad'       , default=0,      type=float)
parser.add_option('-g', '--gun'         ,    dest='gun'                , help='particle to shoot'            , default='e-')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()


#for en in [5,10,25,40,50,60,80,100,150,200,300,400,500,1000]:
#for en in [188,307,503,829]:
#for en in [5,10,20,25,50,75,100,125,150,175,200,300,500]: 
for en in [10,15,18,20,25,30,35,40,45,50,60,80]:
#for en in [100]: 

    nevents=opt.nevts
    if en>150: nevents=nevents/2

    myqueue=opt.lqueue
    if en<100: myqueue=opt.squeue

    outDir='%s/version_%d/%s/e_%d'%(opt.out,opt.version,opt.gun,en)
    eosDir='%s/%s'%(opt.eos,opt.gun)
    if opt.alpha>0 : outDir='%s_%3.3f'%(outDir,opt.alpha) 
    
    os.system('mkdir -p %s'%outDir)

    #wrapper
    scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
    #scriptFile.write('cd %s\n'%(outDir))
    scriptFile.write('cp %s/g4steer.mac .\n'%(outDir))
    scriptFile.write('PFCalEE g4steer.mac %d %d\n'%(opt.version,opt.model))
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('echo "--Local directory is " $localdir > g4.log\n')
    scriptFile.write('ls * >> g4.log\n')
    if len(opt.eos)>0:
        outTag='version%d_e%d'%(opt.version,en)
        if opt.alpha>0 : outTag='%s_alpha%3.3f'%(outTag,opt.alpha) 
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
        scriptFile.write('cmsStage -f PFcal.root %s/HGcal_%s.root\n'%(eosDir,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log\n')
        scriptFile.write('else\n')
        scriptFile.write('eossize=`$myeos ls -l %s/HGcal_%s.root | awk \'{print $5}\'`\n'%(eosDir,outTag))
        scriptFile.write('localsize=`ls -l PFcal.root | awk \'{print $5}\'`\n')
        scriptFile.write('if (( "$eossize" != "$localsize" )); then\n')
        scriptFile.write('echo " --- Copy to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> g4.log\n')
        scriptFile.write('else\n')
        scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> g4.log\n')
        scriptFile.write('echo " --- File PFcal.root successfully copied to EOS: %s/HGcal_%s.root" >> g4.log\n'%(eosDir,outTag))
        scriptFile.write('rm PFcal.root\n')
        scriptFile.write('fi\n')
        scriptFile.write('fi\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()

    #write geant 4 macro
    g4Macro = open('%s/g4steer.mac'%(outDir), 'w')
    g4Macro.write('/control/verbose 0\n')
    g4Macro.write('/control/saveHistory\n')
    g4Macro.write('/run/verbose 0\n')
    g4Macro.write('/event/verbose 0\n')
    g4Macro.write('/tracking/verbose 0\n')
    #g4Macro.write('/N03/det/setField 3.8 T\n')
    g4Macro.write('/N03/det/setModel %d\n'%opt.model)
    g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
    g4Macro.write('/generator/select particleGun\n')
    g4Macro.write('/gun/particle %s\n'%(opt.gun))    
    g4Macro.write('/gun/energy %f GeV\n'%(en))
    g4Macro.write('/gun/direction %f %f %f\n'%(0.,math.sin(opt.alpha),math.cos(opt.alpha)))
    g4Macro.write('/run/beamOn %d\n'%(nevents))
    g4Macro.close()

    #submit
    os.system('chmod u+rwx %s/runJob.sh'%outDir)
    if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(myqueue,outDir)) 
    else: os.system("bsub -q %s \'%s/runJob.sh\'"%(myqueue,outDir))

