#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
#Long Que is 2nw
parser.add_option('-q', '--queue' ,    dest='queue'             , help='batch queue'            , default='1nd')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='hexaV02-01-01')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=1,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=2,      type=int)
parser.add_option('-s', '--signal'      ,    dest='signal'             , help='signal flag'                  , default=0,      type=int)
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
##Location to a txt file containing a csv [particle_energy,dir_x,dir_y,dir_z]
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')

(opt, args) = parser.parse_args()

label=''


nevents=opt.nevts
myqueue=opt.queue

thickness_ = [1]#,2,3,4,5,6,7,8,9,10,12,14,16,18,20]
for thickness in thickness_:
    outDir='%s/git_%s/version_%d/model_%d'%(opt.out,opt.gittag,opt.version,opt.model)
    outDir='%s/%s'%(outDir,label)
    eosDir='%s/git%s'%(opt.eos,opt.gittag)
    if opt.signal>0 : outDir='%s/signal_%3.3f/'%(outDir,opt.signal)
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)
    
    os.system('mkdir -p %s'%outDir)
    
    #wrapper
    scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
    scriptFile.write('cp %s/g4steer.mac .\n'%(outDir))
    outTag='%s_version%d_model%d_thick%s'%(label,opt.version,opt.model,thickness)
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)

    if len(opt.eos)>0:
        scriptFile.write('PFCalEE g4steer.mac %d %d %f %s root://eoscms.cern.ch//eos/cms%s/HGcal_%s.root| tee g4.log\n'%(opt.version,opt.model,opt.signal,thickness,eosDir,outTag))
    if len(opt.eos) == 0:
        scriptFile.write('PFCalEE g4steer.mac %d %d %f %s %s/HGcal_%s.root| tee g4.log\n'%(opt.version,opt.model,opt.signal,thickness,outDir,outTag))

    if (opt.run>=0) : outTag='%s_second'%(outTag)
    scriptFile.write('mv PFcal.root HGcal_%s.root\n'%(outTag))
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('echo "--Local directory is " $localdir >> g4.log\n')
    scriptFile.write('ls * >> g4.log\n')
    if len(opt.eos)>0:
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
        scriptFile.write('$myeos cp HGcal_%s.root /eos/cms%s/HGcal_%s.root\n'%(outTag,eosDir,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log\n')
        scriptFile.write('fi\n')
        scriptFile.write('rm HGcal_%s.root\n'%(outTag))



    scriptFile.write('echo "--deleting core files: too heavy!!"\n')
    scriptFile.write('rm core.*\n')
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
    g4Macro.write('/N03/det/setModel %d\n'%opt.model)
    g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
    g4Macro.write('/run/beamOn %d\n'%(nevents))
    g4Macro.close()
    
    #submit
    os.system('chmod u+rwx %s/runJob.sh'%outDir)
    if opt.nosubmit : os.system('LSB_JOB_REPORT_MAIL=N echo bsub -q %s -N %s/runJob.sh'%(myqueue,outDir))
    else: os.system("LSB_JOB_REPORT_MAIL=N bsub -q %s -N \'%s/runJob.sh\'"%(myqueue,outDir))

