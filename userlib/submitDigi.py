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
parser.add_option('-q', '--queue' ,    dest='queue'             , help='long batch queue'             , default='1nd')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='hexaV02-01-01')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=1,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=2,      type=int)
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()


nSiLayers=2
label=''
INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/"

nPuVtx=0
interCalibList=[3];

granularity='0-27:4'
noise='0-27:0.15'
threshold='0-27:5'


for interCalib in interCalibList:
    if nPuVtx>0 :
        suffix='Pu%d_IC%d'%(nPuVtx,interCalib)
    else :
        suffix='IC%d'%(interCalib)
    myqueue=opt.queue

    if opt.model!=2 : suffix='%s_Si%d'%(suffix,nSiLayers)



    outDir='%s/git_%s/version_%d/model_%d/'%(opt.out,opt.gittag,opt.version,opt.model)
    outDir='%s/%s'%(outDir,label)


    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    if len(opt.eos)>0:
        eosDir='%s/git%s'%(opt.eos,opt.gittag)
        eosDirIn = 'root://eoscms//eos/cms%s/' % (eosDir)


    else:
        eosDir='%s/'%(outDir)
        eosDirIn = '%s/' % (eosDir)

    #Digi%s_%s.root\n'%(eosDir,suffix,outTag)
    outTag='%s_version%d_model%d'%(label,opt.version,opt.model)
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)    
#    eosDirIn = 'root://eoscms//eos/cms%s/HGcal%s_%s.root' % (eosDir,suffix,outTag)
    print 'The eosdir is ' + str(eosDir)
    print 'The eosdirIn is ' + str(eosDirIn)
    outlog='%s/digitizer%s.log'%(outDir,suffix)
    g4log='digijob%s.log'%(suffix)
    os.system('mkdir -p %s'%outDir)

    #wrapper
    scriptFile = open('%s/runDigiJob%s.sh'%(outDir,suffix), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('%s/bin/digitizer %d %s/HGcal_%s.root $localdir/ %s %s %s %d %d %d %s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,outTag,granularity,noise,threshold,interCalib,nSiLayers,nPuVtx,INPATHPU,outlog))
    scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
    scriptFile.write('ls * >> %s\n'%(g4log))
    if len(opt.eos)>0:
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
        scriptFile.write('$myeos cp DigiPFcal.root /eos/cms%s/Digi%s_%s.root\n'%(eosDir,suffix,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally." >> %s\n'%(g4log))
        scriptFile.write('else\n')
        scriptFile.write('eossize=`$myeos ls -l %s/Digi%s_%s.root | awk \'{print $5}\'`\n'%(eosDir,suffix,outTag))
        scriptFile.write('localsize=`ls -l DigiPFcal.root | awk \'{print $5}\'`\n')
        scriptFile.write('if (( "$eossize" != "$localsize" )); then\n')
        scriptFile.write('echo " --- Copy of digi file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> %s\n'%(g4log))
        scriptFile.write('else\n')
        scriptFile.write('echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> %s\n'%(g4log))
        scriptFile.write('echo " --- File DigiPFcal.root successfully copied to EOS: %s/Digi%s_%s.root" >> %s\n'%(eosDir,suffix,outTag,g4log))
        scriptFile.write('rm DigiPFcal.root\n')
        scriptFile.write('fi\n')
        scriptFile.write('fi\n')
    else:
        scriptFile.write('mv DigiPFcal.root Digi%s_%s.root\n'%(suffix,outTag))

    scriptFile.write('echo "--deleting core files: too heavy!!"\n')
    scriptFile.write('rm core.*\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()

    #submit
    os.system('chmod u+rwx %s/runDigiJob%s.sh'%(outDir,suffix))
    if opt.nosubmit : os.system('LSB_JOB_REPORT_MAIL=N echo bsub -q %s -N %s/runDigiJob%s.sh'%(myqueue,outDir,suffix))
    else: os.system("LSB_JOB_REPORT_MAIL=N bsub -q %s -N \'%s/runDigiJob%s.sh\'"%(myqueue,outDir,suffix))



