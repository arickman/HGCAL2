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
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-M', '--macro'       ,    dest='macro'              , help='name of script to run',  default='analyze_hadrons')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')

(opt, args) = parser.parse_args()
label       = ''
suffix      = ''
myqueue=opt.queue
thickness_ = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]
for thickness in thickness_:
    outDir='%s/git_%s/version_%d/model_%d/'%(opt.out,opt.gittag,opt.version,opt.model)
    outDir='%s/%s'%(outDir,label)

    print opt.macro
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    if len(opt.eos)>0:
        eosDir='%s/git%s'%(opt.eos,opt.gittag)
        eosDirIn = 'root://eoscms//eos/cms%s/' % (eosDir)

    else:
        eosDir='%s/'%(outDir)
        eosDirIn = '%s/'%(outDir)
    outTag='%s_version%d_model%d_thick%s'%(label,opt.version,opt.model,thickness)
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)

    print 'The eosdir is ' + str(eosDir)
    print 'The eosdirIn is ' + str(eosDirIn)
    outlog='%s/digitizer%s.log'%(outDir,suffix)
    g4log='digijob%s.log'%(suffix)
    os.system('mkdir -p %s'%outDir)

    #wrapper
    scriptFile = open('%s/runGEN%s.sh'%(outDir,suffix), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/../g4env.sh\n'%(os.getcwd()))
    scriptFile.write('localdir=`pwd`\n')
    scriptFile.write('%s/bin/%s  %s/HGcal_%s.root | tee %s\n'%(os.getcwd(),opt.macro,eosDirIn,outTag,outlog))
    scriptFile.write('echo "--Local directory is " $localdir >> %s\n'%(g4log))
    scriptFile.write('ls * >> %s\n'%(g4log))


    if len(opt.eos)>0:
        scriptFile.write('grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh\n')
        scriptFile.write('source eosenv.sh\n')
        scriptFile.write('$myeos mkdir -p %s\n'%eosDir)
        scriptFile.write('$myeos cp analyzed_tuple.root /eos/cms%s/%s%s_%s.root\n'%(eosDir,opt.macro,suffix,outTag))
        scriptFile.write('if (( "$?" != "0" )); then\n')
        scriptFile.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log\n')
        scriptFile.write('fi\n')
        scriptFile.write('rm analyzed_tuple.root')
    else:
        scriptFile.write('mv analyzed_tuple.root %s/%s%s_%s.root' %(eosDir,opt.macro,suffix,outTag))
        scriptFile.write('rm analyzed_tuple.root')


    scriptFile.write('\nrm core.*\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()

    #submit
    os.system('chmod u+rwx %s/runGEN%s.sh'%(outDir,suffix))
    if opt.nosubmit : os.system('LSB_JOB_REPORT_MAIL=N echo bsub -q %s -N %s/runGEN%s.sh'%(myqueue,outDir,suffix))
    else: os.system("LSB_JOB_REPORT_MAIL=N bsub -q %s -N \'%s/runGEN%s.sh\'"%(myqueue,outDir,suffix))



