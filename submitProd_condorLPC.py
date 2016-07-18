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
# parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
# parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,     type=int)
# parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=110,    type=int)
# parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=5,      type=int)
# parser.add_option('-a', '--eta'         ,    dest='eta'                , help='incidence eta'                , default=0.0,    type=float)
# parser.add_option('-p', '--phi'         ,    dest='phi'                , help='incidence phi angle in pi unit' , default=0.0,  type=float)
# parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0.0,    type=float)
# parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
# parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-N', '--njobs'       ,    dest='njobs'              , help='number of jobs'           , default=1, type=int)
parser.add_option('-n', '--nevtsperjob' ,    dest='nevtsperjob'        , help='number of events'         , default=10, type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default='/store/user/ntran/LDMX/Run_Jun4' )
parser.add_option('-S', '--no-submit'   ,    action="store_true"       ,  dest='nosubmit'           , help='Do not submit batch job.')
parser.add_option('--subdir'            ,    dest='subdir'             , help='directory from which you submit' , default='tmp_condor')
parser.add_option('--partype'           ,    dest='partype'            , help='particle type'             , default='electron' )
parser.add_option('--energy'            ,    dest='parenergy'          , help='particle energy'             , default=4, type=float )
#parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
(opt, args) = parser.parse_args()

def main(): 

	subdir = opt.subdir;
	if not os.path.exists(subdir): os.makedirs(subdir)
	else:  
		os.system('rm -r %s' % subdir)
		os.makedirs(subdir)
	njobs = opt.njobs; 
	nevtsperjob = opt.nevtsperjob;
	partype = opt.partype;
	parenergy = opt.parenergy;
	genfile = partype+"Gen.py";

	os.system("cp ~/geant4_workdir/bin/Linux-g++/PFCalEE ~/geant4_workdir/tmp/Linux-g++/PFCalEE/libPFCalEE.so g4env4lpc.sh userlib/lib/libPFCalEEuserlib.so generators/%s %s/." % (genfile,subdir))
	os.chdir(subdir);
	os.system("tar -cvzf inputs.tar.gz PFCalEE g4env4lpc.sh libPFCalEE.so libPFCalEEuserlib.so %s" % (genfile) );

	for i in range(njobs):

		tag = "gun_%s_%s_job%s" % ( partype, str(round(parenergy,2)), str(i) );

		# files needed: g4steer.mac.tpl, *gun.py, executable, g4env4lpc.sh
		f1n = "tmp_%s.sh" % (tag);
		f1=open(f1n, 'w')
		f1.write('#!/bin/bash \n');
		f1.write("pwd \n");
		f1.write("ls \n");
		f1.write("tar -xvzf inputs.tar.gz \n");
		f1.write('echo  what is in here? \n');
		f1.write("ls -lrt \n");
		f1.write("cat g4env4lpc.sh \n");
		f1.write("source g4env4lpc.sh \n");
		f1.write("python %s -n %s -f \"events\" -e %f \n" % (genfile, opt.nevtsperjob, float(parenergy)) )
		f1.write("ls -lrt \n");
		f1.write("PFCalEE g4steer_%s.mac 1 2 1 \n" % tag)
		f1.write("mv PFcal.root PFcal_%s.root \n" % tag)
		f1.write("xrdcp -f PFcal_%s.root root://cmseos.fnal.gov/%s/PFcal_%s.root \n" % (tag,opt.out,tag))
		f1.close();

		fsn = "g4steer_%s.mac" % tag;
		seed1 = 12345;
		seed2 = 23456;
		fs=open(fsn,'w');
		fs.write("/control/verbose 0 \n");
		fs.write("/control/saveHistory \n");
		fs.write("/run/verbose 0 \n");
		fs.write("/event/verbose 0 \n");
		fs.write("/tracking/verbose 0 \n");
		fs.write("/N03/det/setModel 2 \n");
		fs.write("/filemode/inputFilename events.lhe \n");
		fs.write("/random/setSeeds %i %i \n" % (seed1+i,seed2+i));
		fs.write("/run/initialize \n");
		fs.write("/run/beamOn %i \n" % opt.nevtsperjob);
		fs.close();

		f2n = "tmp_%s.jdl" % (tag);
		outtag = "out_%s_$(Cluster)" % (tag)
		f2=open(f2n, 'w')
		f2.write("universe = vanilla \n");
		f2.write("Executable = %s \n" % (f1n) );
		f2.write('Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )\n');
		f2.write("request_disk = 10000000\n");
		f2.write("request_memory = 10000\n");
		f2.write("Should_Transfer_Files = YES \n");
		f2.write("Transfer_Input_Files = inputs.tar.gz,g4steer_%s.mac \n" % (tag));
		f2.write("WhenToTransferOutput  = ON_EXIT_OR_EVICT \n");
		f2.write("Output = "+outtag+".stdout \n");
		f2.write("Error = "+outtag+".stderr \n");
		f2.write("Log = "+outtag+".log \n");
		f2.write("Notification    = Error \n");
		f2.write("x509userproxy = $ENV(X509_USER_PROXY) \n")
		f2.write("Queue 1 \n");
		f2.close();

		if not opt.nosubmit: os.system("condor_submit %s" % (f2n));

	########## end of loop
	os.chdir("../.");



#----------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main();
#----------------------------------------------------------------------------------------------------------------

