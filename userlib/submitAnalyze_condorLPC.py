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
parser.add_option('-N', '--njobs'       ,    dest='njobs'              , help='number of jobs'           , default=1, type=int)
parser.add_option('-n', '--nevtsperjob' ,    dest='nevtsperjob'        , help='number of events'         , default=10, type=int)
parser.add_option('--subdir'            ,    dest='subdir'             , help='directory from which you submit' , default='tmp_condor')
parser.add_option('--idir'              ,    dest='idir'               , help='directory from which you submit' , default='/store/user/ntran/LDMX/Run_Jun4')
parser.add_option('--odir'              ,    dest='odir'               , help='directory from which you submit' , default='/store/user/ntran/LDMX/Run_Jun4_analyzed')
parser.add_option('--partype'           ,    dest='partype'            , help='particle type'             , default='electron' )
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default='/store/user/ntran/LDMX/Run_Jun4' )
parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
parser.add_option('-g', '--gun'         ,    action="store_true"       ,  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true"       ,  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()



def main(): 

	subdir = opt.subdir;
	if not os.path.exists(subdir): os.makedirs(subdir)
	else:  
			os.system('rm -r %s' % subdir)
			os.makedirs(subdir)

	os.chdir(subdir);
	os.system("cp ../../g4env4lpc.sh ../bin/analyze_sampling ../lib/libPFCalEEuserlib.so .")
	os.system("tar -cvzf inputs.tar.gz g4env4lpc.sh analyze_sampling libPFCalEEuserlib.so")

	listOfFiles = os.listdir("/eos/uscms%s/" % opt.idir);

	ctr = 0;
	for f in listOfFiles:

		f1n = "tmp_%s.sh" % (str(ctr));
		f1=open(f1n, 'w')
		f1.write('#!/bin/bash \n');
		f1.write("tar -xvzf inputs.tar.gz \n");
		f1.write("source g4env4lpc.sh \n");
		f1.write("ls -lrt \n");
		f1.write("xrdcp root://cmseos.fnal.gov/%s/%s . \n" % (opt.idir,f))
		f1.write("./analyze_sampling %s \n" % (f))
		f1.write("rm %s \n" % (f))
		f1.write("mv analyzed_tuple.root analyzed_tuple_%s.root \n" % str(ctr))
		f1.write("xrdcp -f analyzed_tuple_%s.root root://cmseos.fnal.gov/%s/analyzed_tuple_%s.root \n" % (str(ctr),opt.odir,str(ctr)) )
		f1.close();

		f2n = "tmp_%s.jdl" % (str(ctr));
		outtag = "out_%s_$(Cluster)" % (str(ctr))
		f2=open(f2n, 'w')
		f2.write("universe = vanilla \n");
		f2.write("Executable = %s \n" % (f1n) );
		f2.write('Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )\n');
		f2.write("request_disk = 10000000\n");
		f2.write("request_memory = 10000\n");
		f2.write("Should_Transfer_Files = YES \n");
		f2.write("Transfer_Input_Files = inputs.tar.gz \n");
		f2.write("WhenToTransferOutput  = ON_EXIT_OR_EVICT \n");
		f2.write("Output = "+outtag+".stdout \n");
		f2.write("Error = "+outtag+".stderr \n");
		f2.write("Log = "+outtag+".log \n");
		f2.write("Notification    = Error \n");
		f2.write("x509userproxy = $ENV(X509_USER_PROXY) \n")
		f2.write("Queue 1 \n");
		f2.close();

		if not opt.nosubmit: os.system("condor_submit %s" % (f2n));

		ctr += 1;
		# if ctr > 0: break;

	########## end of loop
	os.chdir("../.");


#----------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	main();
#----------------------------------------------------------------------------------------------------------------

