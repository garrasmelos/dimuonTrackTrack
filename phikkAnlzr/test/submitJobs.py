#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

print 
print 'START'
print 
########   YOU ONLY NEED TO FILL THE AREA BELOW   #########
########   customization  area #########
NumberOfJobs= 117              # number of jobs to be submitted
interval = 1                    # number files to be processed in a single job, take care to split your file so that you run on all files. 
                                # The last job might be with smaller number of files (the ones that remain).
OutputFileNames = "dimuonTrackTrack"  # base of the output file name, they will be saved in res directory
ScriptName = "ConfFile_cfg.py" # script to be used with cmsRun
FileList = "List.txt"           # list with all the file directories
#DirInEOS = "root://eosuser.cern.ch//eos/user/g/garamire/Temporal/store/"  #this dir has to exist in eos otherwise set to "0", next linei
DirInEOS = "root://eoscms//eos/cms/store/group/phys_bphys/garamire/MuOnia/dimuonTrackTrack/2015/"
#DirInEOS = "root://eoscms//eos/cms/store/group/dpg_rpc/comm_rpc/Sandbox/garamire/HSCP_MC_GEN-SIM-Feb2017-M247/HSCP_MC_GEN-SIM-Feb2017-M247/170215_203904/0000/"
UseEOS = 0                      # use or not eos, if not used, output will be transfer to current cwd area
queue = "8nh"                   # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
########   customization end   #########

path = str(os.getcwd())
print
print 'do not worry about folder creation:'
os.system("rm -rf tmp")
os.system("mkdir -p tmp")
tmpdir = path+"/tmp/"
if not UseEOS :
   os.system("mkdir -p res")
   resdir = path+"/res/"
print

##### loop for creating and sending jobs #####
for x in range(1, int(NumberOfJobs)+1):
   ##### creates directory and file list for job #######
   os.system("mkdir -p tmp/"+str(x))
   os.chdir("tmp/"+str(x))
   os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' ../../"+FileList+" > list.txt ")
   
   ##### creates jobs #######
   with open('job.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'WORKDIR ' ${PWD}\n")
      fout.write("export SCRAM_ARCH=slc6_amd64_gcc493\n")
      fout.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
      fout.write("cd "+path+"\n")
      fout.write("eval `scram runtime -sh`\n")
      fout.write("cd -\n")
      fout.write("cp -v "+tmpdir+str(x)+"/list.txt .\n")
      fout.write("cp -v "+path+"/ConfFile_cfg.py .\n")
      fout.write("time cmsRun "+ScriptName+" outputFile='rootfile.root' inputFiles_clear inputFiles_load='list.txt' >& log\n")
      fout.write("gzip -9 log\n")
      if UseEOS:
         fout.write("xrdcp -v -f rootfile.root "+DirInEOS+OutputFileNames+"_"+str(x)+".root \n")
         fout.write("xrdcp -v -f log.gz "+DirInEOS+"log_"+str(x)+".log.gz \n")
      else:
         fout.write("cp -v rootfile.root "+resdir+OutputFileNames+"_"+str(x)+".root\n")
         fout.write("cp -v log.gz "+resdir+"log_"+str(x)+".log.gz\n")
      fout.write("echo 'END---------------'\n")
      fout.write("echo\n")
   os.system("chmod 755 job.sh")
   
   ###### sends bjobs ######
   os.system("bsub -q "+queue+" -o logs job.sh")
   time.sleep(0.25) # Time in seconds.
   print "job number " + str(x) + " submitted"
   
   os.chdir("../..")
   
print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
