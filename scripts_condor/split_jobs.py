#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

filesPerJob = 10
tmpJobFileCount = 0
nJobs = 0

inputListFile = open("/uscms_data/d3/dildick/work/DisplacedHeavyNeutralLeptonAnalysis/CMSSW_10_6_20/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p17/MC_Fall18/v1/sixie/HNL_testpoint1.txt", "r")

tmpOutputListFile = open("input_list/" + "input_list_" + str(nJobs) + ".txt","w")
for line in inputListFile:

    #open list file for new job
    if tmpJobFileCount >= filesPerJob:
        tmpOutputListFile.close()
        tmpJobFileCount = 0
        nJobs = nJobs + 1
        tmpOutputListFile = open("input_list/" + "input_list_" + str(nJobs) + ".txt","w")

    #write input file into job list file
    tmpOutputListFile.write(line)
    tmpJobFileCount += 1

tmpOutputListFile.close()

## make tarball with all input files
os.system("cd input_list; tar czf input_list.tgz input_list_*.txt; cp input_list.tgz ../; cd ../")

inputListFile.close()
