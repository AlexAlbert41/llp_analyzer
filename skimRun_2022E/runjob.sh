#!/bin/bash
date
MAINDIR=`pwd`
ls
voms-proxy-info --all
#CMSSW from scratch (only need for root)
export CWD=${PWD}
export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW CMSSW_10_6_20
mv SkimNtuples CMSSW_10_6_20/src
cd CMSSW_10_6_20/src
eval `scramv1 runtime -sh` # cmsenv
echo "Untar JEC:" 
echo "After Untar: "
ls
python ${MAINDIR}/convertList.py -i ${MAINDIR}/tmp_input_list_$1.txt
SkimNtuples .${MAINDIR}/tmp_input_list_$1.txt ${MAINDIR} "skim" ""
echo "coping to eos: "+xrdcp -f ${MAINDIR}/displacedJetMuon_ntupler_*_skim.root root://cmseos.fnal.gov//store/user/amalbert/MDSTriggerEff/./skimRun_2022E  
xrdcp -f ${MAINDIR}/displacedJetMuon_ntupler_*_skim.root root://cmseos.fnal.gov//store/user/amalbert/MDSTriggerEff/./skimRun_2022E
echo "Inside $MAINDIR:"
ls
echo "DELETING..."
rm -rf CMSSW_10_6_20
rm -rf *.pdf *.C core*
cd $MAINDIR  
echo "remove output local file"  
rm -rf *.root 
ls
date
