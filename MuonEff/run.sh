#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh 
cd /cms/ldap_home/jieun/CMSSW_10_1_5/src 
eval `scramv1 runtime -sh` 
cd /cms/ldap_home/jieun/MuonSW/MuonEff

root -l -b < x_file.C >& process.log  
