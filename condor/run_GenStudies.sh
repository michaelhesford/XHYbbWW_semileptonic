#!/bin/bash
echo "Run script starting"
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic.tgz ./
export SCRAM_ARCH=slc7_amd64_gcc820
scramv1 project CMSSW CMSSW_11_1_4
tar -xzvf XHYbbWW_semileptonic.tgz
rm XHYbbWW_semileptonic.tgz
rm *.root

mkdir tardir; cp tarball.tgz tardir/; cd tardir/
tar -xzf tarball.tgz; rm tarball.tgz
cp -r * ../CMSSW_11_1_4/src/semileptonic; cd ../CMSSW_11_1_4/src/
echo 'IN RELEASE'
pwd
ls
eval `scramv1 runtime -sh`
rm -rf timber-env
python -m virtualenv timber-env
source timber-env/bin/activate
cd TIMBER
source setup.sh
cd ../semileptonic

echo python Gen_studies.py $*
python Gen_studies.py $*

xrdcp -f *.root root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic/plots/GenStudies/
