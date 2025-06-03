#!/bin/bash
echo "Run script starting"
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic.tgz ./
export SCRAM_ARCH=el8_amd64_gcc10
scramv1 project CMSSW CMSSW_12_3_5
tar -xzvf XHYbbWW_semileptonic.tgz
rm XHYbbWW_semileptonic.tgz

mkdir tardir; cp tarball.tgz tardir/; cd tardir/
tar -xzf tarball.tgz; rm tarball.tgz
cp -r * ../CMSSW_12_3_5/src/semileptonic/; cd ../CMSSW_12_3_5/src/
echo 'IN RELEASE'
pwd
ls
eval `scramv1 runtime -sh`
rm -rf timber-env
python3 -m venv timber-env

cat <<EOT >> timber-env/bin/activate

export BOOSTPATH=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.78.0-0d68c45b1e2660f9d21f29f6d0dbe0a0/lib
export LIBPATH=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/python3/3.9.6-67e5cf5b4952101922f1d4c8474baa39/lib/
if grep -q '\${BOOSTPATH}' <<< '\${LD_LIBRARY_PATH}'
then
  echo 'BOOSTPATH already on LD_LIBRARY_PATH'
else
  export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${LIBPATH}:\${BOOSTPATH}
  echo 'BOOSTPATH added to PATH'
fi
EOT

source timber-env/bin/activate
cd TIMBER
source setup.sh
cd ../semileptonic

echo python PlotDistributions.py $*
python PlotDistributions.py $*


xrdcp -f *distributions.root root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic/plots/distributions/
