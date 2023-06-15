#create tarball of working directory and store in eos space
cd $CMSSW_BASE/../
tar --exclude-caches-all --exclude-vcs --exclude-caches-all --exclude-vcs -cvzf XHYbbWW_semileptonic.tgz CMSSW_11_1_4 --exclude=tmp --exclude=".scram" --exclude=".SCRAM" --exclude=CMSSW_11_1_4/src/timber-env --exclude=CMSSW_11_1_4/src/semileptonic/logs --exclude=CMSSW_11_1_4/src/TIMBER/docs --exclude=CMSSW_11_1_4/src/semileptonic/plot
xrdcp -f XHYbbWW_semileptonic.tgz root://cmseos.fnal.gov//store/user/$USER/XHYbbWW_semileptonic.tgz
cd $CMSSW_BASE/src/semileptonic/

#--exclude=CMSSW_11_1_4/src/semileptonic/*.root

