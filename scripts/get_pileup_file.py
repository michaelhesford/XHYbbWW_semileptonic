from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
import ROOT


if __name__ == '__main__':
    redirector = 'root://cmseos.fnal.gov/'
    eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/'
    cmd = 'eos {} ls {}'.format(redirector,eos_path)
    
    rawFiles = subprocess.check_output(cmd+'pileup/XHYbbWW_pileup_WJetsToLNu**', shell=True)
    files = rawFiles.split('\n')
    for fName in files:
        ExecuteCmd('xrdcp {}{}pileup/{} /uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_11_1_4/src/semileptonic/pileup/'.format(redirector, eos_path, fName))
    '''
    #now get pileup for 2017 ttbar-allhad, which is stored separately in a bunch of files, due to issues processing the full MC together
    for ttbar_file in subprocess.check_output(cmd+'pileup_ttbar-allhad17/*', shell=True).split('\n'):
        ExecuteCmd('xrdcp {}{}pileup_ttbar-allhad17/{} /uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_11_1_4/src/semileptonic/pileup/'.format(redirector,eos_path,ttbar_file))
    ExecuteCmd('hadd -f pileup/XHYbbWW_pileup_ttbar-allhad_17.root pileup/XHYbbWW_pileup_ttbar-allhad_17*')
    ExecuteCmd('rm pileup/XHYbbWW_pileup_ttbar-allhad_17_*') #remove all individual files
    '''
    #Now make full pileup file 
    ExecuteCmd('hadd -f XHYbbWWpileup.root pileup/*pileup*')


