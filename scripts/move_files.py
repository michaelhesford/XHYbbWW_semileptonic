from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
import ROOT
    
if __name__ == '__main__':

    all_groups = [
    'DataA_18',
    'DataB_16APV',
    'DataB_17',
    'DataB_18',
    'DataC_16APV',
    'DataC_17',
    'DataC_18',
    'DataD_16APV',
    'DataD_17',
    'DataD_18',
    'DataE_16APV',
    'DataE_17',
    'DataF_16',
    'DataF_16APV',
    'DataF_17',
    'DataG_16',
    'DataH_16',
    'ttbar-allhad',
    'ttbar-semilep'
    ]

    for g in all_groups:
        redirector = 'root://cmseos.fnal.gov/'
        eos_path_start = '/store/user/mhesford/XHYbbWW_semileptonic/snapshots/'
        eos_path_end = '/store/user/mhesford/XHYbbWW_semileptonic/snapshots_old/'

        cmd = 'eos {} ls {}XHYbbWWsnapshot_{}*'.format(redirector,eos_path_start,g)
                         
        rawFiles = subprocess.check_output(cmd, shell=True).decode()
        files = rawFiles.split('\n')
        for fName in files:
            ExecuteCmd('xrdcp {}{}{} {}{}'.format(redirector, eos_path_start, fName, redirector, eos_path_end))
            #ExecuteCmd('eos root://cmseos.fnal.gov rm {}{}'.format(eos_path_start,fName))
 
   
