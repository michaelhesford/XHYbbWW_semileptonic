
'''
Gathers all rootfiles from the EOS and dumps it into the rootfiles/ directory 
'''
from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
import ROOT

def CombineCommonSets(groupname): #combines histograms of common root files
    if 'Jets' or 'QCD' in groupname:
        cmd = 'hadd -f plots/{}/{}.root plots/{}/{}*.root'.format(flavor,groupname,flavor,groupname)
    else:
        cmd = 'hadd -f plots/{}/{}.root plots/{}/{}*.root'.format(flavor,groupname,flavor,groupname)

    ExecuteCmd(cmd)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument('-g',type=str,dest='groupnames',nargs='+',
                        action='store',help='data types to combine')
    parser.add_argument('-f',type=str,dest='flavor',
                        action='store',help='name of job - used to tag files')

    args=parser.parse_args()
    groupnames=args.groupnames
    flavor=args.flavor

    for g in groupnames:
        redirector = 'root://cmseos.fnal.gov/'
        eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/plots/{}/'.format(flavor)

        if 'Jets' or 'QCD' in g:
            cmd = 'eos {} ls {}{}*'.format(redirector,eos_path,g)
        else:
            cmd = 'eos {} ls {}{}*'.format(redirector,eos_path,g)
                          
        rawFiles = subprocess.check_output(cmd, shell=True)
        files = rawFiles.split('\n')
        for fName in files:
            ExecuteCmd('xrdcp {}{}{} plots/{}/'.format(redirector, eos_path, fName, flavor))
 
        CombineCommonSets(g)

 
    
