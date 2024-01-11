'''
Grabs specified histograms and dumps them in desired directory
'''
from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
import ROOT
    
def combine_ttbar(flavor):
    variations = ['_nom','_JES_up','_JES_down','_JER_up','_JER_down','_JMS_up','_JMS_down','_JMR_up','_JMR_down']
    for var in variations:
       semilep = 'plots/{}/ttbar-semilep{}_{}.root'.format(flavor,var,flavor)
       allhad = 'plots/{}/ttbar-allhad{}_{}.root'.format(flavor,var,flavor)
       ExecuteCmd('hadd -f plots/{}/ttbar{}_{}.root {} {}'.format(flavor,var,flavor,allhad,semilep))
   
def CombineCommonSets(groupname,make_Run2): #combines histograms of common root files
    if make_Run2:
        cmd = 'hadd -f plots/{}/{}_{}.root plots/{}/{}*.root'.format(flavor,groupname,flavor,flavor,groupname)
        ExecuteCmd(cmd)
    else:   
        for year in ['16','16APV','17','18']:
            cmd = 'hadd -f plots/{}/{}_{}_{}.root plots/{}/{}_{}_*.root'.format(flavor,groupname,year,flavor,flavor,groupname,year)
            ExecuteCmd(cmd)
                          
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument('-g',type=str,dest='groupnames',nargs='+',
                        default='all',action='store',help='data types to combine')
    parser.add_argument('-f',type=str,dest='flavor',required=False,
                        action='store',help='name of job - used to tag files')
    parser.add_argument('--combine_years',action='store_true',help='if used, make full run 2')
    parser.add_argument('--JEC_groups',action='store_true',help='use when processing JEC studies') #DEPRECATED

    all_groups = [
        'XHY-4000-2500',
        'XHY-3500-1800',
        'XHY-3000-1400',
        'XHY-2500-1300',
        'XHY-2000-1000',
        'XHY-1600-700',
        'XHY-1000-500',
        'QCD',
        'WJetsHT', #add "HT" to avoid grabbing semileptonic samples as well
        'ZJetsHT',
        'ttbar-allhad',
        'ttbar-semilep',
        'WJetsToLNu'
        #'WWto4Q',
        #'WWto1L1Nu2Q',
        #'WZto1L1Nu2Q',
        #'WZto2L2Q',
        #'ZZto4Q'
        #'ZZto2L2Q'
    ]

    args = parser.parse_args()
    flavor = args.flavor
    groupnames = args.groupnames if not args.groupnames == 'all' else all_groups
    if args.JEC_groups:
        g_start = groupnames
        groupnames = ['ttbar-allhad_nom','ttbar-semilep_nom']
        for g in g_start:
            for syst in ['JES','JER','JMS','JMR']:
                for var in ['up','down']:
                    groupnames.append('{}_{}_{}'.format(g,syst,var))
    combine_years = args.combine_years    
    
    for g in groupnames:
        redirector = 'root://cmseos.fnal.gov/'
        eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/plots/{}/'.format(flavor)

        cmd = 'eos {} ls {}{}*'.format(redirector,eos_path,g)
                         
        rawFiles = subprocess.check_output(cmd, shell=True)
        files = rawFiles.split('\n')
        for fName in files:
            ExecuteCmd('xrdcp {}{}{} plots/{}/'.format(redirector, eos_path, fName, flavor))
 
        CombineCommonSets(g,combine_years)
    
    if args.JEC_groups:
        combine_ttbar(flavor)

   
