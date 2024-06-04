from TIMBER.Tools.Common import ExecuteCmd
import sys

das = {
    '16APV':{
        'SinglePhotonDataBv1' : '/SinglePhoton/Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SinglePhotonDataBv2' : '/SinglePhoton/Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SinglePhotonDataC' : '/SinglePhoton/Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v4/NANOAOD',
        'SinglePhotonDataD' : '/SinglePhoton/Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SinglePhotonDataE' : '/SinglePhoton/Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SinglePhotonDataF' : '/SinglePhoton/Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD'
    },
    '16':{
        'SinglePhotonDataF' : '/SinglePhoton/Run2016F-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SinglePhotonDataG' : '/SinglePhoton/Run2016G-UL2016_MiniAODv2_NanoAODv9-v2/NANOAOD',
        'SinglePhotonDataH' : '/SinglePhoton/Run2016H-UL2016_MiniAODv2_NanoAODv9-v1/NANOAOD'
    },
    '17':{
        'SinglePhotonDataB' : '/SinglePhoton/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SinglePhotonDataC' : '/SinglePhoton/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SinglePhotonDataD' : '/SinglePhoton/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SinglePhotonDataE' : '/SinglePhoton/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD',
        'SinglePhotonDataF' : '/SinglePhoton/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD'
    }
}

def GetFiles(das_name,setname,year):
    ExecuteCmd('dasgoclient -query "file dataset=%s" > %s_%s_temp.txt'%(das_name,setname,year),'dryrun' in sys.argv)
    f = open('%s_%s_temp.txt'%(setname,year),'r')
    fout = open('raw_nano/%s_%s.txt'%(setname,year),'w')
    for l in f.readlines():
        fout.write('root://cmsxrootd.fnal.gov/'+l)
    f.close()
    fout.close()
    if 'dryrun' not in sys.argv:
        ExecuteCmd('rm %s_%s_temp.txt'%(setname,year),'dryrun' in sys.argv)

latex_lines = {"16":[],"16APV":[],"17":[],"18":[]}
for year in das.keys():
    for setname in das[year].keys():
        GetFiles(das[year][setname],setname,year)
        latex_lines[year].append('| %s | %s |'%(setname,das[year][setname]))
            
for y in sorted(latex_lines.keys()):
    print ('\n20%s'%y)
    print ('| Setname | DAS location |')
    print ('|---------|--------------|')
    for l in latex_lines[y]: print (l)
