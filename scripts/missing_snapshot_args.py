import subprocess, os

''' 
'XHY-1000-500',
'XHY-1600-700',
'XHY-2000-100',
'XHY-3000-1400',
'XHY-3500-1800',
'XHY-4000-2500',
'QCDHT700',
'QCDHT1000',
'QCDHT1500',
'QCDHT2000',
'WJetsHT400',
'WJetsHT600',
'WJetsHT800',
'WJetsToLNuHT200',
'WJetsToLNuHT400',
'WJetsToLNuHT600',
'WJetsToLNuHT800',
'WJetsToLNuHT1200',
'WJetsToLNuHT2500',
'ZJetsHT400',
'ZJetsHT600',
'ZJetsHT800',
'ttbar-allhad',
'ttbar-semilep',
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
'SingleMuonDataA_18',
'SingleMuonDataB_16APV',
'SingleMuonDataB_17',
'SingleMuonDataB_18',
'SingleMuonDataC_16APV',
'SingleMuonDataC_17',
'SingleMuonDataC_18',
'SingleMuonDataD_16APV',
'SingleMuonDataD_17',
'SingleMuonDataD_18',
'SingleMuonDataE_16APV',
'SingleMuonDataE_17',
'SingleMuonDataF_16',
'SingleMuonDataF_16APV',
'SingleMuonDataF_17',
'SingleMuonDataG_16',
'SingleMuonDataH_16',
'EGammaDataA_18',
'EGammaDataB_18',
'EGammaDataC_18',
'EGammaDataD_18',
'SingleElectronDataB_17',
'SingleElectronDataC_17',
'SingleElectronDataD_17',
'SingleElectronDataE_17',
'SingleElectronDataF_17',
'SingleElectronDataF_16',
'SingleElectronDataG_16',
'SingleElectronDataH_16',
'SingleElectronDataB_16APV',
'SingleElectronDataC_16APV',
'SingleElectronDataD_16APV',
'SingleElectronDataE_16APV',
'SingleElectronDataF_16APV'
'''

groups = [
'WJetsHT200',
'WJetsHT400',
'WJetsHT600',
'WJetsHT800',
'WJetsLepHT200',
'WJetsLepHT400',
'WJetsLepHT600',
'WJetsLepHT800',
'WJetsLepHT1200',
'WJetsLepHT2500',
'ZJetsHT200',
'ZJetsHT400',
'ZJetsHT600',
'ZJetsHT800',
'WW',
'WZ',
'ZZ'
]

out = open('condor/missing_snapshot_args.txt','w')
#out = open('condor/missing_PNETsnapshot_args.txt','w')
redirector = 'root://cmseos.fnal.gov/'
for g in groups:
    if not 'Data' in g:
        for year in ['16APV','16','17','18']:
            jobs = []
            eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/snapshots/XHYbbWWsnapshot_{}_{}_*'.format(g,year)
            #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/PNetSnapshots/PNetCRSnapshot_{}_{}_*'.format(g,year)
            cmd = 'eos {} ls {}'.format(redirector,eos_path)
            rawFiles = subprocess.check_output(cmd, shell=True).decode()
            files = rawFiles.split('\n')
            if len(files) == 1:
                print('no files for 20{} {}'.format(year,g))
                continue
            #file name: XHYbbWWsnapshot_ttbar-allhad_17_<JOB>_<NJOBS>.root
            njobs = int(files[0].split('.')[0].split('_')[-1])
            for f in files:
                if f == '': #this string appears at the end of the files list for some reason
                    continue
                ijob = int(f.split('_')[-2]) #this relies on a consistent format for output files
                jobs.append(ijob)
            njob_miss = 0
            for job in range(1,njobs+1):
                if not job in jobs:
                    print('Missing job number {}'.format(job))
                    out.write('-s %s -y %s -j %s -n %s \n'%(g,year,job,njobs)) 
                    njob_miss += 1
            print('Missing {} of {} total jobs for 20{} {}'.format(njob_miss,njobs,year,g))
    else:
        jobs = []
        setname = g.split('_')[0]
        year = g.split('_')[1]
        eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/snapshots/XHYbbWWsnapshot_{}_{}_*'.format(setname,year)
        #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/PNetSnapshots/PNetCRSnapshot_{}_{}_*'.format(setname,year)
        cmd = 'eos {} ls {}'.format(redirector,eos_path)
        rawFiles = subprocess.check_output(cmd, shell=True).decode()
        files = rawFiles.split('\n')
        if len(files) == 1:
            print('no files for 20{} {}'.format(year,setname))
            continue
        #file name: XHYbbWWsnapshot_ttbar-allhad_17_<JOB>_<NJOBS>.root
        njobs = int(files[0].split('.')[0].split('_')[-1])
        for f in files:
            if f == '': #this string appears at the end of the files list for some reason
                continue
            ijob = int(f.split('_')[-2]) #this relies on a consistent format for output files
            jobs.append(ijob)
        njob_miss = 0
        for job in range(1,njobs+1):
            if not job in jobs:
                print('Missing job number {}'.format(job))
                out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,job,njobs))
                njob_miss += 1
        print('Missing {} of {} total jobs for 20{} {}'.format(njob_miss,njobs,year,setname))
out.close()
