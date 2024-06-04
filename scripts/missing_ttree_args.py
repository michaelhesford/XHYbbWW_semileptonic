import subprocess, os
'''
Find snapshot jobs with missing ttrees based on the output of ttree_test.py, then output the missing job arguments to a txt file
'''

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
    'DataH_16'

groups = [
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
]
'''
groups = [
'WJetsLepHT200_16',
'WJetsLepHT200_17',
'WJetsLepHT200_18',
'WJetsLepHT400_16APV',
'WJetsLepHT800_17',
'WZ_16'
]

out = open('condor/missing_ttree_args.txt','w')
redirector = 'root://cmseos.fnal.gov/'
for g in groups:
    '''
    if not 'Data' in g:
        for year in ['16APV','16','17','18']:
            jobs = []
            #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/ttree_test/PNettree_test_{}_{}_*'.format(g,year)
            eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/ttree_test/ttree_test_{}_{}_*'.format(g,year)
            cmd = 'eos {} ls {}'.format(redirector,eos_path)
            rawFiles = subprocess.check_output(cmd, shell=True).decode()
            files = rawFiles.split('\n')
            if len(files) == 1:
                print('no files for 20{} {}'.format(year,g))
                continue
            #file name: ttree_test_ttbar-allhad_17_<JOB>_<NJOBS>.root
            njobs = int(files[0].split('.')[0].split('_')[-1])
            for f in files:
                if f == '': #this string appears at the end of the files list for some reason
                    continue
                ijob = int(f.split('_')[-2]) #this relies on a consistent format for output files
                jobs.append(ijob)
            #snapfile = open('PNetSnapshots/{}_{}_PNetSnapshot.txt'.format(g,year),'r').readlines()
            snapfile = open('snapshots/{}_{}_snapshot.txt'.format(g,year),'r').readlines()
            njobs_snap = snapfile[0].split('_')[-1].split('.')[0]
            njob_miss = 0
            for job in range(1,njobs+1):
                if not job in jobs:
                    print('Missing job number {}'.format(job))
                    print(snapfile[job-1].split('.')[0])
                    job_snap = snapfile[job-1].split('_')[-2]
                    out.write('-s %s -y %s -j %s -n %s \n'%(g,year,job_snap,njobs_snap)) 
                    njob_miss += 1
            print('Missing {} of {} total jobs for 20{} {}'.format(njob_miss,njobs,year,g))
    '''
    #else:
    if not 'Data' in g:
        jobs = []
        setname = g.split('_')[0]
        year = g.split('_')[1]
        #eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/ttree_test/PNettree_test_{}_{}_*'.format(setname,year)
        eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/ttree_test/ttree_test_{}_{}_*'.format(setname,year)
        cmd = 'eos {} ls {}'.format(redirector,eos_path)
        rawFiles = subprocess.check_output(cmd, shell=True).decode()
        files = rawFiles.split('\n')
        if len(files) == 1:
            print('no files for 20{} {}'.format(year,setname))
            continue
        #file name: ttree_test_ttbar-allhad_17_<JOB>_<NJOBS>.root
        njobs = int(files[0].split('.')[0].split('_')[-1])
        for f in files:
            if f == '': #this string appears at the end of the files list for some reason
                continue
            ijob = int(f.split('_')[-2]) #this relies on a consistent format for output files
            jobs.append(ijob)
        #snapfile = open('PNetSnapshots/{}_{}_PNetSnapshot.txt'.format(setname,year),'r').readlines()
        #snapfile = open('snapshots/{}_{}_snapshot.txt'.format(setname,year),'r').readlines()
        snapfile = open('raw_nano/{}_{}.txt'.format(setname,year),'r').readlines()
        njob_miss = 0
        njobs_snap = snapfile[0].split('_')[-1].split('.')[0]
        for job in range(1,njobs+1):
            if not job in jobs:
                print('Missing job number {}'.format(job))
                job_snap = snapfile[job-1].split('_')[-2]
                #out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,job_snap,njobs_snap))
                out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,job,njobs))
                njob_miss += 1
        print('Missing {} of {} total jobs for 20{} {}'.format(njob_miss,njobs,year,setname))
out.close()
