import subprocess, os

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
    'DataH_16'
]
out = open('condor/missing_template_test_args.txt','w')
redirector = 'root://cmseos.fnal.gov/'
for g in groups:
    if not 'Data' in g:
        for year in ['16APV','16','17','18']:
            jobs = []
            eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/plots/PNetTemplates2D_test/PNetTemplates2D_Hbb_{}_{}_*'.format(g,year)
            cmd = 'eos {} ls {}'.format(redirector,eos_path)
            rawFiles = subprocess.check_output(cmd, shell=True).decode()
            files = rawFiles.split('\n')
            if len(files) == 1:
                print('no files for 20{} {}'.format(year,g))
                continue
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
        eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/plots/PNetTemplates2D_test/PNetTemplates2D_{}_{}_*'.format(setname,year)
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
