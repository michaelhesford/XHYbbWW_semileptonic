'''
I mess with this script a lot, but basically this can be used to check if any of the snapshot files are missing a particular branch
'''

import ROOT 
import glob
from TIMBER.Tools.Common import ExecuteCmd

data_files = ['Data_16APV_PNetSnapshot.txt','Data_16_PNetSnapshot.txt','Data_17_PNetSnapshot.txt','Data_18_PNetSnapshot.txt']

out = open('condor/missing_rawFactor.txt','w')
for filename in glob.glob('snapshots/*.txt'):
    if filename.split('/')[-1].startswith('Data') and filename.split('/')[-1] not in data_files: continue # Avoid redundancies
    data_era = filename.split('/')[-1].split('_')[0]
    year = filename.split('/')[-1].split('_')[1]
    #out1 = open(f'PNetSnapshots/{data_era}1_{year}_PNetSnapshot.txt','w')
    #out2 = open(f'PNetSnapshots/{data_era}2_{year}_PNetSnapshot.txt','w')
    #out1 = open(f'snapshots/{data_era}1_{year}_snapshot.txt','w')
    #out2 = open(f'snapshots/{data_era}2_{year}_snapshot.txt','w')
    print(filename)
    for job in open(filename,'r').readlines():
        job = job.strip()
        try: 
            inFile = ROOT.TFile.Open(job,'READ')
        except Exception as e:
            print(e)
        else:
            if inFile.GetListOfKeys().Contains('Events'):
                ttree = inFile.Get('Events')
                if not ttree.GetListOfBranches().FindObject("Jet_rawFactor"): 
                    print(f'Branch not found in {job}')
                    info = job.split('/')[-1].split('.')[0].split('_')
                    setname = info[1]
                    year = info[2]
                    ijob = info[3]
                    njobs = info[4]
                    out.write(f'-s {setname} -y {year} -j {ijob} -n {njobs}\n')
                    #out2.write(job+'\n')
                else:
                    #out1.write(job+'\n')
                    x = 1
            else:
                print("The Events tree is not found in file {}".format(job.split('/')[-1]))
        inFile.Close()
    #out1.close()
    #out2.close()
    #if len(open(f'PNetSnapshots/{data_era}2_{year}_PNetSnapshot.txt','r').readlines()) == 0: # No reason to divide files, just delete out1 and out2
    #    ExecuteCmd(f'rm PNetSnapshots/{data_era}1_{year}_PNetSnapshot.txt')
    #    ExecuteCmd(f'rm PNetSnapshots/{data_era}2_{year}_PNetSnapshot.txt')
    #if len(open(f'snapshots/{data_era}2_{year}_snapshot.txt','r').readlines()) == 0: # No reason to divide files, just delete out1 and out2
    #    ExecuteCmd(f'rm snapshots/{data_era}1_{year}_snapshot.txt')
    #    ExecuteCmd(f'rm snapshots/{data_era}2_{year}_snapshot.txt')
out.close()
