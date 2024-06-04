#script runs through all the files in an eos directory and finds any missing jobs

import subprocess
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-j', type=int, dest='njobs',required=False,
                    action='store', help='number of expected jobs')
parser.add_argument('-p', type=str, dest='path', required=True,
                    action='store', help='eos path to jobs')     
args = parser.parse_args()
eos_path = args.path
#njobs = args.njobs

jobs=[]
redirector = 'root://cmseos.fnal.gov/'
cmd = 'eos {} ls {}*'.format(redirector,eos_path)
rawFiles = subprocess.check_output(cmd, shell=True).decode()
files = rawFiles.split('\n')
#file name: XHYbbWW_pileup_ttbar-allhad_17_<JOB>_<NJOBS>.root
njobs = int(files[0].split('.')[0].split('_')[-1])
for f in files:
    print(f)
    if f == '': #this string appears at the end of the files list for some reason
        continue
    ijob = int(f.split('_')[-2]) #this relies on a consistent format for output files
    jobs.append(ijob)

njob_miss = 0
for job in range(1,njobs+1):
    if not job in jobs:
        print('Missing job number {}'.format(job))
        njob_miss += 1
print('Missing {} of {} total jobs'.format(njob_miss,njobs))


