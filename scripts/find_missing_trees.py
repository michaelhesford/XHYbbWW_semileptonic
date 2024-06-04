import subprocess
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-c', type=str, dest='cluster', required=True,
                    action='store', help='cluster of jobs to check')
args = parser.parse_args()

out = open('condor/missing_snapshot_trees.txt','w')
cmd = "grep -H -R 'WARNING: The following file does NOT contain an Events TTree, skipping.' logs/output_{}*.stdout | cut -d: -f1".format(args.cluster)
cmd2 = "grep -H -R 'WARNING: The following file contains an empty Events TTree, skipping.' logs/output_{}*.stdout | cut -d: -f1".format(args.cluster)
#rawFiles = subprocess.check_output(cmd, shell=True)
files = rawFiles.split('\n')
#rawFiles2 = subprocess.check_output(cmd2, shell=True)
#files.extend(rawFiles2.split('\n'))
unique_files = list(set(files))
for f in unique_files:
    if f == '': #this string appears at the end of the files list for some reason
        continue
    print(f)
    f_open = open(f,'r')
    lines = f_open.readlines()
    for i in range(len(lines)):
        if 'WARNING: The following file does NOT contain an Events TTree, skipping.' in lines[i] or 'WARNING: The following file contains an empty Events TTree, skipping.' in lines[i]:
            line_num = i+1
            line = lines[line_num]
            contents = line.split('/')[-1].split('.')[0].split('_')
            print(contents)
            setname = contents[1]
            year = contents[2]
            ijob = contents[3]
            njobs = contents[4]
            out.write('-s %s -y %s -j %s -n %s \n'%(setname,year,ijob,njobs))            
out.close()

