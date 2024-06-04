from glob import glob
from TIMBER.Tools.Common import ExecuteCmd
import subprocess, os

redirector = 'root://cmseos.fnal.gov/'
eos_path = '/store/user/mhesford/XHYbbWW_semileptonic/event_counts/'
    
cmd = 'eos {} ls {}*count.txt'.format(redirector,eos_path)

rawFiles = subprocess.check_output(cmd, shell=True).decode()
files = rawFiles.split('\n')
for fName in files:
    #Throwing all the output histograms from the selection script into a directory called "raw_selection"
    ExecuteCmd('xrdcp {}{}{} event_counts/'.format(redirector, eos_path, fName)) 
