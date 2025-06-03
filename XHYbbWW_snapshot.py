import ROOT, time
ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
from collections import OrderedDict
from TIMBER.Tools.Common import CompileCpp
from XHYbbWW_class import XHYbbWW

parser = ArgumentParser()
parser.add_argument('-s', type=str, dest='setname',
                    action='store', required=True,
                    help='Setname to process.')
parser.add_argument('-y', type=str, dest='year',
                    action='store', required=True,
                    help='Year of set (16, 17, 18).')
parser.add_argument('-j', type=int, dest='ijob',
                    action='store', default=1,
                    help='Job number')
parser.add_argument('-n', type=int, dest='njobs',
                    action='store', default=1,
                    help='Number of jobs')
args = parser.parse_args()

start = time.time()

setname=args.setname
year=args.year
ijob=args.ijob
njobs=args.njobs

filename='raw_nano/{}_{}.txt'.format(setname,year)
selection = XHYbbWW(filename,ijob,njobs,use_xrd_global=True)
selection.BasicKinematics()
out=selection.ApplyStandardCorrections(snapshot=True,do_ak4JEC=True)
'''
#create cut flow table
outFile = ROOT.TFile.Open('{}_{}_{}_SnapCutFlow.root'.format(setname,year,ijob),'RECREATE')
outFile.cd()

cutflowInfo = OrderedDict([
    ('start',selection.NSTART),
    ('flags',selection.NFLAGS),
    ('lep. pre-sel.',selection.LEPPRE),
    ('jet pre-sel',selection.JETPRE),
    ('MET>25',selection.METPT)	
])

nLabels = len(cutflowInfo)
hCutflow = ROOT.TH1F('cutflow', 'Number of events after each cut', nLabels, 0.5, nLabels+0.5)
nBin = 1
for label, value in cutflowInfo.items():
    hCutflow.GetXaxis().SetBinLabel(nBin, label)
    hCutflow.AddBinContent(nBin, value)
    nBin += 1
hCutflow.Write()
outFile.Close()
'''

selection.Snapshot(out)

print('%s sec'%(time.time()-start))
