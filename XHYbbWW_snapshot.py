import ROOT, time
ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
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
selection = XHYbbWW(filename,ijob,njobs)
selection.BasicKinematics()
out=selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale()) #renormalize events
selection.Snapshot(out)

print('%s sec'%(time.time()-start))
