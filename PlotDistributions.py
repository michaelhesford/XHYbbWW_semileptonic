from XHYbbWW_class import XHYbbWW
#from TIMBER.Analyzer import analyzer
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser

parser = ArgumentParser()


parser.add_argument('-s', type=str, dest='setname',
                    action='store', required=True,
                    help='Setname to process.')
parser.add_argument('-y', type=int, dest='year',
                    action='store', required=True)
parser.add_argument('-j', type=int, dest='ijob',
                    action='store', default=1,
                    help='Job number')
parser.add_argument('-n', type=int, dest='njobs',
                    action='store', default=1,
                    help='Number of jobs')

args = parser.parse_args()

filename = 'raw_nano/Signal/{}_{}.txt'.format(args.setname,args.year)
selection = XHYbbWW(filename, args.ijob, args.njobs)

variables=['MET_pt','MET_phi','nElectrons','nMuons','FatJet_msoftdrop[0]','FatJet_msoftdrop[1]','DeltaPhi_jets','DeltaPhi_Electron_Higgs','DeltaPhi_Muon_Higgs','DeltaPhi_Electron_Wqq','DeltaPhi_Muon_Wqq']

selection.PlotDistributions(variables)

"""
Questions:
-Should I plot delta phi for all leptons or just leading ones?

"""
