import ROOT
from XHYbbWW_class import SplitUp
from TIMBER.Analyzer import analyzer
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser

parser = ArgumentParser()

'''
parser.add_argument('-s', type=str, dest='setname',
                    action='store', required=True,
                    help='Setname to process.')
parser.add_argument('-y', type=int, dest='year',
                    action='store', required=True)
'''
parser.add_argument('-j', type=int, dest='ijob',
                    action='store', default=1,
                    help='Job number')
parser.add_argument('-n', type=int, dest='njobs',
                    action='store', default=1,
                    help='Number of jobs')

args = parser.parse_args()

MxMy=[['300','125'],['500','300'],['1000','400'],['3000','500'],['3500','500'],['4000','1000']]
setnames = ['XHY-'+i[0]+'-'+i[1] for i in MxMy]

variables=['nElectron','nMuon','Electron_pt[0]','Muon_pt[0]','FatJet_msoftdrop[0]','FatJet_msoftdrop[1]','FatJet_pt[0]','FatJet_pt[1]','DeltaPhi_jets','DeltaPhi_Electron_Higgs','DeltaPhi_Muon_Higgs','DeltaPhi_Electron_Wqq','DeltaPhi_Muon_Wqq','MET_pt','MET_phi']

variables_snap=['nElectron','nMuon','Electron_pt','Muon_pt','FatJet_msoftdrop','FatJet_pt','DeltaPhi_jets','DeltaPhi_Electron_Higgs','DeltaPhi_Muon_Higgs','DeltaPhi_Electron_Wqq','DeltaPhi_Muon_Wqq','MET_pt','MET_phi']

var_dict={k:None for k in variables}
for var in variables:
    if 'nElectron' or 'nMuon' in var:
        var_dict[var] = [10,0,10]
    elif 'msoftdrop' in var:
        var_dict[var] = [30,40,200]
    elif 'Electron_pt' or 'Muon_pt' in var:
        var_dict[var] = [30,20,200]
    elif 'DeltaPhi' in var:
        var_dict[var] = [8,0,4]
    else:
        var_dict[var] = [30,200,2000] #Fat jet pt

var_names=var_dict.keys()
binning=var_dict.values()

CompileCpp('delta_phi.cc')

for i in range(len(var_names)):
    var = var_names[i]
    exec("c{} = ROOT.TCanvas('{}','{}')".format(i,var,var))
    exec('c{}.Divide(2,3)'.format(i))

for i in range(len(setnames)):
    name = setnames[i]
    filename = 'raw_nano/Signal/{}_18.txt'.format(name)

    infiles = SplitUp(filename, args.njobs)[args.ijob-1]
    ana = analyzer(infiles)

    ana.Define('DeltaPhi_jets','DeltaPhi(FatJet_phi[0],FatJet_phi[1])')
    ana.Define('DeltaPhi_Electron_Higgs','DeltaPhi(FatJet_phi[0],Electron_phi[0])')
    ana.Define('DeltaPhi_Muon_Higgs','DeltaPhi(FatJet_phi[0],Muon_phi[0])')
    ana.Define('DeltaPhi_Electron_Wqq','DeltaPhi(FatJet_phi[1],Electron_phi[0])')
    ana.Define('DeltaPhi_Muon_Wqq','DeltaPhi(FatJet_phi[1],Muon_phi[1])')

    ana.Snapshot(variables_snap,'PlotData_{}.root'.format(name),'Events', openOption='RECREATE',saveRunChain=True)

    infile=ROOT.TFile.Open('PlotData_{}.root'.format(name),'READ')
    tree=infile.Get('Events')

    for j in range(len(var_names)):
        var = var_names[j]
        bins = binning[j]

        exec('c{}.cd({})'.format(j,i+1))
 
        hist=ROOT.TH1D(name,name,bins[0],bins[1],bins[2])

        nEntries=tree.GetEntries()
        for event in range(nEntries):
            tree.GetEntry(event)
            exec('entry = tree.{}'.format(var))
            hist.Fill(entry)

        hist.Draw()

    infile.Close()

for i in range(len(var_names)):
    exec('c{}.Print(distributions/{}.pdf'.format(i,var_names[i]))
    exec('c{}.Clear()'.format(i))




'''
        hist_tuple = ('{}'.format(name),'{}'.format(name),bins[0],bins[1],bins[2])
        hist = ana.DataFrame.Histo1D(hist_tuple,'{}'.format(var))
        hist.GetValue()
        hist.Draw()
'''

'''
for i in range(len(var_names)):

    var = var_names[i]
    bins = binning[i]

    c = ROOT.TCanvas(var,var)
    c.Divide(2,3)

    for j in range(len(setnames)):
        name = setnames[j]
        filename = 'raw_nano/Signal/{}_18.txt'.format(name)

        infiles = SplitUp(filename, args.njobs)[args.ijob-1]
        ana = analyzer(infiles)

        ana.Define('DeltaPhi_jets','DeltaPhi(FatJet_phi[0],FatJet_phi[1])')
        ana.Define('DeltaPhi_Electron_Higgs','DeltaPhi(FatJet_phi[0],Electron_phi[0])')
        ana.Define('DeltaPhi_Muon_Higgs','DeltaPhi(FatJet_phi[0],Muon_phi[0])')
        ana.Define('DeltaPhi_Electron_Wqq','DeltaPhi(FatJet_phi[1],Electron_phi[0])')
        ana.Define('DeltaPhi_Muon_Wqq','DeltaPhi(FatJet_phi[1],Muon_phi[1])')

        c.cd(j+1)

        hist_tuple = ('{}'.format(name),'{}'.format(name),bins[0],bins[1],bins[2])
        hist = ana.DataFrame.Histo1D(hist_tuple,'{}'.format(var))
        hist.GetValue()
        hist.Draw()

    c.Print('distributions/{}.pdf'.format(var))
    c.Clear()
'''
