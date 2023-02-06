import ROOT
from XHYbbWW_class import SplitUp
from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--MX', type=int, dest='mx',
			action='store', help='X mass')
    parser.add_argument('--MY', type=int, dest='my',
			action='store', help='Y mass')
    parser.add_argument('-y', type=str, dest='year',
			action='store', help='year')
    args = parser.parse_args()

    setname = 'XHY-{}-{}_{}'.format(args.mx, args.my, args.year)
    filename = 'raw_nano/Signal/{}.txt'.format(setname)
    #filename='9D389A37-CDE5-7F42-9458-D1DE3405CAEC.root'
    variables = [
	'nElectron',
	'nMuon',
	'Electron_pt',
	'Muon_pt',
	'FatJet_msoftdrop',
	'FatJet_pt',
	'DeltaPhi_jets',
	'DeltaPhi_Electron_Higgs',
	'DeltaPhi_Muon_Higgs',
	'DeltaPhi_Electron_Wqq',
	'DeltaPhi_Muon_Wqq',
	'MET_pt','MET_phi']

    CompileCpp('delta_phi.cc')

    ana = analyzer(filename)

    # first ensure there are no events with 0 jets/electrons/muons so we don't segfault trying to access nonexistent indices
    before = ana.DataFrame.Sum("genWeight")
    ana.Cut('jetsCut','nFatJet > 2')
    ana.Cut('eleCut','nElectron > 2')
    ana.Cut('muCut','nMuon > 2')
    after = ana.DataFrame.Sum("genWeight") 
    # Get the deltaPhi of the 0th and 1st object in every event 
    ana.Define('DeltaPhi_jets','DeltaPhi(FatJet_phi, FatJet_phi, {0,1})')
    ana.Define('DeltaPhi_Electron_Higgs','DeltaPhi(FatJet_phi, Electron_phi, {0,1})')
    ana.Define('DeltaPhi_Muon_Higgs','DeltaPhi(FatJet_phi, Muon_phi, {0,1})')
    ana.Define('DeltaPhi_Electron_Wqq','DeltaPhi(FatJet_phi, Electron_phi, {0,1})')
    ana.Define('DeltaPhi_Muon_Wqq','DeltaPhi(FatJet_phi, Muon_phi, {0,1})')

    # create a histgroup to store all the histograms of variables easily
    histgroup = HistGroup(setname)
    for varname in variables:
	histname = '{}_{}'.format(setname, varname)
	# Arguments for binning that we pass to the RDataFrame::Histo1D() method
	if varname.startswith('n'):
	    hist_tuple = (histname,histname,10,0,10)
	elif 'pt' in varname:
	    hist_tuple = (histname,histname,40,0,2000)
	elif 'msoftdrop' in varname:
	    hist_tuple = (histname,histname,30,0,3000)
	elif 'DeltaPhi' in varname or varname=='MET_phi':
	    hist_tuple = (histname,histname,30,-3.2,3.2)
        else:
            hist_tuple = (histname,histname,30,40,200)
	print(hist_tuple)
	print(varname)
	hist = ana.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname)
	hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1
	histgroup.Add(varname,hist)

    # save the raw histos to a file
    outFile = ROOT.TFile.Open('{}.root'.format(setname),'RECREATE')
    outFile.cd()
    histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group
    outFile.Close()

    # see how many events were lost with the Jet/electron/muon > 2 cut
    print('electron/muon/Jet > 2 cut efficiency = {}%'.format(after.GetValue()/before.GetValue()*100.))

    '''
    # now make the nice plots
    for varname in variables:
	fName = '{}_{}.png'.format(setname,varname)
    '''	
