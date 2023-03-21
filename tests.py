import ROOT
from XHYbbWW_class import XHYbbWW
from TIMBER.Analyzer import HistGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
from collections import OrderedDict
from argparse import ArgumentParser

def TestPlotter(setname,year,ijob,njobs,variables):

    if 'XHY' in setname:
        signal=True
        filename = 'raw_nano/Signal/{}_{}.txt'.format(setname,year)
    else:
        signal=False
        if 'QCD' in setname:
            filename = 'raw_nano/Background/QCD/{}_{}.txt'.format(setname,year)
        elif 'ttbar' in setname:
            filename = 'raw_nano/Background/ttbar/{}_{}.txt'.format(setname,year)
        elif 'Jets' in setname:
            filename = 'raw_nano/Background/V+Jets/{}_{}.txt'.format(setname,year)

    CompileCpp('HWWmodules.cc')
    CompileCpp('delta_phi.cc')
 
    ana = XHYbbWW(filename,ijob,njobs)

    ana.NSTART = ana.getNweighted()
    ana.AddCutflowColumn(ana.NSTART,'NSTART')

    #first ensure enough leptons/jets exist to avoid referencing indexes which don't exist
    ana.a.Cut('nLepton','nElectron > 0 || nMuon > 0') #NOTE: THIS WILL CAUSE PROBLEMS WITH BASIC DELTA PHI'S
    ana.NLEPTON = ana.getNweighted()
    ana.AddCutflowColumn(ana.NLEPTON,'NLEPTON')

    start = ana.a.Cut('nFatJet','nFatJet > 1') #node to use for later
    ana.NJET = ana.getNweighted()
    ana.AddCutflowColumn(ana.NJET,'NJET')

    #Here is the section which I am mainly editing

    #ana.a.SubCollection('Lead','FatJet',{0},useTake=True) #GOTTA FIGURE OUT WHAT'S WRONG WITH THIS, BUT FOR NOW JUST MANUALLY DEFINE THINGS
    #ana.a.SubCollection('Sublead','FatJet',{1},useTake=True)

    ana.a.Define('Lead_pt','FatJet_pt[0]')
    ana.a.Define('Lead_msoftdrop','FatJet_msoftdrop[0]')
    ana.a.Define('Sublead_pt','FatJet_pt[1]')
    ana.a.Define('Sublead_msoftdrop','FatJet_msoftdrop[1]')    

    ana.SignalLepton()
    ana.Dijets() 

    #ana.a.SubCollection('Wqq','Dijet',{0},useTake=True)
    #ana.a.SubCollection('Higgs','Dijet',{1},useTake=True)

    ana.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_particleNetMD_Xbb/(Dijet_particleNetMD_Xbb+Dijet_particleNetMD_QCD)') #mass-decorrelated H tagger
    ana.a.Define('Dijet_particleNetMD_WqqvsQCD','Dijet_particleNetMD_Xqq/(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_QCD)')

    ana.a.Define('Wqq_msoftdrop','Dijet_msoftdrop[0]')
    ana.a.Define('Wqq_particleNetMD_WqqvsQCD','Dijet_particleNetMD_WqqvsQCD[0]')
    ana.a.Define('Wqq_particleNetMD_HbbvsQCD','Dijet_particleNetMD_HbbvsQCD[0]')
    ana.a.Define('Higgs_msoftdrop','Dijet_msoftdrop[1]')
    ana.a.Define('Higgs_particleNetMD_WqqvsQCD','Dijet_particleNetMD_WqqvsQCD[1]')
    ana.a.Define('Higgs_particleNetMD_HbbvsQCD','Dijet_particleNetMD_HbbvsQCD[1]')

    ana.a.Define('DeltaPhi_Lepton_Wqq','DeltaPhiSingular(Dijet_phi,Lepton_phi,0)')
    dijets = ana.a.Define('DeltaPhi_Lepton_Higgs','DeltaPhiSingular(Dijet_phi,Lepton_phi,1)') #set this as a node

    #ana.JetSelection(mH = [100,150], mW = [60,100], Htag = 0.8, Wtag = 0.8)

    histgroup = HistGroup(setname)
    for varname in variables:
        if 'msoftdrop' in varname:
            binning = [20,0,200]
        elif '_pt' in varname:
            binning = [40,0,2000]
        elif 'particleNet' in varname:
            binning = [10,0,1]
        elif 'DeltaPhi' in varname:
            binning = [30,-3.2,3.2]
        elif 'nEle' or 'nMuo' or 'nFat' in varname:
            binning = [10,0,10]

        if 'Higgs' or 'Wqq' in varname:
            ana.a.SetActiveNode(dijets)
        else:
            ana.a.SetActiveNode(start)

        histname = '{}_{}'.format(setname,varname)
        hist_tuple = (histname,histname,binning[0],binning[1],binning[2])
        print('creating {} histogram'.format(histname))
        hist = ana.a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname)
        hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1
        histgroup.Add(varname,hist)

    # save the raw histos to a file
    outFile = ROOT.TFile.Open('plots/tests/{}_{}_DijetTests.root'.format(setname,year),'RECREATE')
    outFile.cd()
    histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group

    #create cut flow table
    cutflowInfo = OrderedDict([
        ('NSTART',ana.NSTART),
        ('NLEPTON',ana.NLEPTON),
        ('NJET',ana.NJET),
        ('SIGLEP',ana.SIGLEP),
        ('NDIJETS',ana.NDIJETS),
    ])

#        ('JETSELECTION',ana.JETSELECTION)

    nLabels = len(cutflowInfo)
    hCutflow = ROOT.TH1F('cutflow', 'Number of events after each cut', nLabels, 0.5, nLabels+0.5)
    nBin = 1
    for label, value in cutflowInfo.items():
        hCutflow.GetXaxis().SetBinLabel(nBin, label)
        hCutflow.AddBinContent(nBin, value)
        nBin += 1
    hCutflow.Write()

    outFile.Close()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store',help='name of data set to run on')
    parser.add_argument('-y', type=str, dest='year',
                        action='store', help='year',required=False)
    parser.add_argument('-j', type=int, dest='ijob',required=False,
                        action='store', help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
                        action='store', help='number of jobs')
    args = parser.parse_args()

    setname=args.setname
    year=args.year
    ijob=args.ijob
    njobs=args.njobs

    variables = ['DeltaPhi_Lepton_Higgs',
                 'DeltaPhi_Lepton_Wqq',
                 'Lead_pt',
                 'Lead_msoftdrop',
                 'Sublead_pt',
                 'Sublead_msoftdrop',
                 'Higgs_msoftdrop',
                 'Higgs_particleNetMD_WqqvsQCD',
                 'Higgs_particleNetMD_HbbvsQCD',
                 'Wqq_msoftdrop',
                 'Wqq_particleNetMD_WqqvsQCD',
                 'Wqq_particleNetMD_HbbvsQCD',
                 'nElectron',
                 'nMuon',
                 'nFatJet'
    ]
    

    TestPlotter(setname,year,ijob,njobs,variables)









