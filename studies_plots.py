#Script for nice-looking plots from JetToLep selection

import ROOT
from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Plot import *
from collections import OrderedDict
from argparse import ArgumentParser
ROOT.gROOT.SetBatch(True)
import sys
sys.path.insert(0,'./modules/')
from pretty_plots import PrettyPlots

if __name__ == '__main__':

    varnames = {
        'MET_pt':'p_{T}^{miss} (GeV)',
        'Neutrino_eta':'#nu #eta',
        'Electron_pt':'e p_{T} (GeV)',
        'Electron_miniPFRelIso_all':'Electron Iso_{mini}',
        'Muon_pt':'#mu p_{T} (GeV)',
        'Muon_miniPFRelIso_all':'Muon Iso_{mini}',
        'Higgs_msoftdrop_corr':'H_{bb} m_{SD} (GeV)',
        'Wqq_msoftdrop_corr':'W_{qq} m_{SD} (GeV)',
        'Higgs_particleNet_mass':'H_{bb} m_{PNet}',
        'Wqq_particleNet_mass':'W_{qq} m_{PNet}',
        'Higgs_particleNetMD_HbbvsQCD':'H_{bb} particleNetMD HbbvsQCD',
        'Wqq_particleNetMD_WqqvsQCD':'W_{qq} particleNetMD WqqvsQCD',
        'DeltaPhi_Higgs_Wqq':'#Delta#phi H_{bb}, W_{qq}',
        'DeltaPhi_Electron_Higgs':'#Delta#phi e, H_{bb}',
        'DeltaPhi_Electron_Wqq':'#Delta#phi e, W_{qq}',
        'DeltaPhi_Muon_Higgs':'#Delta#phi #mu, H_{bb}',
        'DeltaPhi_Muon_Wqq':'#Delta#phi #mu, W_{qq}',
        'DeltaR_Higgs_Wqq':'#DeltaR H_{bb}, W_{qq}',
        'DeltaR_Electron_Higgs':'#DeltaR e, H_{bb}',
        'DeltaR_Electron_Wqq':'#DeltaR e, W_{qq}',
        'DeltaR_Muon_Higgs':'#DeltaR #mu, H_{bb}',
        'DeltaR_Muon_Wqq':'#DeltaR #mu, W_{qq}',
        'MET_pt_nminus1':'p_{T}^{miss} (GeV) (N-1)',
        'kinElectron_pt_nminus1':'e p_{T} (GeV) (N-1)',
        'kinMuon_pt_nminus1':'#mu p_{T} (GeV) (N-1)',
        'kinElectron_miniPFRelIso_all_nminus1':'e Iso_{mini} (N-1)',
        'kinMuon_miniPFRelIso_all_nminus1':'#mu Iso_{mini} (N-1)',
        'Higgs_msoftdrop_corr_nminus1':'H_{bb} m_{SD} (GeV) (N-1)',
        'Wqq_msoftdrop_corr_nminus1':'W_{qq} m_{SD} (GeV) (N-1)',
        'Higgs_particleNet_mass_nminus1':'H_{bb} m_{PNet} (N-1)',
        'Wqq_particleNet_mass_nminus1':'W_{qq} m_{PNet} (N-1)',
        'Higgs_particleNetMD_HbbvsQCD_nminus1':'H_{bb} particleNetMD HbbvsQCD (N-1)',
        'Wqq_particleNetMD_WqqvsQCD_nminus1':'W_{qq} particleNetMD WqqvsQCD (N-1)',
        'DeltaPhi_Higgs_Wqq_allCuts':'#Delta#phi H_{bb}, W_{qq} (all cuts)',
        'DeltaPhi_Lepton_Wqq_allCuts':'#Delta#phi lepton, W_{qq} (all cuts)',
        'DeltaPhi_Lepton_Higgs_allCuts':'#Delta#phi lepton, H_{bb} (all cuts)',
        'DeltaR_Higgs_Wqq_allCuts':'#DeltaR H_{bb}, W_{qq} (all cuts)',
        'DeltaR_Lepton_Higgs_allCuts':'#DeltaR lepton, H_{bb} (all cuts)',
        'DeltaR_Lepton_Wqq_allCuts':'#DeltaR lepton, W_{qq} (all cuts)',
        'Y_mass_allCuts':'Y mass_{inv} (GeV) (all cuts)',
        'X_mass_allCuts':'X mass_{inv} (GeV) (all cuts)',
        'mjj_allCuts':'m_{jj} (GeV) (all cuts)'
    }

    signals_high = OrderedDict([
        ('XHY-4000-2500','m_{X}=4000, m_{Y}=2500'),
        ('XHY-3500-1800','m_{X}=3500, m_{Y}=1800'),
        ('XHY-3000-1400','m_{X}=3000, m_{Y}=1400')
    ])
    
    signals_low = OrderedDict([
       ('XHY-2000-1000','m_{X}=2000, m_{Y}=1000'),
       ('XHY-1600-700','m_{X}=1600, m_{Y}=700'),
       ('XHY-1000-500','m_{X}=1000, m_{Y}=500')
    ])
   
    signals = OrderedDict([
       ('XHY-4000-2500','m_{X}=4000, m_{Y}=2500'),
       ('XHY-2000-1000','m_{X}=2000, m_{Y}=1000'),
       ('XHY-3000-1400','m_{X}=3000, m_{Y}=1400'),
       ('XHY-1000-500','m_{X}=1000, m_{Y}=500')
    ])

    backgrounds = OrderedDict([
        ('ttbar','ttbar'),
        ('ST','single-top'),
        ('WJets','W+Jets')
        #('QCD','QCD'),
        #('VJets','V+Jets'),
        #('DYJets','Drell-Yan'),
        #('VV','Diboson')
        #('WJets','WJets hadronic'),
        #('WJetsLep','WJets semileptonic'),
        #('ZJets','ZJets')
        #('ttbar-allhad','ttbar-allhad'),
        #('ttbar-semilep','ttbar-semilep'),
        #('WW','WW'),
        #('WZ','WZ'),
        #('ZZ','ZZ')
    ])

    colors_high = OrderedDict([
        ('XHY-4000-2500',ROOT.kBlack),
        ('XHY-3500-1800',ROOT.kBlue),
        ('XHY-3000-1400',ROOT.kRed),
        ('ttbar',ROOT.kGreen),
        ('ST',ROOT.kOrange),
        ('WJets',ROOT.kBlue)
        #('WJets',ROOT.kYellow),
        #('WJetsToLNu',ROOT.kOrange),
        #('ZJets',ROOT.kCyan)
        #('ttbar-allhad',ROOT.kMagenta),
        #('ttbar-semilep',ROOT.kYellow),
        #('QCD',ROOT.kGreen),
        #('WJets',ROOT.kMagenta),
        #('ZJets',ROOT.kGreen)
        #('WWto4Q',ROOT.kRed),
        #('WWto1L1Nu2Q',ROOT.kBlue),
        #('WZto1L1Nu2Q',ROOT.kMagenta),
        #('WZto2L2Q',ROOT.kGreen),
        #('ZZto4Q',ROOT.kCyan),
        #('ZZto2L2Q',ROOT.kYellow)
     ])

    colors_low = OrderedDict([
       ('XHY-2000-1000',ROOT.kBlack),
       ('XHY-1600-700',ROOT.kBlue),
       ('XHY-1000-500',ROOT.kRed),
       ('ttbar',ROOT.kGreen),
       ('ST',ROOT.kOrange),
       ('WJets',ROOT.kBlue)
       #('WJetsToLNu',ROOT.kOrange),
       #('ZJets',ROOT.kCyan)
       #('ttbar-allhad',ROOT.kMagenta),
       #('ttbar-semilep',ROOT.kYellow),
       #('QCD',ROOT.kGreen),
       #('WJets',ROOT.kMagenta),
       #('ZJets',ROOT.kGreen)
     ])

    colors = OrderedDict([
       ('XHY-4000-2500',ROOT.kBlack),
       ('XHY-2000-1000',ROOT.kBlue),
       ('XHY-3000-1400',ROOT.kRed),
       ('XHY-1000-500',ROOT.kGray),
       ('ttbar',ROOT.kGreen),
       ('ST',ROOT.kOrange),
       ('WJets',ROOT.kBlue)
     ])

    PrettyPlots('studies',varnames,'high',signals_high,backgrounds,colors_high,cutflow=True)
    PrettyPlots('studies',varnames,'low',signals_low,backgrounds,colors_low,cutflow=True)
    PrettyPlots('studies',varnames,'range',signals,backgrounds,colors)
                                                                                                         

