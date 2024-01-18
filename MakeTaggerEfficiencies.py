'''
Calculates efficiencies of a chosen tagger (particleNetMD Hbb or Wqq) in selecting a given flavor of jet (Wqq, Hbb, Top, or unmerged) at a given working point
Should work on signal and ttbar
'''

import ROOT
from TIMBER.Analyzer import Correction, HistGroup, CutGroup, VarGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME
from collections import OrderedDict
from XHYbbWW_class import XHYbbWW

def make_ratio_hist(histN,histD):
    ratio_hist = histN.Clone('ratio_hist') #New histogram - ratio of (tagged jets)/(all jets) per bin
    #Loop over all bins, manually calculate ratios + set the bin contents of new histogram (I don't think there is a better way to do this for TH2?)
    for x in range(1,ratio_hist.GetNbinsX()+1):
        for y in range(1,ratio_hist.GetNbinsY()+1):
            val_tagged = histN.GetBinContent(x,y)
            val_all = histD.GetBinContent(x,y)
            if val_all == 0:
                ratio = 0
            else:
                ratio = val_tagged/val_all
            ratio_hist.SetBinContent(x,y,ratio)
    return ratio_hist

def MakeEfficiencyMaps(ana,tagger,wp,flavor,jetId):
    ana.nstart=ana.getNweighted()
    ana.AddCutflowColumn(ana.nstart,"nstart")

    #Define subcollections for all target jets and all tagged target jets separately
    ana.a.SubCollection('{}Jets_all'.format(flavor),'FatJet','FatJet_flavor == {}'.format(jetId)) #may need to update the "truth" criteria
    ana.a.SubCollection('{}Jets_tagged'.format(flavor),'FatJet','FatJet_flavor == {} && FatJet_particleNetMD_{}vsQCD > {}'.format(jetId,tagger,wp))
    ana.a.Cut('n{}Jets'.format(flavor),'n{}Jets_all > 0'.format(flavor))
    ana.NJETFLAV = ana.getNweighted()
    ana.AddCutflowColumn(ana.NJETFLAV, "NJETFLAV")

    #Histgroup 
    hgroup = HistGroup('{}Efficiency'.format(flavor))

    #Make histograms for both groups of jets, binned in pt and eta
    hist_Jets_all = ana.a.GetActiveNode().DataFrame.Histo2D(('Jets all','Jets all',60,0,3000,24,-2.4,2.4),'{}Jets_all_pt_corr'.format(flavor),'{}Jets_all_eta'.format(flavor),'weight__nominal')
    hist_Jets_tagged = ana.a.GetActiveNode().DataFrame.Histo2D(('Jets tagged','Jets tagged',60,0,3000,24,-2.4,2.4),'{}Jets_tagged_pt_corr'.format(flavor),'{}Jets_tagged_eta'.format(flavor),'weight__nominal')  
    hgroup.Add(hist_Jets_all.GetTitle(),hist_Jets_all)
    hgroup.Add(hist_Jets_tagged.GetTitle(),hist_Jets_tagged)
 
    ratio_hist = make_ratio_hist(hist_Jets_tagged, hist_Jets_all)
    ratio_hist.SetName('ratio')
    ratio_hist.SetTitle('ratio')
    hgroup.Add(ratio_hist.GetTitle(),ratio_hist)

    #File format is <SETNAME>_<YEAR>_<JET FLAVOR>_<TAGGER FLAVOR>_Efficiencies.root
    outFile = ROOT.TFile.Open('{}_{}_{}Jets_particleNet{}_Efficiencies.root'.format(ana.setname,ana.year,flavor,tagger),'RECREATE')
    outFile.cd()

    #truth matching test
    #hist_Wmass = ana.a.GetActiveNode().DataFrame.Histo1D(('Wmass','Wmass',40,50,250),'WqqJets_all_msoftdrop','weight__nominal')
    #hgroup.Add(hist_Wmass.GetTitle(),hist_Wmass)

    #Save histogram to file    
    hgroup.Do('Write')
    outFile.Close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store', required=True,
                        help='Setname to process.')
    parser.add_argument('-y', type=str, dest='year',
                        action='store', required=True,
                        help='Year of set (16APV, 16, 17, 18).')
    parser.add_argument('-f', type=str, dest='flavor',
                        action='store', required=True,
                        help='Flavor of jets to target for efficiencies: "Hbb", "Wqq", "Top", or "UM" (unmerged)')                                                    
    parser.add_argument('-t', type=str, dest='tagger',
                        help='flavor of jet tagger - "Hbb" or "Wqq"',
                        action='store', required=True)
    parser.add_argument('-w', type=float, dest='wp', 
                        help='tagger working point',
                        action='store', required=True)
    args = parser.parse_args()

    filename = 'snapshots/{}_{}_snapshot.txt'.format(args.setname,args.year)
    ana = XHYbbWW(filename)
    ana.ApplyJMECorrections('None')
    ana.a.MakeWeightCols(extraNominal='' if ana.a.isData else 'genWeight*%s'%ana.GetXsecScale())
 
    jetIds = {'Top':0, 'Wqq': 1, 'Hbb': 2, 'UM': 3}
    jetId = jetIds[args.flavor]

    #Define jetId's for all jets and particleNet taggers
    ana.a.Define('FatJet_flavor','JetFlavor_{}(FatJet_eta, FatJet_phi, GenPart_eta, GenPart_phi, GenPart_pdgId, GenPart_genPartIdxMother)'.format('signal' if 'XHY' in ana.setname else 'ttbar'))
    ana.a.Define('FatJet_particleNetMD_HbbvsQCD','FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD)')
    ana.a.Define('FatJet_particleNetMD_WqqvsQCD','(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc)/(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc+FatJet_particleNetMD_QCD)')

    '''
    hist = ana.a.GetActiveNode().DataFrame.Histo1D(('test','test',4,0,4),'FatJet_flavor','weight__nominal')
    outFile = ROOT.TFile.Open('test.root','RECREATE')
    outFile.cd()
    hist.Write()
    '''

    MakeEfficiencyMaps(ana, args.tagger, args.wp, args.flavor, jetId)
