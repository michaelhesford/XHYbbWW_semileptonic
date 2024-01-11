import ROOT
from TIMBER.Analyzer import Correction, HistGroup, CutGroup, VarGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME
from collections import OrderedDict
from XHYbbWW_class import XHYbbWW

def MakeHbbEfficiencyMaps(ana,wp):
    ana.a.Define('FatJet_particleNetMD_HbbvsQCD','FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD)') #Define tagger
    #Define subcollections for all Hbb jets and all tagged Hbb jets separately
    ana.a.SubCollection('HbbJets_all','FatJet','FatJet_hadronFlavour == 5 && FatJet_nBHadrons == 2') #may need to update the "truth" criteria
    ana.a.SubCollection('HbbJets_tagged','FatJet','FatJet_hadronFlavour == 5 && FatJet_nBHadrons == 2 && FatJet_particleNetMD_HbbvsQCD > {}'.format(wp))
    ana.a.Cut('nHbbJets','nHbbJets_all > 0')

    #Histgroup 
    hgroup = HistGroup('HbbEfficiency')

    #Make histograms for both groups of jets, binned in pt and eta
    hist_HbbJets_all = ana.a.GetActiveNode().DataFrame.Histo2D(('HbbJets all','HbbJets all',60,0,3000,24,-2.4,2.4),'HbbJets_all_pt_corr','HbbJets_all_eta','weight__nominal')
    hist_HbbJets_tagged = ana.a.GetActiveNode().DataFrame.Histo2D(('HbbJets tagged','HbbJets tagged',60,0,3000,24,-2.4,2.4),'HbbJets_tagged_pt_corr','HbbJets_tagged_eta','weight__nominal')  
    hgroup.Add(hist_HbbJets_all.GetTitle(),hist_HbbJets_all)
    hgroup.Add(hist_HbbJets_tagged.GetTitle(),hist_HbbJets_tagged)

    ratio_hist = hist_HbbJets_tagged.Clone('ratio_hist') #New histogram - ratio of (tagged Hbb jets)/(all Hbb jets) per bin
    
    #Loop over all bins, manually calculate ratios + set the bin contents of new histogram (I don't think there is a better way to do this for TH2?)
    for x in range(1,ratio_hist.GetNbinsX()+1):
        for y in range(1,ratio_hist.GetNbinsY()+1):
            val_Hbb_tagged = hist_HbbJets_tagged.GetBinContent(x,y)
            val_Hbb_all = hist_HbbJets_all.GetBinContent(x,y)
            if val_Hbb_all == 0:
                ratio = 0
            else:
                ratio = val_Hbb_tagged/val_Hbb_all
            ratio_hist.SetBinContent(x,y,ratio)    

    ratio_hist.SetName('HbbJets ratio')
    ratio_hist.SetTitle('HbbJets ratio')
    hgroup.Add(ratio_hist.GetTitle(),ratio_hist)

    outFile = ROOT.TFile.Open('{}_{}_HbbEfficiencies.root'.format(ana.setname,ana.year),'RECREATE')
    outFile.cd()

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
                        help='Year of set (16, 17, 18).')                                                          
    args = parser.parse_args()

    filename = 'snapshots/{}_{}_snapshot.txt'.format(args.setname,args.year)
    ana = XHYbbWW(filename)
    ana.ApplyJMECorrections('None')
    ana.a.MakeWeightCols(extraNominal='' if ana.a.isData else 'genWeight*%s'%ana.GetXsecScale())
    
    MakeHbbEfficiencyMaps(ana,0.98)
    
