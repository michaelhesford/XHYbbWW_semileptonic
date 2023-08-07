#Final selection script, outputs templates for background estimation

import ROOT, time
from TIMBER.Analyzer import HistGroup, Correction
from TIMBER.Tools.Common import CompileCpp
from collections import OrderedDict
import TIMBER.Tools.AutoJME as AutoJME
from XHYbbWW_class import XHYbbWW

def adjust_negative_bins(template_hists):
    #currently have some negative bins in some VV MC histos - temporary fix is just to zero them out, or (for better systematics) set them to a slightly positive value
    for temp in template_hists.keys():
        hist = template_hists[temp]
        print('Adjusting negative bins for histogram {}'.format(temp))
        for x in range(0,hist.GetNbinsX()):
            for y in range(0,hist.GetNbinsY()):
                if hist.GetBinContent(x,y) < 0.0:
                    hist.SetBinContent(x,y,0.002)

def loosen_config_cuts(self):
    #values can be adjusted manually here
    self.cuts['JET']['particleNetMD_HbbvsQCD'] = 0.6
    self.cuts['JET']['particleNetMD_WqqvsQCD'] = [0,0.8]
    self.cuts['JET']['btag'] = 0.6
    self.cuts['JET']['mh'] = [60,150]
    self.cuts['JET']['mw'] = [50,130]
    self.cuts['ELECTRON']['RelIso'] = [0.5,1]
    self.cuts['ELECTRON']['mva80'] = 'Electron_mvaFall17V2noIso_WPL'
    self.cuts['MUON']['RelIso'] = [0.5,1]
    self.cuts['MUON']['mediumId'] = 'Muon_looseId'

def rescale_hists(self,region,binsX,binsY,node,template_hists):
    stdInts = self.get_standard_int(region,binsX,binsY,node)
    for temp in template_hists.keys():
        hist = template_hists[temp]
        looseInt = hist.Integral()
        ratio = stdInts[temp]/looseInt
        hist.Scale(ratio)

def XHYbbWW_selection(self,variation,nom_hists):

    ##### NEW VARIABLES FOR LATER USE #####

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(MET_pt,0,MET_phi,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)') #use updated mass/pt after JME corrections
    self.a.Define('Hbb_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)')

    #Invariant masses of "Y"/"X"
    self.a.Define('mwlv','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')
    self.a.Define('mhwlv','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect,Hbb_vect})')

    taggers = ['particleNetMD']

    #######################################

    if not nom_hists
    outFile = ROOT.TFile.Open('XHYbbWWselection_{}_{}.root'.format(self.setname, variation if variation != 'None' else ''),'UPDATE')
    outFile.cd()

    # now we want to plot mX vs mY for QCD, ttbar, and signal
    for t in taggers:
        pre_mcut=self.a.GetActiveNode()

        if 'to' in self.setname:
            loosen_config_cuts(self)
    
        self.ApplyMassCuts()
        post_mcut=self.a.GetActiveNode()

        # We use Wqq tagging scores and lepton cuts to divide data into two regions: signal (enriched in signal) and control (enriched in background) - our control region is designed yield plenty of ttbar, which is the dominant background in this search
        #    - Signal:    Wqq > 0.8, miniRelIso < 0.2, pass lepton medium ID
        #    - Control:   Wqq < 0.8, miniRelIso < 0.2, pass lepton medium ID, b-tagged AK4 jet exists
        # We define a pass/fail criteria for the Hbb score within each region 
        #   - Region 1 (fail):      Hbb < 0.94
        #   - Region 2 (pass):      Hbb > 0.94

        #signal region
        SR=self.ApplySRorCR('SR',t)
        SR_FP=self.ApplyPassFail('SR',t)

        #control region - ttbar enriched
        self.a.SetActiveNode(post_mcut)
        ttCR=self.ApplySRorCR('ttCR',t)
        ttCR_FP=self.ApplyPassFail('ttCR',t)

        nodes=OrderedDict()
        nodes.update(SR_FP)
        nodes.update(ttCR_FP)

        binsX = [90,0,4500]
        binsY = [90,0,4500]      

        for region in nodes.keys():
            self.a.SetActiveNode(nodes[region])
            print('MX vs MY: Evaluating for {}'.format(region))
            if nom_hists:
                #just make nominal histograms for each region to store for later use
                self.a.GetActiveNode().DataFrame.Histo2D(('MXvMY_{}_{}__nominal'.format(region,self.year,'nominal'),'X vs Y Invariant Mass - {} {}'.format(region.split('_')[1],region.split('_')[0]),binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),'mhwlv','mwlv','weight__nominal') 

            templates = selection.a.MakeTemplateHistos(ROOT.TH2F('MXvMY_{}_{}'.format(region,year), 'X vs Y Invariant Mass - {} {}'.format(region.split('_')[1],region.split('_')[0]), binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),['mhwlv','mwlv'])
            if 'to' in self.setname:
                rescale_hists(self,region,binsX,binsY,pre_mcut,templates)
            adjust_negative_bins(templates)
            templates.Do('Write')

    cutflowInfo = OrderedDict([
       ('nDijets',selection.nDijets),
       ('nkinLep',selection.nkinLep),
       ('nHiggs',selection.nHiggs),
       ('nWqq',selection.nWqq),
       ('nWtag_SR',selection.nWtag_SR),
       ('nlepIso_SR',selection.nlepIso_SR),
       ('nlepQ_SR',selection.nlepQ_SR),
       ('nWtag_ttCR',selection.nWtag_ttCR),
       ('nlepIso_ttCR',selection.nlepIso_ttCR),
       ('nlepQ_ttCR',selection.nlepQ_ttCR),
       ('nJetB_ttCR',selection.nJetB_ttCR)
    ])

    nLabels = len(cutflowInfo)
    hCutflow = ROOT.TH1F('cutflow', 'Number of events after each cut', nLabels, 0.5, nLabels+0.5)
    nBin = 1
    for label, value in cutflowInfo.items():
	hCutflow.GetXaxis().SetBinLabel(nBin, label)
	hCutflow.AddBinContent(nBin, value)
	nBin += 1
    hCutflow.Write()

    if not selection.a.isData:
        scale = ROOT.TH1F('scale','xsec*lumi/genEventSumw',1,0,1)
        scale.SetBinContent(1,selection.GetXsecScale())
        scale.Write()

    #self.a.PrintNodeTree('NodeTree.pdf',verbose=True)
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
    parser.add_argument('-v', type=str, dest='variation',
                        action='store', default='None',
                        help='JES_up, JES_down, JMR_up,...')
    parser.add_argument('--make_nominal_hists',action='store_true',help='if used, just make nominal histograms')
    args = parser.parse_args()
    nom_hists = args.make_nominal_hists
    setname=args.setname
    year=args.year
    variation=args.variation if not nom_hists else 'None'

    if 'to' in setname:
        #using a looser snapshot to increase statistics for VV samples
        filename='VVsnapshots_loose/{}_{}_VVsnapshot_loose.txt'.format(setname,year)    
    else:
        filename='snapshots/{}_{}_snapshot.txt'.format(setname,year)

    selection = XHYbbWW(filename,1,1)
    selection.ApplyStandardCorrections(snapshot=False)
    selection.ApplyJMECorrections(variation)
    selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale())

    selection.Dijets()
    selection.KinematicLepton()

    XHYbbWW_selection(selection,variation,nom_hists)
   
