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
    self.cuts['JET']['btag'] = 0.1
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

def XHYbbWW_selection(self,variation):

    ##### NEW VARIABLES FOR LATER USE #####

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(ChsMET_pt,0,ChsMET_phi,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)') #use updated mass/pt after JME corrections
    self.a.Define('Hbb_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)')

    #Invariant masses of "Y"/"X"
    self.a.Define('mwlv','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')
    self.a.Define('mhwlv','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect,Hbb_vect})')

    taggers = ['particleNetMD']

    #######################################

    outFile = ROOT.TFile.Open('XHYbbWWselection_{}_{}{}.root'.format(self.setname, self.year, '_'+variation if variation != 'None' else ''),'UPDATE')
    outFile.cd()

    # now we want to plot mX vs mY for QCD, ttbar, and signal
    for t in taggers:

        if 'to' in self.setname:
            loosen_config_cuts(self)

        start=self.a.GetActiveNode()
        self.ApplyMassCuts()
        post_mcut=self.a.GetActiveNode()

        #signal region
        SR=self.ApplySRorCR('SR',t)
        SR_FP=self.ApplyPassFail('SR',t)
       
        #control region - ttbar enriched
        self.a.SetActiveNode(post_mcut)
        ttCR=self.ApplySRorCR('ttCR_HF',t)
        #ttCR=self.ApplySRorCR('ttCR_dPhi',t)
        ttCR_FP=self.ApplyPassFail('ttCR',t)

        nodes=OrderedDict()
        nodes.update(SR_FP)
        nodes.update(ttCR_FP)

        binsX = [90,0,4500]
        binsY = [90,0,4500]      

        for region in nodes.keys():
            self.a.SetActiveNode(nodes[region])
            print('MX vs MY: Evaluating for {}'.format(region))
            templates = HistGroup('MXvMY')
            plot_vars = ['mhwlv','mwlv']
            templates = selection.a.MakeTemplateHistos(ROOT.TH2F('MXvMY_{}'.format(region), 'X vs Y Invariant Mass - {} {}'.format(region.split('_')[1],region.split('_')[0]), binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),plot_vars)
             
            if 'to' in self.setname:
                rescale_hists(self,region,binsX,binsY,pre_mcut,templates) #should find a way to do this post combining histos
            #adjust_negative_bins(templates) #Sholud be saved for after histograms are combined for all years
            templates.Do('Write')
           
    cutflowInfo = OrderedDict([
       ('start',self.nstart),
       ('nTrigs',self.nTrigs),
       ('nkinLep',self.nkinLep),
       ('nDijets',self.nDijets),
       ('nHiggs',self.nHiggs),
       ('nWqq',self.nWqq),
       ('nHtag_SR',self.nHtag_SR),
       ('nDPhiLH_SR',self.nDPhiLH_SR),
       ('nWP_SR',self.nWP_SR),
       ('nWF_SR',self.nWF_SR),
       ('nHtag_ttCR',self.nHtag_ttCR),
       #('nDPhiLH_ttCR',self.nDPhiLH_ttCR),
       ('nJetB_ttCR_HF',self.nJetB_ttCR_HF),
       #('nJetB_ttCR_dPhi',self.nJetB_ttCR_dPhi),
       ('nWP_ttCR',self.nWP_ttCR),
       ('nWF_ttCR',self.nWF_ttCR)
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
        scale = ROOT.TH1F('scale','genEventSumw',1,0,1)
        scale.SetBinContent(1,selection.a.genEventSumw)
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
    parser.add_argument('-j', type=int, dest='ijob',required=False,
        action='store', default=1, help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
        action='store', default=1, help='number of jobs')
    args = parser.parse_args()
    setname=args.setname
    year=args.year
    ijob=args.ijob
    njobs=args.njobs
    variation = args.variation

    if 'to' in setname:
        #using a looser snapshot to increase statistics for VV samples
        filename='VVsnapshots_loose/{}_{}_VVsnapshot_loose.txt'.format(setname,year)    
    else:
        filename='snapshots/{}_{}_snapshot.txt'.format(setname,year)

    selection = XHYbbWW(filename,ijob,njobs)
    selection.nstart = selection.getNweighted()
    selection.AddCutflowColumn(selection.nstart,'nstart')

    selection.ApplyTrigs()
    selection.ApplyJMECorrections(variation)
    selection.Dijets()
    selection.KinematicLepton()
    selection.ApplyStandardCorrections(snapshot=False)
    selection.a.MakeWeightCols(extraNominal='' if selection.a.isData else 'genWeight*%s'%selection.GetXsecScale())

    XHYbbWW_selection(selection,variation)
   
