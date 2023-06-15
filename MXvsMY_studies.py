#right now, I am using this script to play around w/ differnt definitions of signal and control regions

import ROOT
from TIMBER.Analyzer import HistGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from XHYbbWW_class import XHYbbWW
from collections import OrderedDict

def KinematicLepton(self): #bringing this function in here so I can select lepton w/o quality/isolation cuts

    self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_phi,Higgs_phi,Wqq_phi)')
    self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_phi,Higgs_phi,Wqq_phi)')
    self.a.Cut('kinLepton_cut','kinEleIdx != -1 || kinMuIdx != -1') #at least one good lepton
    self.a.Define('LeptonType','LeptonIdx(kinEleIdx,kinMuIdx,Electron_pt,Muon_pt)') #picks higher pt signal lepton - output = 0 (lepton is electron) or 1 (lepton is muon)

    self.SIGLEP = self.getNweighted()
    self.AddCutflowColumn(self.SIGLEP,'SIGLEP')

    #For ease, merge some lepton columns that will be useful later (for lepton-type specific variables, use LeptonType to determine if electron or muon)
    self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : Electron_pt[kinEleIdx]')
    self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[kinMuIdx] : Electron_eta[kinEleIdx]')
    self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[kinMuIdx] : Electron_phi[kinEleIdx]')
    self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[kinMuIdx] : Electron_mass[kinEleIdx]')

    return self.a.GetActiveNode()

def MXvsMY_studies(self):

    ##### NEW VARIABLES FOR LATER USE #####

    #W_leptonic transverse mass
    self.a.Define('W_massTran','TransverseMass(MET_pt,Lepton_pt,MET_phi,Lepton_phi)') #Transverse W mass
#   self.a.Define('W_massTran_genMET','TransverseMass(MET_fiducialGenPt,Lepton_pt,MET_fiducialGenPhi,Lepton_phi)') #using generator-level MET variables

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(MET_pt,0,MET_phi,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt,Wqq_eta,Wqq_phi,Wqq_msoftdrop)')
    self.a.Define('Hbb_vect','hardware::TLvector(Higgs_pt,Higgs_eta,Higgs_phi,Higgs_msoftdrop)')

    #Invariant masses of W/Y/X
    self.a.Define('W_massInv','hardware::InvariantMass({MET_vect,Lepton_vect})') #full invariant mass
    self.a.Define('Y_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')
    self.a.Define('X_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect,Hbb_vect})')

    studiesPlots = HistGroup('studiesPlots')

    #######################################
   
    #First lets make some plots examining the lepton quality cuts in the different MC samples
    #Muon_mediumId, Electron_mvaFall17V2noIso vs eta, Electron_mvaFall17V2noIso_WP80, Electron_mvaFall17V2noIso_WP90, Electron_mvaFall17V2noIso_WPL
    
    start=self.a.GetActiveNode()
    muonEvents=self.a.Cut('Muon_events','LeptonType == 1')
    self.a.SetActiveNode(muonEvents)
    self.a.ObjectFromCollection('kinMu','Muon','kinMuIdx')	
   
    #studiesPlots.Add('kinMu_mediumId',self.a.GetActiveNode().DataFrame.Histo1D(('kinMu_mediumId','kinMu_mediumId',2,0,2),'kinMu_mediumId','weight__nominal')) #bins may not work
    
    self.a.SetActiveNode(start)
    electronEvents=self.a.Cut('Electron_events','LeptonType == 0')
    self.a.SetActiveNode(electronEvents)
    self.a.ObjectFromCollection('kinEle','Electron','kinEleIdx')

    #studiesPlots.Add('kinEle_mvaFall17V2noIso vs eta',self.a.DataFrame.Histo2D(('kinEle_mvaFall17V2noIso vs eta','kinEle_mvaFall17V2noIso vs eta',1000,0,1,250,0,2.5),'kinEle_mvaFall17V2noIso', 'kinEle_eta','weight__nominal'))
    
    #Make three plots for electron mva (for different etas/ECAL regions)     
    no_eta = self.a.GetActiveNode()
    inner_barrel = self.a.Cut('inner_barrel','abs(kinEle_eta) < 0.8')
    self.a.SetActiveNode(inner_barrel)
    studiesPlots.Add('kinEle_mvaFall17V2noIso (inner barrel)',self.a.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso (|eta| < 0.8)','kinEle_mvaFall17V2noIso (inner barrel - |eta| < 0.8)',100,0,1),'kinEle_mvaFall17V2noIso', 'weight__nominal'))

    self.a.SetActiveNode(no_eta)
    outer_barrel = self.a.Cut('outer_barrel','abs(kinEle_eta) > 0.8 && abs(kinEle_eta) < 1.479')
    self.a.SetActiveNode(outer_barrel)
    studiesPlots.Add('kinEle_mvaFall17V2noIso (outer barrel)',self.a.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso (0.8 < |eta| < 1.479)','kinEle_mvaFall17V2noIso (outer barrel - 0.8 < |eta| < 1.479)',100,0,1),'kinEle_mvaFall17V2noIso', 'weight__nominal'))

    self.a.SetActiveNode(no_eta)
    endcap = self.a.Cut('endcap','abs(kinEle_eta) > 1.479 && abs(kinEle_eta) < 2.5')
    self.a.SetActiveNode(endcap)
    studiesPlots.Add('kinEle_mvaFall17V2noIso (endcap)',self.a.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso (1.479 < |eta| < 2.5)','kinEle_mvaFall17V2noIso (endcap - 1.479 < |eta| < 2.5)',100,0,1),'kinEle_mvaFall17V2noIso', 'weight__nominal'))

 
    ''' 
    studiesPlots.Add('kinEle_mvaFall17V2noIso_WP80',self.a.GetActiveNode().DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_WP80','kinEle_mvaFall17V2noIso_WP80',2,0,2),'kinEle_mvaFall17V2noIso_WP80','weight__nominal').GetValue())
    print('kinele_mvaWP80 plot made')
    studiesPlots.Add('kinEle_mvaFall17V2noIso_WP90',self.a.GetActiveNode().DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_WP90','kinEle_mvaFall17V2noIso_WP90',2,0,2),'kinEle_mvaFall17V2noIso_WP90','weight__nominal').GetValue())
    print('kinele_mvaWP90 plot made')
    studiesPlots.Add('kinEle_mvaFall17V2noIso_WPL',self.a.GetActiveNode().DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_WPL','kinEle_mvaFall17V2noIso_WPL',2,0,2),'kinEle_mvaFall17V2noIso_WPL','weight__nominal').GetValue())
    print('kinele_mvaWPL plot made')
    '''

    self.a.SetActiveNode(start)

    taggers = ['particleNetMD']

    # now we want to plot mX vs mY for QCD, ttbar, and signal
    for t in taggers:
        self.ApplyMassCuts()
        start=self.a.GetActiveNode()

        # We use Wqq tagging scores to divide data into two regions: signal (enriched in signal) and control (enriched in background)
        #    - Signal:    Wqq > 0.8, pass lepton medium ID
        #    - Control:    Wqq < 0.8, fail lepton medium ID
        # We define a pass/fail criteria for the Hbb score within each region 
        #   - Region 1 (fail):      Hbb < 0.94
        #   - Region 2 (pass):      Hbb > 0.94

        SR=self.ApplySRorCR('SR',t)
        SR_FP=self.ApplyPassFail('SR',t)
            
        self.a.SetActiveNode(start)
        CR=self.ApplySRorCR('CR',t)
        CR_FP=self.ApplyPassFail('CR',t)
  
        nodes=OrderedDict()
        nodes.update(SR_FP)
        nodes.update(CR_FP)

        bins = [80,0,4500] 

        for node in nodes.keys():
            self.a.SetActiveNode(nodes[node])
            print('MX vs MY: Plotting for {}'.format(node))
            studiesPlots.Add('MXvsMY_{}'.format(node), self.a.DataFrame.Histo2D(('MXvsMY_{}'.format(node), 'X vs Y Invariant Mass - {} {}'.format(node.split('_')[1],node.split('_')[0]), bins[0], bins[1], bins[2], bins[0], bins[1], bins[2]), 'X_mass', 'Y_mass', 'weight__nominal'))           
    outFile = ROOT.TFile.Open('{}_{}_{}_MXvsMYstudies.root'.format(self.setname,self.year,self.ijob),'RECREATE')
    outFile.cd()
    studiesPlots.Do('Write') 
    #self.a.PrintNodeTree('NodeTree.pdf',verbose=True)
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

    filename='snapshots/{}_{}_snapshot.txt'.format(setname,year)
    ana = XHYbbWW(filename,ijob,njobs)

#   ana.ApplyStandardCorrections(post_snapshot=True)
    ana.Dijets()
    KinematicLepton(ana)

    MXvsMY_studies(ana)   
