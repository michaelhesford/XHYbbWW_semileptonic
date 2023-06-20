'''
Script w/ functions to produce basic plots of MC along with N-1 plots - consolidated from previous work. 
'''

import ROOT
from TIMBER.Analyzer import HistGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from XHYbbWW_class import XHYbbWW
from collections import OrderedDict

def KinPlots(self):

    ### DEFINITIONS ###

    #For starters, we will take the two highest pt jets to be our objects of interest
    #self.a.Define('DijetIdxs','ROOT::VecOps::RVec({0,1})')
    #self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True))

    self.a.ObjectFromCollection('Lead','FatJet','0')
    self.a.ObjectFromCollection('Sublead','FatJet','1')

    #Define mass-decorrelated Higgs/Wqq taggers
    self.a.Define('FatJet_particleNetMD_HbbvsQCD','FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD)')
    self.a.Define('FatJet_particleNetMD_WqqvsQCD','FatJet_particleNetMD_Xqq/(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_QCD)')

    #mark one jet as Higgs and other as Wqq
    self.a.Define('HiggsIdx','FatJet_particleNetMD_HbbvsQCD[0] > FatJet_particleNetMD_HbbvsQCD[1] ? 0 : 1')
    self.a.Define('WqqIdx','HiggsIdx == 0 ? 1 : 0')

    self.a.ObjectFromCollection('Higgs','FatJet','HiggsIdx')
    self.a.ObjectFromCollection('Wqq','FatJet','WqqIdx')

    #We will also take only the highest pt lepton
    self.a.Define('LeptonType','leptonType_studies(Electron_pt, Muon_pt)') #output: 0 - electron event OR 1 - muon event   
    self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[0] : Electron_pt[0]')
    self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[0] : Electron_eta[0]')
    self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[0] : Electron_phi[0]')
    self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[0] : Electron_mass[0]')
    self.a.Define('Lepton_miniPFRelIso_all','LeptonType == 1 ? Muon_miniPFRelIso_all[0] : Electron_miniPFRelIso_all[0]')

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(MET_pt,0,MET_phi,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt,Wqq_eta,Wqq_phi,Wqq_msoftdrop)')
    self.a.Define('Higgs_vect','hardware::TLvector(Higgs_pt,Higgs_eta,Higgs_phi,Higgs_msoftdrop)')

    #Define some delta phi variables
    self.a.Define('DeltaPhi_Higgs_Wqq','hardware::DeltaPhi(Higgs_phi,Wqq_phi)')    
    self.a.Define('DeltaPhi_Lepton_MET','hardware::DeltaPhi(Lepton_phi,MET_phi)')
    self.a.Define('DeltaPhi_Lepton_Higgs','hardware::DeltaPhi(Lepton_phi,Higgs_phi)')
    self.a.Define('DeltaPhi_Lepton_Wqq','hardware::DeltaPhi(Lepton_phi,Wqq_phi)')

    #Now look at Delta R
    self.a.Define('DeltaR_Higgs_Wqq','hardware::DeltaR(Higgs_vect,Wqq_vect)')
    self.a.Define('DeltaR_Lepton_Higgs','hardware::DeltaR(Lepton_vect,Higgs_vect)')
    self.a.Define('DeltaR_Lepton_Wqq','hardware::DeltaR(Lepton_vect,Wqq_vect)')

    ###################

    start = self.a.GetActiveNode()
    kinPlots=HistGroup('kinPlots')

    kinPlots.Add('MET_pt',self.a.GetActiveNode().DataFrame.Histo1D(('MET_pt','MET_pt',40,0,1000),'MET_pt','weight__nominal'))
    kinPlots.Add('Lead_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Lead_pt','Lead_pt',34,300,2000),'Lead_pt','weight__nominal'))
    kinPlots.Add('Sublead_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Sublead_pt','Sublead_pt',36,200,2000),'Sublead_pt','weight__nominal'))
    kinPlots.Add('DeltaPhi_Higgs_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Higgs_Wqq','DeltaPhi_Higgs_Wqq',30,-3.2,3.2),'DeltaPhi_Higgs_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Higgs_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Higgs_Wqq','DeltaR_Higgs_Wqq',25,0,5),'DeltaR_Higgs_Wqq','weight__nominal'))
    kinPlots.Add('Higgs_msoftdrop',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_msoftdrop','Higgs_msoftdrop',40,50,250),'Higgs_msoftdrop','weight__nominal'))
    kinPlots.Add('Wqq_msoftdrop',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_msoftdrop','Wqq_msoftdrop',40,50,250),'Wqq_msoftdrop','weight__nominal'))

    #Make separate kinematic plots for electron/muon events
    eleEvents=self.a.Cut('eleEvents','LeptonType == 0')
    self.a.SetActiveNode(eleEvents)

    kinPlots.Add('Electron_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_pt','Electron_pt',50,10,1000),'Lepton_pt','weight__nominal'))
    kinPlots.Add('DeltaPhi_Electron_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Electron_Higgs','DeltaPhi_Electron_Higgs',30,-3.2,3.2),'DeltaPhi_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaPhi_Electron_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Electron_Wqq','DeltaPhi_Electron_Wqq',30,-3.2,3.2),'DeltaPhi_Lepton_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Electron_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Electron_Higgs','DeltaR_Electron_Higgs',25,0,5),'DeltaR_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaR_Electron_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Electron_Wqq','DeltaR_Electron_Wqq',25,0,5),'DeltaR_Lepton_Wqq','weight__nominal'))

    self.a.SetActiveNode(start)
    muEvents=self.a.Cut('muEvents','LeptonType == 1')
    self.a.SetActiveNode(muEvents)

    kinPlots.Add('Muon_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_pt','Muon_pt',50,10,1000),'Lepton_pt','weight__nominal')) 
    kinPlots.Add('DeltaPhi_Muon_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Muon_Higgs','DeltaPhi_Muon_Higgs',30,-3.2,3.2),'DeltaPhi_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaPhi_Muon_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Muon_Wqq','DeltaPhi_Muon_Wqq',30,-3.2,3.2),'DeltaPhi_Lepton_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Muon_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Muon_Higgs','DeltaR_Muon_Higgs',25,0,5),'DeltaR_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaR_Muon_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Muon_Wqq','DeltaR_Muon_Wqq',25,0,5),'DeltaR_Lepton_Wqq','weight__nominal'))
    
    self.a.SetActiveNode(start)

    return kinPlots

def Nminus1Group(t):
    NCuts=CutGroup('NCuts')
    #Cuts to apply for all plots
    NCuts.Add('Lepton_pt','Lepton_pt > 25')
    NCuts.Add('Lepton_eta','LeptonType == 0? Lepton_eta < 2.5 : Lepton_eta < 2.4')
    NCuts.Add('Dijet_pt','Lead_pt > 300 && Sublead_pt > 200')
    NCuts.Add('Dijet_eta','abs(Higgs_eta) < 2.4 && abs(Wqq_eta) < 2.4')
    #Now for the cuts we actually want to plot
    NCuts.Add('leptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[0] : Muon_mediumId[0]')
    NCuts.Add('Lepton_miniPFRelIso_all','Lepton_miniPFRelIso_all < 0.2')    
    NCuts.Add('Higgs_msoftdrop','Higgs_msoftdrop > 100 && Higgs_msoftdrop < 150')
    NCuts.Add('Higgs_{}_HbbvsQCD'.format(t),'Higgs_{}_HbbvsQCD > 0.8'.format(t))
    NCuts.Add('Wqq_msoftdrop','Wqq_msoftdrop > 60 && Wqq_msoftdrop < 100')
    NCuts.Add('Wqq_{}_WqqvsQCD'.format(t),'Wqq_{}_WqqvsQCD > 0.8'.format(t))    
    NCuts.Add('DeltaPhi_Lepton_Higgs','abs(DeltaPhi_Lepton_Higgs) > 1.571') #pi/2
    NCuts.Add('DeltaPhi_Lepton_Wqq','abs(DeltaPhi_Lepton_Wqq) < 1.571')
    NCuts.Add('DeltaPhi_Higgs_Wqq','abs(DeltaPhi_Higgs_Wqq) > 1.571')

    return NCuts

def Nminus1Plots(self):

    ##### DEFINITIONS #####

    #W_leptonic transverse mass
    self.a.Define('W_massTran','TransverseMass(MET_pt,Lepton_pt,MET_phi,Lepton_phi)') #Transverse W mass

    #Invariant masses of W/Y/X
    self.a.Define('W_massInv','hardware::InvariantMass({MET_vect,Lepton_vect})') #full invariant mass
    self.a.Define('Y_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')

    nminus1Plots = HistGroup('nminus1Plots')

    #######################################

    taggers = ['particleNetMD']

    for t in taggers:

        nminusGroup = Nminus1Group(t) # dictionary of nodes
        nminusNodes = self.a.Nminus1(nminusGroup) # TIMBER can now be fed the dict and automatically do N-1

        binnings = {
            'Lepton_miniPFRelIso_all':[20,0,1],
            'leptonQuality':[100,0,1],
            'Higgs_msoftdrop':[40,50,250],
            'Higgs_{}_HbbvsQCD'.format(t):[50,0,1],
            'Wqq_msoftdrop':[40,50,250],
            'Wqq_{}_WqqvsQCD'.format(t):[50,0,1],
            'DeltaPhi_Lepton_Higgs':[30,-3.2,3.2],
            'DeltaPhi_Lepton_Wqq':[30,-3.2,3.2],
            'DeltaPhi_Higgs_Wqq':[30,-3.2,3.2],
            'full':{'W_massTran':[30,0,500],'W_massInv':[30,0,500],'Y_mass':{'bkg':[60,100,3000],'sig_low':[30,100,1600],'sig_high':[60,100,3000]}}
        }

        startNminus = self.a.GetActiveNode()
        for n in binnings.keys():
            if n == 'full':
                for var in binnings['full'].keys(): #want to use different binnings for different masses
                    bins = binnings['full'][var]
                    if var == 'Y_mass':
                        if 'XHY' in self.setname:
                            if int(self.setname.split('_')[0].split('-')[2]) > 1000: #expected y mass > 1000   
                                bins = bins['sig_high']
                            else:
                                bins = bins['sig_low']
                        else:
                            bins = bins['bkg']
                    print('N-1: Plotting {}'.format(var))
                    nminus1Plots.Add(var+'_allCuts',nminusNodes[n].DataFrame.Histo1D((var+'_allCuts',var+'_allCuts',bins[0],bins[1],bins[2]),var,'weight__nominal'))
         
            elif n == 'leptonQuality': #use this node to plot electron mva values
                bins = binnings[n]
                self.a.SetActiveNode(nminusNodes[n])
                self.a.Cut('LeptonType','LeptonType == 0')
                self.a.ObjectFromCollection('kinEle','Electron','0')
                no_eta = self.a.GetActiveNode()

                print('N-1: Plotting Electron_mvaFall17V2noIso')
                inner_barrel = self.a.Cut('inner_barrel','abs(kinEle_eta) < 0.8')
                self.a.SetActiveNode(inner_barrel)
                nminus1Plots.Add('kinEle_mvaFall17V2noIso_inner_nminus1',inner_barrel.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_inner_nminus1','kinEle_mvaFall17V2noIso (inner barrel - |eta| < 0.8) - Nminus1',bins[0],bins[1],bins[2]),'kinEle_mvaFall17V2noIso', 'weight__nominal'))

                self.a.SetActiveNode(no_eta)
                outer_barrel = self.a.Cut('outer_barrel','abs(kinEle_eta) > 0.8 && abs(kinEle_eta) < 1.479')
                self.a.SetActiveNode(outer_barrel)
                nminus1Plots.Add('kinEle_mvaFall17V2noIso_outer_nminus1',outer_barrel.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_outer_nminus1','kinEle_mvaFall17V2noIso (outer barrel - 0.8 < |eta| < 1.479) - Nminus1',bins[0],bins[1],bins[2]),'kinEle_mvaFall17V2noIso', 'weight__nominal'))

                self.a.SetActiveNode(no_eta)
                endcap = self.a.Cut('endcap','abs(kinEle_eta) > 1.479 && abs(kinEle_eta) < 2.5')
                self.a.SetActiveNode(endcap)
                nminus1Plots.Add('kinEle_mvaFall17V2noIso_endcap_nminus1',endcap.DataFrame.Histo1D(('kinEle_mvaFall17V2noIso_endcap_nminus1','kinEle_mvaFall17V2noIso (endcap - 1.479 < |eta| < 2.5) - Nminus1',bins[0],bins[1],bins[2]),'kinEle_mvaFall17V2noIso', 'weight__nominal'))
 
            else:
                var = n
                bins = binnings[n]
                print('N-1: Plotting {}'.format(var))
                nminus1Plots.Add(var+'_nminus1',nminusNodes[n].DataFrame.Histo1D((var+'_nminus1',var+'_nminus1',bins[0],bins[1],bins[2]),var,'weight__nominal'))

    return nminus1Plots

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

    filename = 'snapshots/{}_{}_snapshot.txt'.format(setname,year) 

    ana = XHYbbWW(filename,ijob,njobs)

    #ana.ApplyStandardCorrections(snapshot=False)
    #ana.a.MakeWeightCols(extraNominal='' if ana.a.isData else 'genWeight*%s'%ana.GetXsecScale()) #renormalize events 

    outFile = ROOT.TFile.Open('{}_{}_{}_studies.root'.format(ana.setname,ana.year,ana.ijob),'RECREATE')
    outFile.cd()

    kinPlots = KinPlots(ana)
    kinPlots.Do('Write')
    nminus1Plots = Nminus1Plots(ana)  
    nminus1Plots.Do('Write')

    outFile.Close()





 
