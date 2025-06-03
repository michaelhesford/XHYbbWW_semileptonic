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

    self.nStart = self.getNweighted()
    self.AddCutflowColumn(self.nStart,'nStart')

    #### DEFINITIONS ####

    self.a.Define('FatJet_particleNetMD_HbbvsQCD','FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD)')
    self.a.Define('FatJet_particleNetMD_WqqvsQCD','(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc)/(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc+FatJet_particleNetMD_QCD)')
    self.a.Define('MET_pt_corr','MET_pt')
    self.a.Define('MET_phi_corr','MET_phi')
    self.ApplyMETShiftXY() # Define XY-shifted MET variables

    #### LEPTON SELECTION ####

    # NOTE: I've loosened the lepton pt/isolation cuts from normal inside HWWmodules.cc (pt > 25, iso < 0.5)
    self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_miniPFRelIso_all,Electron_dxy,Electron_dz,Electron_sip3d)')
    self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_miniPFRelIso_all,Muon_dxy,Muon_dz,Muon_sip3d)')
    self.a.Cut('kinLepton_cut','kinEleIdx != -1 || kinMuIdx != -1') # At least one good lepton
    self.a.Define('LeptonType','LeptonIdx(kinEleIdx,kinMuIdx,Electron_pt,Muon_pt)') # Picks higher pt signal lepton - output = 0 (lepton is electron) or 1 (lepton is muon)
   
    self.nkinLep = self.getNweighted()
    self.AddCutflowColumn(self.nkinLep,'nkinLep')

    self.a.Cut('LeptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP90[kinEleIdx] : Muon_mediumId[kinMuIdx]') # Quality cut

    self.nlepQ = self.getNweighted()
    self.AddCutflowColumn(self.nlepQ,'nlepQ')
    
    self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : Electron_pt[kinEleIdx]')
    self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[kinMuIdx] : Electron_eta[kinEleIdx]')
    self.a.Define('Lepton_abseta','abs(Lepton_eta)')
    self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[kinMuIdx] : Electron_phi[kinEleIdx]')
    self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[kinMuIdx] : Electron_mass[kinEleIdx]')
    self.a.Define('Lepton_miniPFRelIso_all','LeptonType == 1 ? Muon_miniPFRelIso_all[kinMuIdx] : Electron_miniPFRelIso_all[kinEleIdx]')
   
    self.a.Define('kinMuon_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : -10')
    self.a.Define('kinMuon_miniPFRelIso_all','LeptonType == 1 ? Muon_miniPFRelIso_all[kinMuIdx] : -10')
    self.a.Define('kinElectron_pt','LeptonType == 0 ? Electron_pt[kinEleIdx] : -10')
    self.a.Define('kinElectron_miniPFRelIso_all','LeptonType == 0 ? Electron_miniPFRelIso_all[kinEleIdx] : -10')

    self.a.Define('Neutrino_eta','NeutrinoEta(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass,MET_pt_XYshifted,MET_phi_XYshifted)') # Calculate neutrino eta

    #### JET SELECTION ####

    self.a.Define('DijetIdxs','PickDijets(FatJet_pt_corr,FatJet_eta,FatJet_particleNet_mass_corr,FatJet_jetId)')
    self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')
    self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)
 
    self.nDijets = self.getNweighted()
    self.AddCutflowColumn(self.nDijets,'nDijets')

    self.a.Define('HiggsIdx','pick_Higgs_from_dijets(Dijet_particleNetMD_HbbvsQCD, Dijet_particleNetMD_WqqvsQCD, Dijet_phi, Lepton_phi)') # Experimental: select Higgs/W jets
    self.a.Define('WqqIdx','HiggsIdx == 0 ? 1 : 0')

    self.a.ObjectFromCollection('Higgs','Dijet','HiggsIdx')
    self.a.ObjectFromCollection('Wqq','Dijet','WqqIdx')

    #### MORE CORRECTIONS ####

    self.ApplyStandardCorrections(snapshot=False)
    self.ApplyTopPtReweight('Higgs','Wqq',scale = 0.5)
    self.ApplyLeptonCorrections()
    print('Tracking corrections: \n%s'%('\n\t- '.join(list(self.a.GetCorrectionNames())))) # Check which corrections are being tracked
    uncerts_to_corr = {
        'MuonIDWeight' : ['MuonIDWeight_uncert_stat','MuonIDWeight_uncert_syst'],
        'MuonRecoWeight' : ['MuonRecoWeight_uncert_stat','MuonRecoWeight_uncert_syst'],
    }
    self.a.MakeWeightCols(correctionNames=list(self.a.GetCorrectionNames()), uncerts_to_corr=uncerts_to_corr , extraNominal='' if self.a.isData else str(self.GetXsecScale()))
    
    #### LORENTZ 4-VECTORS ####

    self.a.Define('MET_vect','hardware::TLvector(MET_pt_XYshifted,Neutrino_eta,MET_phi_XYshifted,0)')
    self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)')
    self.a.Define('Higgs_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)')
    self.a.Define('Jet_vect','hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass)')
    
    #### APPLY THE USUAL B-TAGGING VETO ####
    
    btag_wp = self.cuts['JET']['btagDeepJet']['T'][self.year]
    self.a.Define('kinJet_status','KinJetIdxs(Jet_vect, Higgs_vect, Jet_jetId)') #returns vector of indices correponding to -1 (if jet fails kinematic cuts) or index in full full jet collection
    self.a.SubCollection('kinJet','Jet','kinJet_status != -1') # Sub-collection of jets to consider for b-tagging
    self.a.Define('bjet_exists',f'does_bjet_exist(kinJet_btagDeepFlavB,{btag_wp})')
    self.a.Cut('btag_veto','!bjet_exists')

    self.nbVeto = self.getNweighted()
    self.AddCutflowColumn(self.nbVeto,'nbVeto')

    #### DEFINE SOME VARIABLES FOR PLOTTING ####

    self.a.Define('DeltaPhi_Higgs_Wqq','abs(hardware::DeltaPhi(Higgs_phi,Wqq_phi))')    
    self.a.Define('DeltaPhi_Lepton_MET','abs(hardware::DeltaPhi(Lepton_phi,MET_phi))')
    self.a.Define('DeltaPhi_Lepton_Higgs','abs(hardware::DeltaPhi(Lepton_phi,Higgs_phi))')
    self.a.Define('DeltaPhi_Lepton_Wqq','abs(hardware::DeltaPhi(Lepton_phi,Wqq_phi))')

    self.a.Define('DeltaR_Higgs_Wqq','hardware::DeltaR(Higgs_vect,Wqq_vect)')
    self.a.Define('DeltaR_Lepton_Higgs','hardware::DeltaR(Lepton_vect,Higgs_vect)')
    self.a.Define('DeltaR_Lepton_Wqq','hardware::DeltaR(Lepton_vect,Wqq_vect)')
    
    self.a.Define('Y_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect})')
    self.a.Define('X_mass','hardware::InvariantMass({Lepton_vect,MET_vect,Wqq_vect,Higgs_vect})')
    self.a.Define('mjj','hardware::InvariantMass({Higgs_vect,Wqq_vect})')
    #self.a.Define('W_massTran','TransverseMass(MET_pt,Lepton_pt,MET_phi,Lepton_phi)') #Transverse W mass
    #self.a.Define('W_massInv','hardware::InvariantMass({MET_vect,Lepton_vect})') #full invariant mass
    #self.a.Define('W_massInv_noMETeta','hardware::InvariantMass({MET_vect_noEta,Lepton_vect})') #full invariant mass
    #self.a.Define('Y_mass_noMETeta','hardware::InvariantMass({Lepton_vect,MET_vect_noEta,Wqq_vect})')
    #self.a.Define('X_mass_noMETeta','hardware::InvariantMass({Lepton_vect,MET_vect_noEta,Wqq_vect,Higgs_vect})')

    #### BEGIN PLOTTING SEQUENCE ####

    start = self.a.GetActiveNode()
    kinPlots=HistGroup('kinPlots')
    
    kinPlots.Add('MET_pt',self.a.GetActiveNode().DataFrame.Histo1D(('MET_pt','MET_pt',40,0,1000),'MET_pt_corr','weight__nominal'))
    kinPlots.Add('Neutrino_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Neutrino_eta','Neutrino_eta',100,-2.5,2.5),'Neutrino_eta','weight__nominal'))
    kinPlots.Add('DeltaPhi_Higgs_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Higgs_Wqq','DeltaPhi_Higgs_Wqq',30,0,3.2),'DeltaPhi_Higgs_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Higgs_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Higgs_Wqq','DeltaR_Higgs_Wqq',50,0,5),'DeltaR_Higgs_Wqq','weight__nominal'))
    kinPlots.Add('Higgs_msoftdrop_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_msoftdrop_corr','Higgs_msoftdrop_corr',40,50,250),'Higgs_msoftdrop_corr','weight__nominal'))
    kinPlots.Add('Wqq_msoftdrop_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_msoftdrop_corr','Wqq_msoftdrop_corr',40,50,250),'Wqq_msoftdrop_corr','weight__nominal'))
    kinPlots.Add('Higgs_particleNet_mass',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_particleNet_mass','Higgs_particleNet_mass',40,50,250),'Higgs_particleNet_mass','weight__nominal'))
    kinPlots.Add('Wqq_particleNet_mass',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_particleNet_mass','Wqq_particleNet_mass',40,50,250),'Wqq_particleNet_mass','weight__nominal'))
    kinPlots.Add('Higgs_particleNetMD_HbbvsQCD',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_particleNetMD_HbbvsQCD','Higgs_particleNetMD_HbbvsQCD',50,0,1),'Higgs_particleNetMD_HbbvsQCD','weight__nominal'))
    kinPlots.Add('Wqq_particleNetMD_WqqvsQCD',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_particleNetMD_WqqvsQCD','Wqq_particleNetMD_WqqvsQCD',50,0,1),'Wqq_particleNetMD_WqqvsQCD','weight__nominal'))

    # Make separate kinematic plots for electron/muon events
    eleEvents=self.a.Cut('eleEvents','LeptonType == 0')
    self.a.SetActiveNode(eleEvents)
    
    kinPlots.Add('Electron_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_pt','Electron_pt',50,10,1000),'Lepton_pt','weight__nominal'))
    kinPlots.Add('Electron_miniPFRelIso_all',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_miniPFRelIso_all','Electron_miniPFRelIso_all',20,0,1),'Lepton_miniPFRelIso_all','weight__nominal'))
    kinPlots.Add('DeltaPhi_Electron_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Electron_Higgs','DeltaPhi_Electron_Higgs',30,0,3.2),'DeltaPhi_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaPhi_Electron_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Electron_Wqq','DeltaPhi_Electron_Wqq',30,0,3.2),'DeltaPhi_Lepton_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Electron_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Electron_Higgs','DeltaR_Electron_Higgs',50,0,5),'DeltaR_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaR_Electron_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Electron_Wqq','DeltaR_Electron_Wqq',50,0,5),'DeltaR_Lepton_Wqq','weight__nominal'))
    
    self.a.SetActiveNode(start)
    muEvents=self.a.Cut('muEvents','LeptonType == 1')
    self.a.SetActiveNode(muEvents)

    kinPlots.Add('Muon_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_pt','Muon_pt',50,10,1000),'Lepton_pt','weight__nominal')) 
    kinPlots.Add('Muon_miniPFRelIso_all',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_miniPFRelIso_all','Muon_miniPFRelIso_all',20,0,1),'Lepton_miniPFRelIso_all','weight__nominal'))
    kinPlots.Add('DeltaPhi_Muon_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Muon_Higgs','DeltaPhi_Muon_Higgs',30,0,3.2),'DeltaPhi_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaPhi_Muon_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaPhi_Muon_Wqq','DeltaPhi_Muon_Wqq',30,0,3.2),'DeltaPhi_Lepton_Wqq','weight__nominal'))
    kinPlots.Add('DeltaR_Muon_Higgs',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Muon_Higgs','DeltaR_Muon_Higgs',50,0,5),'DeltaR_Lepton_Higgs','weight__nominal'))
    kinPlots.Add('DeltaR_Muon_Wqq',self.a.GetActiveNode().DataFrame.Histo1D(('DeltaR_Muon_Wqq','DeltaR_Muon_Wqq',50,0,5),'DeltaR_Lepton_Wqq','weight__nominal'))

    self.a.SetActiveNode(start)

    cutflowInfo = OrderedDict([
        ('nStart',self.nStart),
        ('nkinLep',self.nkinLep),
        ('nlepQ',self.nlepQ),
        ('nDijets',self.nDijets),
        ('nbVeto',self.nbVeto)
    ])
    nLabels = len(cutflowInfo)
    hCutflow = ROOT.TH1F('cutflow', 'Number of events after each cut', nLabels, 0.5, nLabels+0.5)
    nBin = 1
    for label, value in cutflowInfo.items():
        hCutflow.GetXaxis().SetBinLabel(nBin, label)
        hCutflow.AddBinContent(nBin, value)
        nBin += 1

    return kinPlots, hCutflow

def Nminus1Group(self):

    WtagLP = {'16APV':0.637,'16':0.642,'17':0.579,'18':0.59}
    WtagMP = {'16APV':0.845,'16':0.842,'17':0.810,'18':0.82}

    NCuts=CutGroup('NCuts')
    #Cuts to apply for all plots - will make N-1 for all of these (except lepton quality)

    NCuts.Add('Lepton_pt','Lepton_pt > 35')
    NCuts.Add('Lepton_miniPFRelIso_all','Lepton_miniPFRelIso_all < 0.1')    
    NCuts.Add('Higgs_particleNet_mass','Higgs_particleNet_mass > 100 && Higgs_particleNet_mass < 150')
    NCuts.Add('Higgs_particleNetMD_HbbvsQCD','Higgs_particleNetMD_HbbvsQCD > 0.98')
    NCuts.Add('Wqq_particleNet_mass','Wqq_particleNet_mass > 55 && Wqq_particleNet_mass < 105')
    NCuts.Add('Wqq_particleNetMD_WqqvsQCD',f'Wqq_particleNetMD_WqqvsQCD > {WtagLP[self.year]}')    
    NCuts.Add('MET_pt','MET_pt_XYshifted > 25')
    #NCuts.Add('DeltaR_Lepton_Higgs','abs(DeltaR_Lepton_Higgs) > 0.8')
    #NCuts.Add('DeltaR_Higgs_BJet','!ak4_exists')
    #NCuts.Add('DeltaEta_Higgs_Wqq','DeltaEta_Higgs_Wqq < 1.3') #experimental, increases efficiency of JetHT triggers
    #NCuts.Add('DeltaPhi_Lepton_Higgs','abs(DeltaPhi_Lepton_Higgs) > 0.5')
    #NCuts.Add('DeltaPhi_Lepton_Wqq','abs(DeltaPhi_Lepton_Wqq) > 1 && abs(DeltaPhi_Lepton_Wqq) < 2.5')

    return NCuts

def Nminus1Plots(self):

    nminus1Plots = HistGroup('nminus1Plots')
    nminusGroup = Nminus1Group(self) # dictionary of nodes
    nminusNodes = self.a.Nminus1(nminusGroup) # TIMBER can now be fed the dict and automatically do N-1
    
    binnings = {
        'Lepton_pt': { 
            'kinElectron_pt':[50,10,1000],
            'kinMuon_pt':[50,10,1000]
        },
        'Lepton_miniPFRelIso_all': {
            'kinElectron_miniPFRelIso_all':[20,0,1],
            'kinMuon_miniPFRelIso_all':[20,0,1]
        },
        'Higgs_particleNet_mass': {
            'Higgs_msoftdrop_corr': [40,50,250],
            'Higgs_particleNet_mass': [40,50,250]
        },
        'Wqq_particleNet_mass': {
            'Wqq_msoftdrop_corr': [40,50,250],
            'Wqq_particleNet_mass': [40,50,250]
        },
        'Higgs_particleNetMD_HbbvsQCD':[50,0,1],
        'Wqq_particleNetMD_WqqvsQCD':[50,0,1],
        'MET_pt':[40,0,1000],
        'full': {
             'X_mass':[76,600,4400],
             'Y_mass':[60,100,3000],
             'mjj':[76,600,4400],
             'DeltaPhi_Lepton_Wqq':[30,0,3.2],
             'DeltaPhi_Lepton_Higgs':[30,0,3.2],
             'DeltaPhi_Higgs_Wqq':[30,0,3.2],
             'DeltaR_Higgs_Wqq':[50,0,5],
             'DeltaR_Lepton_Wqq':[50,0,5],
             'DeltaR_Lepton_Higgs':[50,0,5],
        }
    }

    startNminus = self.a.GetActiveNode()
    for n in binnings.keys():
        if n == 'full':
            for var in binnings['full'].keys():
                bins = binnings['full'][var]
                print('N-1: Plotting {}'.format(var))
                nminus1Plots.Add(var+'_allCuts',nminusNodes[n].DataFrame.Histo1D((var+'_allCuts',var+'_allCuts',bins[0],bins[1],bins[2]),var,'weight__nominal'))
        elif isinstance(binnings[n],dict):
            for var in binnings[n].keys():
                bins = binnings[n][var]
                print('N-1: Plotting {}'.format(var))
                nminus1Plots.Add(var+'_nminus1',nminusNodes[n].DataFrame.Histo1D((var+'_nminus1',var+'_nminus1',bins[0],bins[1],bins[2]),var,'weight__nominal'))
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

    ana.ApplyTrigs # Apply triggers
    ana.ApplyJMECorrections(variation='None',jets=['FatJet']) # Apply JME corrections to AK8 jets

    kinPlots, hCutflow = KinPlots(ana)
    nminus1Plots = Nminus1Plots(ana) 

    outFile = ROOT.TFile.Open('{}_{}_studies.root'.format(setname,year),'RECREATE')
    outFile.cd()

    hCutflow.Write()
    kinPlots.Do('Write')
    nminus1Plots.Do('Write')

    outFile.Close()





 
