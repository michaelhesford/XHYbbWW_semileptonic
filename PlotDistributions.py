'''
Plot basic kinematic distributions after pre-selection
'''

import ROOT
from TIMBER.Analyzer import Correction, HistGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp
from argparse import ArgumentParser
from XHYbbWW_class import XHYbbWW
from collections import OrderedDict
from math import pi
import time

def KinPlots(self):

    start = self.a.GetActiveNode()
    self.nStart = self.getNweighted()
    self.AddCutflowColumn(self.nStart,'nStart')
    kinPlots=HistGroup('kinPlots')
    
    #Make separate kinematic plots for electron/muon events
    eleEvents=self.a.Cut('eleEvents','LeptonType == 0')
    self.a.SetActiveNode(eleEvents)

    kinPlots.Add('Electron_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_pt','Electron_pt',23,35,1200),'Lepton_pt','weight__nominal'))
    kinPlots.Add('Electron_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_eta','Electron_eta',40,-3,3),'Lepton_eta','weight__nominal'))
    kinPlots.Add('Electron_phi',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_phi','Electron_phi',42,-pi,pi),'Lepton_phi','weight__nominal'))
    kinPlots.Add('Electron_miniPFRelIso_all',self.a.GetActiveNode().DataFrame.Histo1D(('Electron_miniPFRelIso_all','Electron_miniPFRelIso_all',100,0,1),'Lepton_miniPFRelIso_all','weight__nominal'))

    self.a.SetActiveNode(start)
    muEvents=self.a.Cut('muEvents','LeptonType == 1')
    self.a.SetActiveNode(muEvents)

    print('Plotting muon distributions')
    kinPlots.Add('Muon_pt',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_pt','Muon_pt',23,35,1200),'Lepton_pt','weight__nominal'))
    kinPlots.Add('Muon_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_eta','Muon_eta',40,-3,3),'Lepton_eta','weight__nominal'))
    kinPlots.Add('Muon_phi',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_phi','Muon_phi',42,-pi,pi),'Lepton_phi','weight__nominal'))
    kinPlots.Add('Muon_miniPFRelIso_all',self.a.GetActiveNode().DataFrame.Histo1D(('Muon_miniPFRelIso_all','Muon_miniPFRelIso_all',100,0,1),'Lepton_miniPFRelIso_all','weight__nominal'))
    
    #Make the rest of the kinematic distributions
    self.a.SetActiveNode(start)

    print('Plotting MET/neutrino distributions')
    print('Number of events: {}'.format(self.a.DataFrame.Count().GetValue()))
    kinPlots.Add('MET_pt',self.a.GetActiveNode().DataFrame.Histo1D(('MET_pt','MET_pt',28,0,1400),'MET_pt','weight__nominal'))
    kinPlots.Add('MET_phi',self.a.GetActiveNode().DataFrame.Histo1D(('MET_phi','MET_phi',42,-pi,pi),'MET_phi','weight__nominal'))
    
    self.ApplyMETShiftXY() # MET XY shift correction 
    kinPlots.Add('MET_pt_XYshifted',self.a.GetActiveNode().DataFrame.Histo1D(('MET_pt_XYshifted','MET_pt_XYshifted',28,0,1400),'MET_pt_XYshifted','weight__nominal'))
    kinPlots.Add('MET_phi_XYshifted',self.a.GetActiveNode().DataFrame.Histo1D(('MET_phi_XYshifted','MET_phi_XYshifted',42,-pi,pi),'MET_phi_XYshifted','weight__nominal'))
    
    self.a.Define('Neutrino_eta','NeutrinoEta(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass,MET_pt,MET_phi)')
    self.a.Define('Neutrino_eta_XYshifted','NeutrinoEta(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass,MET_pt_XYshifted,MET_phi_XYshifted)')
    kinPlots.Add('Neutrino_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Neutrino_eta','Neutrino_eta',40,-3,3),'Neutrino_eta','weight__nominal'))
    kinPlots.Add('Neutrino_eta_XYshifted',self.a.GetActiveNode().DataFrame.Histo1D(('Neutrino_eta_XYshifted','Neutrino_eta_XYshifted',40,-3,3),'Neutrino_eta_XYshifted','weight__nominal'))
    
    print('Plotting Higgs jet distributions')
    kinPlots.Add('Higgs_pt_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_pt_corr','Higgs_pt_corr',44,300,2500),'Higgs_pt_corr','weight__nominal'))
    kinPlots.Add('Higgs_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_eta','Higgs_eta',40,-3,3),'Higgs_eta','weight__nominal'))
    kinPlots.Add('Higgs_phi',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_phi','Higgs_phi',42,-pi,pi),'Higgs_phi','weight__nominal'))
    kinPlots.Add('Higgs_msoftdrop_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_msoftdrop_corr','Higgs_msoftdrop_corr',50,50,250),'Higgs_msoftdrop_corr','weight__nominal'))
    kinPlots.Add('Higgs_particleNet_mass',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_particleNet_mass','Higgs_particleNet_mass',50,50,250),'Higgs_particleNet_mass','weight__nominal'))
    kinPlots.Add('Higgs_particleNetMD_HbbvsQCD',self.a.GetActiveNode().DataFrame.Histo1D(('Higgs_particleNetMD_HbbvsQCD','Higgs_particleNetMD_HbbvsQCD',50,0,1),'Higgs_particleNetMD_HbbvsQCD','weight__nominal'))

    print('Plotting W jet distributions')
    kinPlots.Add('Wqq_pt_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_pt_corr','Wqq_pt_corr',44,300,2500),'Wqq_pt_corr','weight__nominal'))
    kinPlots.Add('Wqq_eta',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_eta','Wqq_eta',40,-3,3),'Wqq_eta','weight__nominal'))
    kinPlots.Add('Wqq_phi',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_phi','Wqq_phi',42,-pi,pi),'Wqq_phi','weight__nominal'))
    kinPlots.Add('Wqq_msoftdrop_corr',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_msoftdrop_corr','Wqq_msoftdrop_corr',50,50,250),'Wqq_msoftdrop_corr','weight__nominal'))
    kinPlots.Add('Wqq_particleNet_mass',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_particleNet_mass','Wqq_particleNet_mass',50,50,250),'Wqq_particleNet_mass','weight__nominal'))
    kinPlots.Add('Wqq_particleNetMD_WqqvsQCD',self.a.GetActiveNode().DataFrame.Histo1D(('Wqq_particleNetMD_WqqvsQCD','Wqq_particleNetMD_WqqvsQCD',50,0,1),'Wqq_particleNetMD_WqqvsQCD','weight__nominal'))
    
    return kinPlots

if __name__ == "__main__":

    t_start = time.time()

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

    ########## Perform object pre-selection + apply correction ##########

    ana.ApplyJMECorrections(variation='None',jets=['FatJet','Jet'])

    ana.KinematicLepton()
    ana.ApplyTrigs()
    
    ana.Dijets()
    ana.ApplyStandardCorrections(snapshot=False)
    ana.ApplyTopPtReweight('Higgs','Wqq',scale = 0.5, isJet1AK4 = False)
    ana.ApplyLeptonCorrections()

    ana.a.MakeWeightCols(correctionNames=list(ana.a.GetCorrectionNames()), extraNominal='' if ana.a.isData else str(ana.GetXsecScale()))
    #####################################################################

    kinPlots = KinPlots(ana)
    
    #Make HEM plot for 2018
    if ana.year == '18':
        # Affected regions are -1.57 < phi < 0.87, -2.5 < eta < -1.3
        # Binnin constructed so that edges of affected regions correspond to a bin edge
        HEM_hbb = ana.a.GetActiveNode().DataFrame.Histo2D(('HEM_hbb','HEM_hbb',38,-3.32,3.33,40,-2.95,3.05),'Higgs_phi','Higgs_eta','weight__nominal')
        HEM_wqq = ana.a.GetActiveNode().DataFrame.Histo2D(('HEM_wqq','HEM_hbb',38,-3.32,3.33,40,-2.95,3.05),'Wqq_phi','Wqq_eta','weight__nominal')
    outFile = ROOT.TFile.Open('{}_{}_distributions.root'.format(setname,year),'RECREATE')
    outFile.cd()

    print('Writing histograms to file {}_{}_distributions.root'.format(setname,year))
    if ana.year == '18':
        print('Writing HEM')
        HEM_hbb.Write()
        HEM_wqq.Write()
    
    kinPlots.Do('Write')

    outFile.Close()

    print('%s sec'%(time.time()-t_start))
