'''
Perform selection and make template histograms used for QCD background estimate
'''

import ROOT, time
from TIMBER.Analyzer import HistGroup, Correction
from TIMBER.Tools.Common import CompileCpp, ExecuteCmd
from collections import OrderedDict
import TIMBER.Tools.AutoJME as AutoJME
from XHYbbWW_class import XHYbbWW
from modules.Btag_weight import BTagSF
from modules.ApplyMETUncertainties import ApplyMETUncertainties

def QCD_rpf(self,variation):

    ##### NEW VARIABLES FOR LATER USE #####
    
    #Neutrino eta
    self.a.Define('Neutrino_eta','NeutrinoEta(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass,MET_pt_corr,MET_phi_corr)')

    #Lorentz 4-vectors
    self.a.Define('MET_vect','hardware::TLvector(MET_pt_corr,Neutrino_eta,MET_phi_corr,0)') #neutrino mass negligable, for now assuming MET_eta = 0 (p_z = 0)
    self.a.Define('Wqq_vect','hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)') #use updated mass/pt after JME corrections
    #self.a.Define('Higgs_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)')

    taggers = ['particleNetMD']

    #######################################

    outFile = ROOT.TFile.Open('QCDrpf_{}_{}{}.root'.format(self.setname, self.year, '_'+variation if variation != 'None' else ''),'RECREATE')
    outFile.cd()
    
    # now we want to plot mX vs mY for QCD, ttbar, and signal
    for t in taggers:

        start=self.a.GetActiveNode()
        
        btag_wp = self.cuts['JET']['btagDeepJet']['T'][self.year]
        self.a.Define('bjet_exists',f'does_bjet_exist(kinJet_btagDeepFlavB,{btag_wp})')

        self.ApplyWMass()
        sideband = self.ApplyHiggsMass(sidebands=True) # Control region inside the Higgs mass
        
        #signal region
        self.a.SetActiveNode(sideband)
        SR_node = self.ApplySRorCR('SR',t)
      
        #ttbar-enriched control region
        self.a.SetActiveNode(sideband)
        ttCR_node = self.ApplySRorCR('ttCR',t)

        nodes=OrderedDict([
            ('SR',SR_node),
            ('ttCR',ttCR_node)
        ])

        binsX = [50,0,1]
        #binsY = [1,0,1]
        binsY = [5,55,105]
            
        plot_vars = ['Wqq_particleNetMD_WqqvsQCD','Wqq_particleNet_mass_corr']

        for region in nodes.keys():
            self.a.SetActiveNode(nodes[region])
            print('W tag: Evaluating for {}'.format(region))
            templates = self.a.MakeTemplateHistos(ROOT.TH2F(f'Wtag_{region}', f'Wtag_{region}',binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),plot_vars)
            templates.Do('Write')
    
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
    
    filename='snapshots/{}_{}_snapshot.txt'.format(setname,year)

    selection = XHYbbWW(filename,ijob,njobs)
    selection.nstart = selection.getNweighted()
    selection.AddCutflowColumn(selection.nstart,'nstart')

#####

    selection.ApplyJMECorrections(variation,jets=['FatJet','Jet'])

    selection.KinematicLepton() #NOTE: I changed the isolation from <1 to < 0.1

    selection.ApplyTrigs()

    selection.Dijets()
    selection.ApplyStandardCorrections(snapshot=False)
    selection.ApplyTopPtReweight('Higgs','Wqq',scale = 0.5, isJet1AK4 = False)

    Hbb_wp = selection.cuts['JET']['particleNetMD_HbbvsQCD'][1]
    Wqq_wp = selection.cuts['JET']['particleNetMD_WqqvsQCD']

    selection.ApplyLeptonCorrections()

    # Apply b-tagging scale factors
    selection.a.Define('Jet_vect','hardware::TLvector(Jet_pt_corr,Jet_eta,Jet_phi,Jet_mass)')
    selection.a.Define('Higgs_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)') # Need this here
    selection.a.Define('kinJet_status','KinJetIdxs(Jet_vect, Higgs_vect, Jet_jetId)') #returns vector of indices correponding to -1 (if jet fails kinematic cuts) or index in full full jet collection
    selection.a.SubCollection('kinJet','Jet','kinJet_status != -1') # Sub-collection of jets to consider for b-tagging
    if not selection.a.isData:
        BTagSF(selection,'kinJet',wp='T',mujets_or_comb='comb') # run the function to create the variosu b-tagging corrections/uncertainties

    # MET corrections
    selection.ApplyMETShiftXY()
    ApplyMETUncertainties(selection,variation)

    print('CORRECTIONS:')
    for correction in selection.a.GetCorrectionNames():
        print(correction)

    uncerts_to_corr = {
        'Btag' : ['Btag_bc_correlated','Btag_bc_uncorrelated','Btag_light_correlated','Btag_light_uncorrelated'],
        'MuonIDWeight' : ['MuonIDWeight_uncert_stat','MuonIDWeight_uncert_syst'],
        'MuonRecoWeight' : ['MuonRecoWeight_uncert_stat','MuonRecoWeight_uncert_syst'],
    }

    selection.a.MakeWeightCols(correctionNames=list(selection.a.GetCorrectionNames()), uncerts_to_corr=uncerts_to_corr , extraNominal='' if selection.a.isData else str(selection.GetXsecScale()))

    QCD_rpf(selection,variation)
