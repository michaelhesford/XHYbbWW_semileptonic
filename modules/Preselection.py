'''
Module used to apply standard pre-selection and MC corrections
'''

from modules.ApplyMETUncertainties import ApplyMETCorrections
from modules.Btag_weight import BTagSF


def ApplyPreselection(ana,variation):
    
    # Initial event count
    ana.nStart = ana.getNweighted()
    ana.AddCutflowColumn(ana.nStart,'nStart')

    # Lepton selection/scale factors
    ana.KinematicLepton()
    ana.ApplyLeptonCorrections()
    
    # Apply triggers and corresponding scale factors
    ana.ApplyTrigs()

    # Apply JME corrections, identify Higgs/W candidate jets
    ana.ApplyJMECorrections(variation,jets=['FatJet','Jet'])
    ana.Dijets()
   
    # MET corrections/cuts
    ApplyMETCorrections(ana,variation)
    ana.a.Cut('MET_pt','MET_pt_corr > 30')
    ana.nMET30 = ana.getNweighted()
    ana.AddCutflowColumn(ana.nMET30, 'nMET30')

    # Apply b-tagging scale factors
    ana.a.Define('Jet_vect','hardware::TLvector(Jet_pt_corr,Jet_eta,Jet_phi,Jet_mass)')
    ana.a.Define('Higgs_vect','hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)') # Need this here
    ana.a.Define('bCandIdxs','KinJetIdxs(Jet_vect, Higgs_vect, Jet_jetId)') #returns vector of indices correponding to jets passing the kinematic cuts
    ana.a.SubCollection('kinJet','Jet','bCandIdxs',useTake=True) # Sub-collection of jets to consider for b-tagging
    if not ana.a.isData:
        BTagSF(ana,'kinJet',wp='T',mujets_or_comb='comb') # run the function to create the various b-tagging corrections/uncertainties
    
    # Apply additional standard corrections + top-pt reweight
    ana.ApplyStandardCorrections(snapshot=False)
    ana.ApplyTopPtReweight('Higgs','Wqq',scale = 0.5)
