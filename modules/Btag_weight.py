'''
Script to run b-tagging corrections in timber with correctionlib 
'''
import ROOT
from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *
import correctionlib

def BTagSF(self,jet_coll,wp='T',mujets_or_comb='comb'):
    # self: XHYbbWW object on which to perform SF application 
    # jet_coll: name of jet collection to consider for b-tagging (i.e. "Jet", "kinJet",etc...) - this constitues all jets for which we query the b-tagging status in the analysis
    # wp: working point for tagger used ('L','M', or 'T')
    # mujets_or_comb: choose from "mujets" or "comb" variant of SF's for b/c jets

    # need to do this for it to work
    correctionlib.register_pyroot_binding() 

    # Compile the C++ functions/classes we need 
    CompileCpp('modules/Btag_weight.cc') 

    # load the btagging correction set
    print('Loading BTV SF correction set')
    if '16' in self.year:
        if self.year == '16APV' : cyear = '16preVFP'
        else: cyear = '16postVFP'
    else:
        cyear = self.year
    ROOT.gInterpreter.Declare(f'auto btag_sf_set = correction::CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/20{cyear}_UL/btagging.json.gz");')

    # For DeepJet we use the following scale factors:
    #   - SFbc:
    #       * either "mujets" or "comb" SFs are to be applied 
    #         mujets: SFs derived from QCD-enriched regions; 
    #         comb: additionally also SFs from ttbar enriched regions
    #   - SFlight:
    #       * for light jets (hadron flavor 0) "incl" SFs
    # One common multiplicative SF is to be applied to the nominal prediction.
    # Separate uncertainties are applied for b/c jets and light jets.

    # Get all of the relevant SFs
    ROOT.gInterpreter.Declare(f'auto deepJet_bc_sf = btag_sf_set->at("deepJet_{mujets_or_comb}");') # b/c SFs 
    ROOT.gInterpreter.Declare('auto deepJet_incl_sf = btag_sf_set->at("deepJet_incl");')     # light SFs

    '''
    # Rebin efficiency map
    rebinFactor = 6 if 'XHY' in self.setname else 2
    effFile = ROOT.TFile.Open(f"{path}{self.setname}_{self.year}_bJets_btagEfficiencies.root")
    hAll = effFile.Get('all').RebinY(rebinFactor)
    hTagged = effFile.Get('tagged').RebinY(rebinFactor)
    eff = ROOT.TEfficiency(hTagged,hAll)
    eff.SetName('eff_rebin')
    eff.Write()
    effFile.Close()
    '''

    # Setnames for efficiency maps 
    combined_procs = ['DYJets','ST','WJetsHT','WJetsLepHT','ZJets']
    effname = self.setname
    for name in combined_procs:
        if name in self.setname:
            effname = name
            break


    # Load efficiency maps
    print("Loading b-tagging efficiency maps")
    path = 'root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic/plots/btagEfficiencies/'
    print(f'{path}{effname}_{self.year}_bJets_btagEfficiencies.root')
    print(f'{path}{effname}_{self.year}_cJets_btagEfficiencies.root')
    print(f'{path}{effname}_{self.year}_lightJets_btagEfficiencies.root')
    ROOT.gInterpreter.Declare(f'TEfficiency* eff_b = (TEfficiency*)TFile::Open("{path}{effname}_{self.year}_bJets_btagEfficiencies.root")->Get("eff");')
    ROOT.gInterpreter.Declare(f'TEfficiency* eff_c = (TEfficiency*)TFile::Open("{path}{effname}_{self.year}_cJets_btagEfficiencies.root")->Get("eff");')
    ROOT.gInterpreter.Declare(f'TEfficiency* eff_light = (TEfficiency*)TFile::Open("{path}{effname}_{self.year}_lightJets_btagEfficiencies.root")->Get("eff");')

    # Value to cut on for b-tagging
    wp_val = self.cuts['JET']['btagDeepJet'][wp][self.year]

    # Evaluate nominal b-tagging event weights using method 1A described here:
    # https://btv-wiki.docs.cern.ch/PerformanceCalibration/fixedWPSFRecommendations/#scale-factor-recommendations-for-event-reweighting
    self.a.Define('Btag_weight_nom',f'evaluate_corr(deepJet_bc_sf,deepJet_incl_sf,eff_b,eff_c,eff_light,{jet_coll}_vect,{jet_coll}_hadronFlavour,{jet_coll}_btagDeepFlavB,"{wp}",{wp_val})')
    #self.a.DataFrame.Display([f'Btag_weight_nom']).Print()  # debug

    #Evaluate up/down uncertainty weights
    # A breakdown of full Run2 uncertainties can be found here: 
    # https://btv-wiki.docs.cern.ch/PerformanceCalibration/SFUncertaintiesAndCorrelations/#working-point-based-sfs-fixedwp-sfs
    for corr in ['correlated','uncorrelated']: 
        # Apply a fully correlated uncertainty across the four data taking eras and a separate uncorrelated one
        self.a.Define(f'Btag_weight_uncert_bc_{corr}',f'evaluate_uncert(deepJet_bc_sf,"{corr}",{{5,4}},{{eff_b,eff_c}},{jet_coll}_vect,{jet_coll}_hadronFlavour,{jet_coll}_btagDeepFlavB,"{wp}",{wp_val})') #b/c
        self.a.Define(f'Btag_weight_uncert_light_{corr}',f'evaluate_uncert(deepJet_incl_sf,"{corr}",{{0}},{{eff_light}},{jet_coll}_vect,{jet_coll}_hadronFlavour,{jet_coll}_btagDeepFlavB,"{wp}",{wp_val})') #light
        
        #self.a.DataFrame.Display([f'Btag_weight_uncert_bc_{corr}']).Print()  # debug
        #self.a.DataFrame.Display([f'Btag_weight_uncert_light_{corr}']).Print()  # debug

    # Now use the TIMBER BranchCorrection module to add the weight columns as corrections to be tracked
    btag_weight_nom = Correction('Btag','TIMBER/Framework/src/BranchCorrection.cc',corrtype='corr',mainFunc="evalCorrection") #nominal correction
    self.a.AddCorrection(btag_weight_nom, evalArgs={'val':'Btag_weight_nom'}) 
    # Up/down uncertainties
    for corr in ['correlated','uncorrelated']:
        for bc_l in ['bc','light']:
            btag_weight_uncert = btag_weight_nom.Clone(f'Btag_{bc_l}_{corr}',newMainFunc="evalUncert",newType="uncert")
            self.a.AddCorrection(btag_weight_uncert, evalArgs={'valUp':f'Btag_weight_uncert_{bc_l}_{corr}[0]','valDown':f'Btag_weight_uncert_{bc_l}_{corr}[1]'})

