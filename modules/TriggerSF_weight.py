'''
Script to appy trigger corrections in timber
Muon SF's are taken from POG measurements and applied using correctionlib. See: 
EGamma SF's are taken from HIG-24-008 and accessed from a .pkl file
'''
import ROOT
from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *
import correctionlib
import pickle

correctionlib.register_pyroot_binding() 
CompileCpp('modules/TriggerSF_weight.cc') # Do this here to avoid compiling twice
corrTemp = Correction(f'corrTemp','TIMBER/Framework/src/BranchCorrection.cc',corrtype='weight',mainFunc='evalWeight') # make dummy correction 

def MuonSF(self):

    # load the btagging correction set
    print('Loading MUO SF correction set')
    if self.year == '16APV':
        hipm = 'HIPM_'
        fyear = '16'
    else:
        hipm = ''
        fyear = self.year

    ROOT.gInterpreter.Declare(f'auto muo_sf_set = correction::CorrectionSet::from_file("corrections/Efficiencies_muon_generalTracks_Z_Run20{fyear}_UL_{hipm}SingleMuonTriggers_schemaV2.json");')

    for pt_cat in ['low','high']:

        if '16' in self.year:
            if pt_cat == 'low':
                corrName = 'NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium'
            elif pt_cat == 'high':
                corrName = 'NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose'
        elif self.year == '17':
            if pt_cat == 'low':
                corrName = 'NUM_IsoMu27_DEN_CutBasedIdMedium_and_PFIsoMedium'
            elif pt_cat == 'high':
                corrName = 'NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose'
        elif self.year == '18':
            if pt_cat == 'low':
                corrName = 'NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium'
            elif pt_cat == 'high':
                corrName = 'NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose'

        ROOT.gInterpreter.Declare(f'auto muo_sf_{pt_cat} = muo_sf_set->at("{corrName}");')
        
        self.a.Define(f'MuonTriggerSF_{pt_cat}Pt_vec',f'make_muo_sf_vec(muo_sf_{pt_cat}, Lepton_abseta, Lepton_pt, LeptonType, "{pt_cat}")')

        muon_trigger_weight = corrTemp.Clone(f'MuoTrigSF_{pt_cat}Pt')
        self.a.AddCorrection(muon_trigger_weight, evalArgs = {'val':f'MuonTriggerSF_{pt_cat}Pt_vec[0]', 'valUp':f'MuonTriggerSF_{pt_cat}Pt_vec[1]', 'valDown':f'MuonTriggerSF_{pt_cat}Pt_vec[2]'})


def ElectronSF(self):
    
    with open('corrections/EGMtrig_SF.pkl','rb') as infile:
        sf_set = pickle.load(infile)

    cyear = 'UL'+self.year if 'APV' not in self.year else 'UL16'

    # This is not very elegant, but it works...
    vals_nom = str(sf_set[cyear]['nominal']).replace('\n ',',').replace('[','{').replace(']','}').replace('  ',' ').replace(' ',',')
    vals_up = str(sf_set[cyear]['up']).replace('\n ',',').replace('[','{').replace(']','}').replace('  ',' ').replace(' ',',')
    vals_down = str(sf_set[cyear]['down']).replace('\n ',',').replace('[','{').replace(']','}').replace('  ',' ').replace(' ',',')

    ROOT.gInterpreter.Declare(f'float sf_vals_nom[3][5] = {vals_nom};')
    ROOT.gInterpreter.Declare(f'float sf_vals_up[3][5] = {vals_up};')
    ROOT.gInterpreter.Declare(f'float sf_vals_down[3][5] = {vals_down};')

    self.a.Define('ElectronTriggerSF_vec','make_egm_sf_vec(LeptonType, Lepton_eta, Lepton_pt, sf_vals_nom, sf_vals_up, sf_vals_down)')

    egamma_trigger_weight = corrTemp.Clone('EgmTrigSF')
    self.a.AddCorrection(egamma_trigger_weight, evalArgs = {'val':f'ElectronTriggerSF_vec[0]', 'valUp':f'ElectronTriggerSF_vec[1]', 'valDown':f'ElectronTriggerSF_vec[2]'})   

   


