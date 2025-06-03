'''
Script to propagate uncertainties on the JEC's and unclustered energy to MET
'''

import ROOT
from TIMBER.Analyzer import *
from TIMBER.Tools.Common import *

def ApplyMETCorrections(self,variation,applyXYshift=True):

    CompileCpp('modules/ApplyMETShiftXY.cc')
    CompileCpp('modules/ApplyMETUncertainties.cc')
    
    # Apply the XY-shift to the nominal MET
    if applyXYshift:
        self.a.Define('Met_MetPhi_XYShifted',f'METXYCorr_Met_MetPhi("{self.setname}","{self.year}",!{str(self.a.isData).lower()},MET_pt,MET_phi,PV_npvs)')
        self.a.Define('MET_pt_XYshifted','Met_MetPhi_XYShifted.first')
        self.a.Define('MET_phi_XYshifted','Met_MetPhi_XYShifted.second')
    else:
        self.a.Define('MET_pt_XYshifted','MET_pt')
        self.a.Define('MET_phi_XYshifted','MET_phi')

    # Now propagate the JEC's to MET, based on the specified variation
    # Variations are : UE_up/down (uncertainty in the unclustered energy propagated to MET), JES_up/down, JER_up/down
    # Jet mass scale/resolution corrections do not affect the pT of jets, so do not need to be accounted for in the MET
    if 'JES' not in variation and 'JER' not in variation:
        # Propagate nominal JES and JER (for MC) to MET
        if self.a.isData:
            self.a.Define('Jet_JER_nom','Jet_pt')

        self.a.Define('MET_corrected',f'propagate_jec_met(Jet_pt,Jet_JES_nom,Jet_JER_nom,Jet_phi,Jet_rawFactor,Jet_muonSubtrFactor,Jet_neEmEF,Jet_chEmEF,MET_pt_XYshifted,MET_phi_XYshifted,{str(self.a.isData).lower()})')

        self.a.Define('MET_px','MET_corrected[0] * cos(MET_corrected[1])')
        self.a.Define('MET_py','MET_corrected[0] * sin(MET_corrected[1])')
        
        if variation == 'UE_up':
            self.a.Define('MET_px_corr','MET_px + MET_MetUnclustEnUpDeltaX')
            self.a.Define('MET_py_corr','MET_py + MET_MetUnclustEnUpDeltaY')

        elif variation == 'UE_down':
            self.a.Define('MET_px_corr','MET_px - MET_MetUnclustEnUpDeltaX')
            self.a.Define('MET_py_corr','MET_py - MET_MetUnclustEnUpDeltaY')

        else:
            self.a.Define('MET_px_corr','MET_px')
            self.a.Define('MET_py_corr','MET_py')

        self.a.Define('MET_pt_corr','sqrt(MET_px_corr*MET_px_corr + MET_py_corr*MET_py_corr)')
        self.a.Define('MET_phi_corr','atan2(MET_py_corr, MET_px_corr)')

    else: # One of the JES/JER variations
        if 'JES' in variation:
            jes = variation
            jer = 'JER_nom'
        elif 'JER' in variation:
            jes = 'JES_nom'
            jer = variation
        # This code is slightly wrong, as I have not yet implemented L1 JEC's to AK4 jets
        self.a.Define('MET_corrected',f'propagate_jec_met(Jet_pt,Jet_{jes},Jet_{jer},Jet_phi,Jet_rawFactor,Jet_muonSubtrFactor,Jet_neEmEF,Jet_chEmEF,MET_pt_XYshifted,MET_phi_XYshifted,{str(self.a.isData).lower()})')

        self.a.Define('MET_pt_corr','MET_corrected[0]')
        self.a.Define('MET_phi_corr','MET_corrected[1]')


