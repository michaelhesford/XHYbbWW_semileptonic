// Function to propagate JES/JER corrections to MET
// See: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py#L694


RVec<float> propagate_jec_met(RVec<float> jet_pt, RVec<float> jet_jes, RVec<float> jet_jer, RVec<float> jet_phi, RVec<float> jet_rawFactor, 
		              RVec<float> jet_muSubtFactor, RVec<float> neEmEF, RVec<float> chEmEF, float MET_pt, float MET_phi, bool isData) {

    float MET_px = MET_pt * cos(MET_phi);
    float MET_py = MET_pt * sin(MET_phi);

    for (int i=0; i<jet_pt.size(); i++) {
	// Calculated corrected pt (muon-subtracted)
        float jet_pt_noMu = jet_pt[i] * (1 - jet_rawFactor[i]) * (1 - jet_muSubtFactor[i]);
	float muon_pt = jet_pt[i] * (1 - jet_rawFactor[i]) * jet_muSubtFactor[i]; // pt of muons subtracted from the jet

        float jet_pt_noMuL1L2L3 = jet_pt_noMu * jet_jes[i];
        float jet_pt_L1L2L3 = jet_pt_noMuL1L2L3 + muon_pt;

	float jet_pt_fullCorr; // Fully corrected jet pt
        if (isData) {
	    jet_pt_fullCorr = jet_pt_L1L2L3; // No JER smearing applied to MET in data
	}
	else {
	    jet_pt_fullCorr = jet_pt_L1L2L3 * jet_jer[i];
        }

        if (jet_pt_noMuL1L2L3 > 15 && (neEmEF[i] + chEmEF[i]) < 0.9) { // Only apply to MET if corrected pt (w/o muon) is above threshold
            float jet_px_fullCorr = jet_pt_fullCorr * cos(jet_phi[i]);
            float jet_py_fullCorr = jet_pt_fullCorr * sin(jet_phi[i]);

            float jet_px_L1 = jet_pt[i] * cos(jet_phi[i]);
            float jet_py_L1 = jet_pt[i] * sin(jet_phi[i]);

            MET_px = MET_px - (jet_px_fullCorr - jet_px_L1);
            MET_py = MET_py - (jet_py_fullCorr - jet_py_L1);
	}
    }

    float MET_pt_corr = sqrt(MET_px*MET_px + MET_py*MET_py);
    float MET_phi_corr = atan2(MET_py,MET_px);

    return {MET_pt_corr,MET_phi_corr};
}
