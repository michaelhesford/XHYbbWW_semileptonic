#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <stdio.h>
#include <cmath>

RVec<int> PickDijets(RVec<float> pt, RVec<float> eta, RVec<float> mass) {
    // We are looking for a lead jet separated by a sublead jet by at least 90 degrees
    int jet0Idx = -1;
    int jet1Idx = -1;
    // Since the vectors are ordered by pT, the first jet should roughly be our Higgs. 
    for (int ijet=0; ijet<eta.size(); ijet++) {
        if (std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50 && pt[ijet] > 300) {
            // we've found our lead jet, break loop
            jet0Idx = ijet;
            break;
        }
    }
    // if no lead jet found somehow, return 
    if (jet0Idx == -1) {
        return {-1, -1};
    }
    // now loop over the remaining jets,
    for (int ijet=jet0Idx+1; ijet<eta.size(); ijet++) {
        if (std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50 && pt[ijet] > 300) { //&& std::abs(eta[ijet] - eta[jet0Idx]) < 1.3) { //experimental delta eta cut
            jet1Idx = ijet;
            break;
        }
    }
    if (jet1Idx == -1) {
        return {-1,-1};
    }
    else {
        return {jet0Idx,jet1Idx};
    }
}

int PickJetB(RVec<ROOT::Math::PtEtaPhiMVector> Jet_vect, ROOT::Math::PtEtaPhiMVector Lepton_vect, RVec<int> jetID, RVec<float> btag, float wp) {
    // Pick index of highest pt b-tagged ak4 jet with deltaR (jet,lepton) < 2
    int Jet_idx = -1;
    for (int i=0; i<Jet_vect.size(); i++) {
        if (Jet_vect[i].Pt() > 25 && std::abs(Jet_vect[i].Eta()) < 2.4 && btag[i] > wp && jetID[i] == 6 && abs(hardware::DeltaPhi(Jet_vect[i].Phi(),Lepton_vect.Phi())) < 2) {
            Jet_idx = i;
	    break;
	}
    }
    return Jet_idx;
}

int PickAK8(RVec<ROOT::Math::PtEtaPhiMVector> FatJet_vect, ROOT::Math::PtEtaPhiMVector BJet_vect, ROOT::Math::PtEtaPhiMVector Lepton_vect) {
    int Jet_idx = -1;
    for (int idx=0; idx<FatJet_vect.size(); idx++) {
        if (FatJet_vect[idx].Pt() > 300 && std::abs(FatJet_vect[idx].Eta()) < 2.4 && abs(hardware::DeltaPhi(FatJet_vect[idx].Phi(),Lepton_vect.Phi())) > 2 && hardware::DeltaR(FatJet_vect[idx],BJet_vect) > 0.8) {
            Jet_idx = idx;
            break;
        }
    }
    return Jet_idx;
}

bool does_bjet_exist(ROOT::Math::PtEtaPhiMVector Hbb_vec, RVec<ROOT::Math::PtEtaPhiMVector> ak4_vec, RVec<float> ak4_btag, float wp) {
    // Identify b-tagged AK4 jets which are outside the radius of the AK8 Higgs jets
    std::size_t NJ = ak4_vec.size();
    bool ak4 = false;
    for (int i=0; i<NJ; i++) {
	if (ak4_btag[i] > wp && hardware::DeltaR(ak4_vec[i], Hbb_vec) > 1.2) {  
	    ak4 = true;
	    break;
	}
    }
    return ak4;
}

bool isLeptonPreselected(int nElectron, RVec<float> elePt, RVec<float> eleEta, RVec<float> eleIso, int nMuon, RVec<float> muPt, RVec<float> muEta, RVec<float> muIso){
    if (nElectron < 1 && nMuon < 1) {return false;}
    bool isLepPre = false;
    for (int idx=0; idx<elePt.size(); idx++){
        if (elePt[idx] > 20 && std::abs(eleEta[idx]) < 2.5 && eleIso[idx] < 1){
            isLepPre = true;
            return isLepPre;
        }
    }
    for (int idx=0; idx<muPt.size(); idx++){
        if (muPt[idx] > 20 && std::abs(muEta[idx]) < 2.4 && muIso[idx] < 1){
            isLepPre = true;
            break;
        }
    }   
    return isLepPre;
}

bool isJetPreselected(int nFatJet, RVec<float> FatJet_pt, RVec<float> FatJet_eta, RVec<float> FatJet_mass) {
    bool isJetPre = false;
    // perform standard jet pre-selection
    if (nFatJet < 2) {return isJetPre;}
    int idx1 = -1;
    for (int idx=0; idx<FatJet_pt.size(); idx++) {
        //check if first ak8 jets exists
        if (FatJet_pt[idx] > 200 && std::abs(FatJet_eta[idx]) < 2.4 && FatJet_mass[idx] > 50) {
            idx1 = idx;
            break;
        }
    }
    if (idx1 == -1) {return isJetPre;}
    for (int idx=idx1+1; idx<FatJet_pt.size(); idx++) {
        if (FatJet_pt[idx] > 200 && std::abs(FatJet_eta[idx]) < 2.4 && FatJet_mass[idx] > 50) {
            isJetPre = true;
            break;
        }
    }
    return isJetPre;
}

int kinElectron(RVec<float> pt, RVec<float> eta, RVec<float> iso){
    // Select highest pt signal-like electron
    int eleIdx = -1;
    for (int idx=0; idx<pt.size(); idx++) {
        if (pt[idx] > 25 && std::abs(eta[idx]) < 2.5) { //  && iso[idx] < 0.1) {
            eleIdx = idx;
            break;
        }
    }
    return eleIdx;
}

int kinMuon(RVec<float> pt, RVec<float> eta, RVec<float> iso){
    // Select highest pt signal-like muon
    int muIdx = -1;
    for (int idx=0; idx<pt.size(); idx++) {
        if (pt[idx] > 25 && std::abs(eta[idx]) < 2.4) { // && iso[idx] < 0.1) {
            muIdx = idx;
            break;
        }
    }
    return muIdx;
}

int kinLepton(RVec<float> pt, RVec<float> eta){
    // Select highest pt signal-like lepton
    int Idx = -1;
    for (int idx=0; idx<pt.size(); idx++) {
        if (pt[idx] > 25 && std::abs(eta[idx]) < 2.5) {
            Idx = idx;
            break;
        }
    }
    return Idx;
}

int PNetLepton(RVec<float> pt, RVec<float> eta, RVec<float> iso){
    // Select highest pt signal-like lepton
    int Idx = -1;
    for (int idx=0; idx<pt.size(); idx++) {
        if (pt[idx] > 30 && std::abs(eta[idx]) < 2.4 && iso[idx] < 0.1) {
            Idx = idx;
            break;
        }
    }
    return Idx;
}

int leptonType_studies(RVec<float> electron_pt, RVec<float> muon_pt) {
    // Find out if highest lepton in event is election or muon (return 0 if election, 1 if muon)
    // Used in XHYbbWW_studies.py
    if (electron_pt.size() == 0) {
        return 1;
    }
    else if (muon_pt.size() == 0) {
        return 0;
    }
    else if (muon_pt[0] > electron_pt[0]) {
        return 1;
    }
    else {
        return 0;
    }
}

int LeptonIdx(int electronIdx, int muonIdx, RVec<float> electron_pt, RVec<float> muon_pt) {
    // Decide whether an event is a lepton or muon type
    // Pick highest pt lepton amongst previously selected signal-like electrons and muons
    int lepIdx = -1;
    if (electronIdx == -1) {
        lepIdx = 1; //its a muon
    }
    else if (muonIdx == -1) {
        lepIdx = 0; //its an electron
    }
    else {
        if (electron_pt[electronIdx] > muon_pt[muonIdx]) {
            lepIdx = 0;
        }
        else {
            lepIdx = 1;
        }
    }
    return lepIdx;
}

RVec<int> FindMothersPdgId(RVec<int> genpart_id, RVec<int> selected_genpart_mother_indices){
    // For gen particles - return a vector of pdg Id's corresponding to the mother particles of an input list of particles 
    std::size_t Ni = selected_genpart_mother_indices.size();
    RVec<int> mother_pdgids(Ni);
    for(std::size_t i=0; i<Ni; i++) {
        mother_pdgids[i] = genpart_id[selected_genpart_mother_indices[i]];
    }
    return mother_pdgids;

}

float TransverseMass(float MET_pt, float obj_pt, float MET_phi, float obj_phi) {
    // Calculates the transverse mass of a composite object (typically MET + something else)
    // Alternative to TIMBER's hardware::TransverseMass()
    return sqrt(2.0*MET_pt*obj_pt*(1-cos(hardware::DeltaPhi(MET_phi,obj_phi))));
}

float DeltaR(float eta1, float phi1, float eta2, float phi2) {
    // Calculates the difference in the R parameter between two objects
    // Alternative to TIMBER's hardware::DeltaR(ROOT::Math::PtEtaPhiMVector v1, ROOT::Math::PtEtaPhiMVector v2)
    // DeltaR = sqrt[(eta1 - eta2)^2 + (phi1 - phi2)^2]
    float deta = eta2 - eta1;
    float dphi = phi2 - phi1;
    float deltaR = std::abs(sqrt(deta*deta + dphi*dphi));
    return deltaR;
}

int JetFlavor_ttbar(float Jet_eta, float Jet_phi, RVec<float> Gen_eta, RVec<float> Gen_phi, RVec<int> pdgId, RVec<int> motherIdx) {
    // EXPERIMENTAL: Assigns true flavor to ttbar MC AK8 jets based on gen particles inside the jet cone (R=0.8):
    // To be used on a single jet (not an RVec collection)
    // Output: a number corresponding to jet flavor
    // 0: merged top jet
    // 1: merged W jet
    // 2: bq jet
    // 3: unmerged jet
    std::vector<int> b, q; // bottom quarks from top decays, quarks from W decays
    for (int igen = 0; igen<Gen_eta.size(); igen++) { // loop over gen particles
        int genId = pdgId[igen];
        int mother_genId = pdgId[motherIdx[igen]];
        int grandmother_genId = pdgId[motherIdx[motherIdx[igen]]];
        // Pull out all b quarks from top decays, quarks from W decays (w/ grandmother top)
        if (std::abs(genId) == 5 && std::abs(mother_genId) == 6) {b.push_back(igen);}
        else if (std::abs(genId) < 6 && std::abs(genId) > 0 && std::abs(mother_genId) == 24) {
            // make sure mother of W is not just another W
            int greatgrandmother_genId = pdgId[motherIdx[motherIdx[motherIdx[igen]]]];
            if (std::abs(grandmother_genId) == 6) {q.push_back(igen);}
            else if (std::abs(grandmother_genId) == 24 && std::abs(greatgrandmother_genId) == 6) {q.push_back(igen);}
        }
    }
    int jetId = 3; // assume unmerged jet
    int n_bfromt = 0; // b quarks from top decay inside jet cone
    int n_qfromW = 0; // quarks form W decay insdie jet cone
    // Impose delta R requirements
    for (int i : b) {
        if (DeltaR(Jet_eta, Jet_phi, Gen_eta[i], Gen_phi[i]) < 0.8) {n_bfromt += 1;}
    }
    for (int i : q) {
        if (DeltaR(Jet_eta, Jet_phi, Gen_eta[i], Gen_phi[i]) < 0.8) {n_qfromW += 1;}
    }
    if (n_qfromW >= 1) {
        if (n_qfromW >= 2) {
            if (n_bfromt >= 1) {jetId = 0;} // merged top jet
            else {jetId = 1;} // merged W jet
        }
        else if (n_bfromt >= 1) {jetId = 2;} // bq jet        
    }
    return jetId;
}

RVec<int> JetFlavor_ttbar_vec(RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Gen_eta, RVec<float> Gen_phi, RVec<int> pdgId, RVec<int> motherIdx) {
    // EXPERIMENTAL: Assigns true flavor to ttbar MC AK8 jets based on gen particles inside the jet cone (R=0.8):
    // To be used on an RVec collection of jets (like FatJet)
    // Output: a number corresponding to jet flavor
    // 0: merged top jet
    // 1: merged W jet
    // 2: bq jet
    // 3: unmerged jet
    std::vector<int> b, q; // bottom quarks from top decays, quarks from W decays
    for (int igen = 0; igen<Gen_eta.size(); igen++) { // loop over gen particles
        int genId = pdgId[igen];
        int mother_genId = pdgId[motherIdx[igen]];
        int grandmother_genId = pdgId[motherIdx[motherIdx[igen]]];
        // Pull out all b quarks from top decays, quarks from W decays (w/ grandmother top)
        if (std::abs(genId) == 5 && std::abs(mother_genId) == 6) {b.push_back(igen);} 
        else if (std::abs(genId) < 6 && std::abs(genId) > 0 && std::abs(mother_genId) == 24) {
            // make sure mother of W is not just another W
            int greatgrandmother_genId = pdgId[motherIdx[motherIdx[motherIdx[igen]]]];
            if (std::abs(grandmother_genId) == 6) {q.push_back(igen);}
            else if (std::abs(grandmother_genId) == 24 && std::abs(greatgrandmother_genId) == 6) {q.push_back(igen);}           
        }
    }
    RVec<int> jetId(Jet_eta.size());
    for (int ijet = 0; ijet<Jet_eta.size(); ijet++) { // Loop through all jets in the event
        int flav = 3; // assume unmerged jet
        int n_bfromt = 0; // b quarks from top decay inside jet cone
        int n_qfromW = 0; // quarks form W decay insdie jet cone
        // Impose delta R requirements
        for (int i : b) {
            if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Gen_eta[i], Gen_phi[i]) < 0.8) {n_bfromt += 1;}
        }
        for (int i : q) {
            if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Gen_eta[i], Gen_phi[i]) < 0.8) {n_qfromW += 1;}
        }
        if (n_qfromW >= 1) {
            if (n_qfromW >= 2) {
                if (n_bfromt >= 1) {flav = 0;} // merged top jet
                else {flav = 1;} // merged W jet
            }
            else if (n_bfromt >= 1) {flav = 2;} // bq jet        
        }
        jetId[ijet] = flav;
    }
    return jetId;
}

RVec<int> JetFlavor_signal_vec(RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Gen_eta, RVec<float> Gen_phi, RVec<int> pdgId, RVec<int> motherIdx) {
    // EXPERIMENTAL: Assigns true flavor to signal MC AK8 jets based on gen particles inside the jet cone (R=0.8): 
    // Output: an RVec of numbers corresponding to jet flavors
    //     1: merged Wqq jet
    //     4: merged Hbb jet
    //     3: unmerged jet
    std::vector<int> b, q; // bottom quarks from H decays, quarks from W decays
    RVec<int> jetId(Jet_eta.size()); // vector of jet flavors
    for (int igen = 0; igen<Gen_eta.size(); igen++) { // loop over gen particles
        int genId = pdgId[igen];
        int mother_genId = pdgId[motherIdx[igen]];
        // First search for bottom quarks from Higgs decay
        if (std::abs(genId) == 5 && std::abs(mother_genId) == 25) {b.push_back(igen);}
        // Then search for quarks from W decay
        else if (std::abs(genId) < 6 && std::abs(genId) > 0 && std::abs(mother_genId) == 24) {q.push_back(igen);}
    }
    for (int ijet = 0; ijet<Jet_eta.size(); ijet++) { // loop through all jets in the event
        int flav = 3; // assume unmerged jet
        int n_bfromH = 0; // number of b quarks from H decay inside jet cone
        int n_qfromW = 0; // number of quarks from W decay inside jet cone
        // Impose delta R requirements
        for (int i : b) {
            if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Gen_eta[i], Gen_phi[i]) < 0.8) {n_bfromH += 1;}
        }
        for (int i : q) {
            if (DeltaR(Jet_eta[ijet], Jet_phi[ijet], Gen_eta[i], Gen_phi[i]) < 0.8) {n_qfromW += 1;}
        }
        if (n_bfromH >= 2) {flav = 4;} // First check if merged Higgs
        else if (n_qfromW >= 2) {flav = 1;} // Merged W
        jetId[ijet] = flav;
    }
    return jetId;
}

RVec<float> deltaR_from_reference(RVec<ROOT::Math::PtEtaPhiMVector> jet_vec, ROOT::Math::PtEtaPhiMVector ref_vec) {
    // Takes a collection of objects and returns a vector of DeltaR values between those objects and some reference object
    // i.e. calculates DeltaR between a collection of b-tagged AK4 jets and the Higgs candidate jet
    std::size_t nJets = jet_vec.size();
    RVec<float> DeltaR(nJets);
    for (int i=0; i<nJets; i++) {
        float dR = hardware::DeltaR(ref_vec,jet_vec[i]);
        DeltaR[i] = dR;	
    }
    return DeltaR;
}

RVec<float> deltaPhi_from_reference(RVec<float> jet_phi, float ref_phi) {
    // Takes a collection of objects and returns a vector of DeltaPhi values between those objects and some reference object
    std::size_t nJets = jet_phi.size();
    RVec<float> DeltaPhi(nJets);
    for (int i=0; i<nJets; i++) {
        float dPhi = abs(hardware::DeltaPhi(ref_phi,jet_phi[i]));
        DeltaPhi[i] = dPhi;
    }
    return DeltaPhi;
}

RVec<bool> is_b_in_AK4(RVec<int> GenPart_pdgId, RVec<ROOT::Math::PtEtaPhiMVector> GenPart_vect, RVec<ROOT::Math::PtEtaPhiMVector> AK4_vect) {
    // Identifies which gen particles are b quarks inside one constituent of an input jet collection
    std::size_t nPart = GenPart_pdgId.size();
    RVec<bool> is_b(nPart);
    for (int i=0; i<nPart; i++) {
	bool b = false;
        if (std::abs(GenPart_pdgId[i]) == 5) {
    	    for (int j=0; j<AK4_vect.size(); j++) {
	        if (hardware::DeltaR(GenPart_vect[i], AK4_vect[j]) < 0.8) { 	
	            b = true;
		    break;
		}
            }
	}
	is_b[i] = b;
    }
    return is_b;
}

RVec<int> does_AK4_have_b(RVec<int> GenPart_pdgId, RVec<ROOT::Math::PtEtaPhiMVector> GenPart_vect, RVec<ROOT::Math::PtEtaPhiMVector> AK4_vect, RVec<float> AK4_btag, float wp) {
    // Identifies which gen particles are b quarks inside one constituent of an input jet collection, and returns the matching jets
    std::size_t nJet = AK4_btag.size();
    RVec<int> is_bjet(nJet);
    for (int j=0; j<nJet; j++) {
	int bjet = -1;
	if (AK4_btag[j] > wp) {
            for (int i=0; i<GenPart_pdgId.size(); i++) {
                if (std::abs(GenPart_pdgId[i]) == 5 && hardware::DeltaR(GenPart_vect[i],AK4_vect[j]) < 0.8) {
                    bjet = j;
                    break;
                }
            }
	}
        is_bjet[j] = bjet;
    }
    return is_bjet;
}

RVec<int> PickBfromT(RVec<int> pdgId, RVec<int> mother_idx) {
    RVec<int> b_idxs(pdgId.size());
    for (int i=0; i<pdgId.size(); i++) {
	int b_from_t = -1;
        if (std::abs(pdgId[i]) == 5 && std::abs(pdgId[mother_idx[i]]) == 6) {
            b_from_t = i;
	}
    }
    return b_idxs;
}
