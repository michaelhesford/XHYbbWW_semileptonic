#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <stdio.h>
#include <cmath>

RVec<int> PickDijets(RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    // We are looking for a lead jet separated by a sublead jet by at least 90 degrees
    int jet0Idx = -1;
    int jet1Idx = -1;
    // Since the vectors are ordered by pT, the first jet should roughly be our Higgs. 
    for (int ijet=0; ijet<eta.size(); ijet++) {
        if (std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
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
    for (int ijet=jet0Idx; ijet<eta.size(); ijet++) {
        if (std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50 && std::abs(hardware::DeltaPhi(phi[jet0Idx],phi[ijet])) > M_PI/2) {
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

int SignalLepton(RVec<float> lep_pt, RVec<float> lep_eta, RVec<float> lep_phi, RVec<float> lep_iso, float H_phi, float W_phi) {
// function for lepton preselction (will expand/make separate for electrons/muons in future)
    int leptonIdx = -1;
    for (int il=0; il<lep_pt.size(); il++) {
        if (lep_pt[il] > 25 && std::abs(lep_eta[il]) < 2.4 && lep_iso[il] < 0.2 && std::abs(hardware::DeltaPhi(lep_phi[il],W_phi)) < 2.5 && std::abs(hardware::DeltaPhi(lep_phi[il],W_phi)) > 1 && std::abs(hardware::DeltaPhi(lep_phi[il],H_phi)) > M_PI/2) {
            leptonIdx = il;
            break;
        }
    }
    return leptonIdx;
}

RVec<int> LeptonIdx(int electronIdx, int muonIdx, RVec<float> electron_pt, RVec<float> muon_pt) {
// Decide whether an event is a lepton or muon type, collects lepton indices together
    int e_index = -1;
    int m_index = -1;
    if (electronIdx == -1) {
        m_index = muonIdx;
    }
    else if (muonIdx == -1) {
        e_index = electronIdx;
    }
    else {
        if (electron_pt[electronIdx] > muon_pt[muonIdx]) {
            e_index = electronIdx;
        }
        else {
            m_index = muonIdx;
        }
    }
    return {e_index,m_index};
}

RVec<int> JetMasses(RVec<float> mass) { // Selects jets with pt > 50 GeV
    int jet0Idx = -1;
    int jet1Idx = -1;
    int jet2Idx = -1;
    for (int ijet = 0; ijet<mass.size(); ijet++) {   
        if (mass[ijet] > 50) {
            jet0Idx = ijet;
            break;
        }
    }
    if (jet0Idx == -1) {
        return {-1,-1,-1};
    }
    for (int ijet = jet0Idx+1; ijet<mass.size(); ijet++) {
        if (mass[ijet] > 50) {
            jet1Idx = ijet;
            break;
        }
    }
    if (jet1Idx == -1) {
        return {jet0Idx,-1,-1};
    }
    for (int ijet = jet1Idx+1; ijet<mass.size(); ijet++) {
        if (mass[ijet] > 50) {
            jet2Idx = ijet;
            break;
        }
    }
    if (jet2Idx == -1) {
        return {jet0Idx,jet1Idx,-1};
    }
    else {
        return {jet0Idx,jet1Idx,jet2Idx};
    }
}

RVec<int> FindMothersPdgId(RVec<int> genpart_id, RVec<int> selected_genpart_mother_indices){

    std::size_t Ni = selected_genpart_mother_indices.size();
    RVec<int> mother_pdgids(Ni);
    for(std::size_t i=0; i<Ni; i++) {
        mother_pdgids[i] = genpart_id[selected_genpart_mother_indices[i]];
    }
    return mother_pdgids;

}

int HighestPtIdx(RVec<float> lep_pt){
    int idx = -1;
    float pt = 0;
    for (int id = 0; id<lep_pt.size(); id++) {
        if (lep_pt[id] > pt) {
            idx = id;
            pt = lep_pt[id];
        }
    }
    return idx;
}

float TransverseMass(float MET_pt, float obj_pt, float MET_phi, float obj_phi) {
    return sqrt(2.0*MET_pt*obj_pt*(1-cos(hardware::DeltaPhi(MET_phi,obj_phi))));
}

/*
int GenSignalLepton(RVec<float> lep_pt, RVec<float> lep_eta, RVec<float> lep_phi, float H_phi, float W_phi) {
//function for lepton preselction (will expand/make separate for electrons/muons in future)
    RVec<int> leptonIdx = -1;
    for (int il=0; il<lep_pt.size(); il++) {
        if (lep_pt[il] > 25 && std::abs(lep_eta[il]) < 2.4 && std::abs(hardware::DeltaPhi(lep_phi[il],W_phi)) < 2.5 && std::abs(hardware::DeltaPhi(lep_phi[il],W_phi)) > 1 && std::abs(hardware::DeltaPhi(lep_phi[il],H_phi)) > M_PI/2) {
            leptonIdx = il
            break;
        }
    }
     return leptonIdx;
}
*/
