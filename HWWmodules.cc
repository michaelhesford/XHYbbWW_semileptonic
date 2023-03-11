#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <stdio.h>
#include <cmath>

using namespace ROOT::VecOps;

RVec<int> PickDijets(RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> jet_mass, RVec<int> leptonIdx, RVec<float> electron_phi, RVec<float> muon_phi) {
    // Looking for the closest jet to our lepton (the W candidate) as well as a jet separated by both by Pi/2 -- update minimum delta phi requirement for W/lepton
    int jet0Idx = -1;
    int jet1Idx = -1;
    RVec<float> lepton_phi = muon_phi; // lets assume its a muon...
    int lepId = leptonIdx[1];
    if (lepId == -1) { // Our good lepton is an electron
        lepton_phi = electron_phi;
        lepId = leptonIdx[0];
    }
    for (int ijet=0; ijet<jet_mass.size(); ijet++) {
        if (jet0Idx == -1) { // no jet selected yet
            if (jet_mass[ijet] > 50 && jet_eta[ijet] < 2.4 && std::abs(hardware::DeltaPhi(jet_phi[ijet],lepton_phi[lepId])) < M_PI/2) {
                jet0Idx = ijet;
            }
        }
        else {
            if (jet_mass[ijet] > 50 && jet_eta[ijet] < 2.4 && std::abs(hardware::DeltaPhi(jet_phi[ijet],lepton_phi[lepId])) < std::abs(hardware::DeltaPhi(jet_phi[jet0Idx],lepton_phi[lepId]))) {
                jet0Idx = ijet;
            }
        }
    }
    
    for (int ijet=0; ijet<jet_mass.size(); ijet++) {
        if (jet_mass[ijet] > 50 && jet_eta[ijet] < 2.4 && std::abs(hardware::DeltaPhi(jet_phi[ijet],jet_phi[jet0Idx])) > M_PI/2) {
            jet1Idx = ijet;
        }
    }
    return {jet0Idx,jet1Idx};
}

int LeptonPre(RVec<float> pt, RVec<float> eta, RVec<float> iso) {
// function for lepton preselction (will expand/make separate for electrons/muons in future)
    int leptonIdx = -1;
    for (int il=0; il<pt.size(); il++) {
        if (pt[il] > 10 && std::abs(eta[il]) < 2.5 && iso[il] < 0.3) {
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

int GetClosestJet(RVec<int> leptonIdx, RVec<float> electron_phi, RVec<float> muon_phi, RVec<float> jet_phi) {
    RVec<float> lepton_phi = muon_phi; // assume lepton a muon...
    int lepId = leptonIdx[1];
    if (lepId == -1) { // Our good lepton is an electron
        lepton_phi = electron_phi;
        lepId = leptonIdx[0];
    }
    int jet_index = 0;
    for (int ij=0; ij<jet_phi.size(); ij++) {
        if (abs(hardware::DeltaPhi(lepton_phi[lepId],jet_phi[ij])) < abs(hardware::DeltaPhi(lepton_phi[lepId],jet_phi[jet_index]))) {
            jet_index = ij;
        }
    }
    return jet_index;
} 
   
float GetMETmass(float pt, float Et) {
    float c = 2.998*pow(10,8);
    float exp = 0.5;
    float mass = pow((Et*Et-pt*pt*pow(c,2))/pow(c,4),exp);
    return mass;
}

float GetMETeta() {
    float eta = 0;
    return eta;
}

float TransverseMass(float lep_pt, float MET) {
    float mass_tran = pow(2*lep_pt*MET,0.5);
    return mass_tran;
}








