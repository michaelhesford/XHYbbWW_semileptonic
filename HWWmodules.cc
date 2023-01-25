#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <stdio.h>

using namespace ROOT::VecOps;

RVec<int> PickDijets(RVec<float> pt, RVec<float> eta, RVec<float> jet_phi, RVec<float> mass, RVec<int> leptonIdx, RVec<float> electron_phi, RVec<float> muon_phi) {
// We are looking for a lead jet separated by the two sublead jets by at least 90 degrees
    int jet0Idx = -1;
    int jet1Idx = -1;
    // Since the vectors are ordered by pT, the first jet should roughly be our Higgs. 
    for (int ijet=0; ijet<pt.size(); ijet++) {
        if (leptonIdx[1] == -1) { // Our good lepton is an electron
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50 && hardware::DeltaPhi(jet_phi[ijet],electron_phi[leptonIdx[0]]) > M_PI/2) {
            // we've found our lead jet, break loop
                jet0Idx = ijet;
                break;
            } 
        }
        else { // good lepton is a muon
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50 && hardware::DeltaPhi(jet_phi[ijet],muon_phi[leptonIdx[1]]) > M_PI/2) {
                jet0Idx = ijet;
                break;
            }
        }
    }
    // if no lead jet found somehow, return 
    if (jet0Idx == -1) {
        return {-1, -1};
    }
    // now loop over the remaining jets, starting from the index of the lead jet
    for (int ijet=jet0Idx; ijet<pt.size(); ijet++) {
        if (pt[ijet] > 200 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 40 && hardware::DeltaPhi(jet_phi[jet0Idx],jet_phi[ijet]) > M_PI/2) {
            jet1Idx = ijet;
            break;
        }
    }
    return {jet0Idx,jet1Idx};
}

RVec<int> SignalElectrons(RVec<float> pt, RVec<float> eta, RVec<float> dxy, RVec<float> dz, RVec<float> sip3d,
RVec<float> iso_rel, RVec<float> hoe, RVec<float> eInvMinusPInv,RVec<float> deltaEtaSC, RVec<float> sieie) {
// Pick "good" electrons for preselection based on pT, eta, impact parameter, isolation requirements, and other variables
    int electronIdx = -1;
    for (int ie=0; ie<pt.size(); ie++) {
        if (pt[ie] > 30 && std::abs(eta[ie]) < 1.479 && dxy[ie] < 0.045 && dz[ie] < 0.2 && sip3d[ie] < ## && 
        iso_rel[ie] < 0.3 && hoe[ie] < ## &&  eInvMinusPInv[ie] < 0.184 && sieie[ie] < 0.0106 && deltaEtaSC < ##) {
            electronIdx = ie;
            break;
        }
    }
    return {electronIdx};
}

RVec<int> SignalMuons(RVec<float> pt, RVec<float> eta, RVec<float> dxy, RVec<float> dz, RVec<float> sip3d,RVec<float> iso_rel, RVec<float> mediumID, RVec<float> segmentComp, RVec<float> isGlobal) {
// Now pick signal-like muons
    int muonIdx = -1;
    for (int im=0; im<pt.size(); im++) {
        if isGlobal[im] == 1 {   // different requirments for segment compatibility for global vs tracker muons
            if (pt[im] > 27 && std::abs(eta[im]) < 2.4 && dxy[im] < 0.045 && dz[im] < 0.2 && sip3d[im] < ## && iso_rel[im] < 0.3 && mediumID[im] == 1 && segmentComp[im] > 0.303) {
                muonIdx = im;
                break;
            }
        }
        else {
            if (pt[im] > 27 && std::abs(eta[im]) < 2.4 && dxy[im] < 0.045 && dz[im] < 0.2 && sip3d[im] < ## && iso_rel[im] < 0.3 mediumID[im] == 1 && segmentComp[im] > 0.451) {
                muonIdx = im;                                                                                                  
                break; 
            }
        }
    }
    return {muonIdx};
}

RVec<int> LeptonIdx(RVec<int> electronIdx, RVec<int> muonIdx, RVec<float> electron_pt, RVec<float> muon_pt) {
// Decide whether an event is a lepton or muon type, collects lepton indices together
    int e_index == -1;
    int m_index == -1;
    if (electronIdx == -1) {
        m_index == muonIdx[0];
    }
    else if (muonIdx == -1) {
        e_index == electronIdx[0];
    }
    else {
        if (electron_pt[electronIdx[0]] > muon_pt[muonIdx[0]]) {
            e_index == electronIdx[0];
        }
        else {
            m_index == muonIdx[0];
        }
    }
    return {e_index,m_index};
 






