#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <stdio.h>

using namespace ROOT::VecOps;

RVec<int> PickDijets(RVec<float> pt, RVec<float> eta, RVec<float> jet_phi, RVec<float> mass) {
// We are looking for a lead jet separated by the two sublead jets by at least 90 degrees
    int jet0Idx = -1;
    int jet1Idx = -1;
    // Since the vectors are ordered by pT, the first jet should roughly be our Higgs. 
    for (int ijet=0; ijet<pt.size(); ijet++) {
        if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
        // we've found our lead jet, break loop
            jet0Idx = ijet;
            break;
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

RVec<int> PickIsolatedLeptons(RVec<float> jet_phi, RVec<float> phi, RVec<float> pt, RVec<float> eta, RVec<float> dxy, Vec<float> dz, RVec<float> iso_rel) {
// Pick "good" leptons for preselection based on pT, eta, dxyz, and isolation requirements
    int leptonIdx = -1;
    for (int il=0; il<phi.size(); il++) {
        if (pt[il] > 10 && std::abs(eta[il]) < 2.4 && dxy[il] < 0.045 && dz[il] < 0.2 && iso_rel < 0.3 && hardware::DeltaPhi(phi[il],jet_phi[0]) > M_PI/2) {
            leptonIdx = il;
            break;
        }
    }
    return {leptonIdx};
}
