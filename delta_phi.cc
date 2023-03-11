#include "TIMBER/Framework/include/common.h"
#include "ROOT/RVec.hxx"
#include <stdio.h>

using namespace ROOT::VecOps;

float DeltaPhi(RVec<float> phi1, RVec<float> phi2, RVec<int> idxs){
    // Get the indices
    int idx1 = idxs[0];
    int idx2 = idxs[1];
    // Calculate delta phi between two objects
    float delta_phi = hardware::DeltaPhi(phi1[idx1],phi2[idx2]);
    return delta_phi;
}

float DeltaPhiSingular(RVec<float> phi1, float phi2, int idx){
    // Calculate delta phi between two objects
    float delta_phi = hardware::DeltaPhi(phi1[idx],phi2);
    return delta_phi;
}

float DeltaPhiSingular2(float phi1, float phi2) {
    float delta_phi = hardware::DeltaPhi(phi1,phi2);
    return delta_phi;
}
