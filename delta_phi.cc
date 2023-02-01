#include "TIMBER/Framework/include/common.h"
#include "ROOT/RVec.hxx"
#include <stdio.h>

using namespace ROOT::VecOps;

RVec<float> DeltaPhi(RVec<float> phi1, RVec<float> phi2){
    // Calculate delta phi between two objects
    float delta_phi = hardware::DeltaPhi(phi1[0],phi2[0]);
    return {delta_phi};
}
