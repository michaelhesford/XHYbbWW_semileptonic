#include <iostream>
#include <TH1.h>
#include <string>
#include <TEfficiency.h>
#include <ROOT/RVec.hxx>
#include "correction.h"
#include <cstdlib>
#include <stdio.h>

using namespace ROOT::VecOps;

RVec<float> make_muo_sf_vec(std::shared_ptr<const correction::Correction> sf_set, float abseta, float pt, int type, std::string cat) {
    // cat: "high" or "low" (specifies which triggers we are applying, i.e. which leptons to consider for SF application)
   
    if (type == 0) {
	// Electron event, return no SF
	return {1.0,1.0,1.0};
    }
    if (cat == "high" && (pt < 55 || pt > 200)) { // SF's don't exist for pt > 200
	return {1.0,1.0,1.0};
    }
    else if (cat == "low" && (pt < 29 || pt >= 55)) {
	return {1.0,1.0,1.0};
    }
    else {
	float nom = (float)sf_set->evaluate({abseta, pt, "nominal"});
        float stat = (float)sf_set->evaluate({abseta, pt, "stat"});
	float syst = (float)sf_set->evaluate({abseta, pt, "syst"});
        float uncert = stat + syst;
        return {nom, nom+uncert, nom-uncert};
    }
}

RVec<float> make_egm_sf_vec(int type, float eta, float pt, float sf_set_nom[3][5], float sf_set_up[3][5], float sf_set_down[3][5]) {

    if (type == 1) {
        // Muon event, return no SF
	return {1.0,1.0,1.0};
    }
    else if (pt < 30 || pt > 2000) {
        return {1.0,1.0,1.0};
    }
    else {
        // First determine sf bin
	int pt_bin;
	int eta_bin;

	if (pt >= 30 && pt < 120) {pt_bin = 0;}
	else if (pt >= 120 && pt < 200) {pt_bin = 1;}
	else if (pt >= 200 && pt <= 2000) {pt_bin = 2;}

	if (eta >= -2.5 &&eta < -1.5) {eta_bin = 0;}
	else if (eta >= -1.5 && eta < -0.5) {eta_bin = 1;}
	else if (eta >= -0.5 && eta < 0.5) {eta_bin = 2;}
	else if (eta >= 0.5 && eta < 1.5) {eta_bin = 3;}
	else if (eta >= 1.5 && eta <= 2.5)  {eta_bin = 4;}

	float sf_nom = sf_set_nom[pt_bin][eta_bin];
	float sf_up = sf_set_up[pt_bin][eta_bin];
	float sf_down = sf_set_down[pt_bin][eta_bin];

	return {sf_nom,sf_up,sf_down};
    }
}
