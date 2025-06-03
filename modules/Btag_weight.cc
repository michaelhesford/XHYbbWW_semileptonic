#include <iostream>
#include <TH1.h>
#include <string>
#include <TEfficiency.h>
#include <ROOT/RVec.hxx>
#include "correction.h"

using namespace ROOT::VecOps;

float evaluate_corr(std::shared_ptr<const correction::Correction> bc_sf, std::shared_ptr<const correction::Correction> incl_sf, TEfficiency* eff_b, TEfficiency* eff_c, TEfficiency* eff_light, RVec<ROOT::Math::PtEtaPhiMVector> Jet_vec, RVec<int> Jet_hadronFlavour, RVec<float> Jet_btagDisc, std::string wp, float wp_val) {
    // Calculate nominal event reweight - handles scale factors for b/c/light jets in one function
    float val_tag   = 1.0;
    float val_notag = 1.0;
    for (size_t ijet = 0; ijet < Jet_vec.size(); ijet++) {
        // find MC efficiency/scale factor
	float efficiency;
	float sf;
	float pt = Jet_vec[ijet].Pt();
	float abseta = std::abs(Jet_vec[ijet].Eta());
        if (Jet_hadronFlavour[ijet] == 5) { // b-jet
	    int globalbin = eff_b->FindFixBin(pt,abseta);
	    efficiency = eff_b->GetEfficiency(globalbin); 
    	    sf = (float)bc_sf->evaluate({"central", wp.c_str(), Jet_hadronFlavour[ijet], abseta, pt});
        }
	else if (Jet_hadronFlavour[ijet] == 4) { // c-jet
            int globalbin = eff_c->FindFixBin(pt,abseta);
            efficiency = eff_c->GetEfficiency(globalbin);
            sf = (float)bc_sf->evaluate({"central", wp.c_str(), Jet_hadronFlavour[ijet], abseta, pt});
        } 
	else { // light jet
	    int globalbin = eff_light->FindFixBin(pt,abseta);
            efficiency = eff_light->GetEfficiency(globalbin);
            sf = (float)incl_sf->evaluate({"central", wp.c_str(), Jet_hadronFlavour[ijet], abseta, pt});
	}

        // Calculate contribution to event weight
        if (Jet_btagDisc[ijet] > wp_val) {
	    float valTag = (sf * efficiency) / (efficiency);
            val_tag *= valTag;
        }
	else {
            float valNoTag = (1. - sf*efficiency) / (1. - efficiency);
            val_notag *= valNoTag;
	}
    }
    return val_tag * val_notag; // overall event weight
}

RVec<float> evaluate_uncert(std::shared_ptr<const correction::Correction> sf_set, std::string corr, std::vector<int> hadron_flavours, std::vector<TEfficiency*> eff_maps, RVec<ROOT::Math::PtEtaPhiMVector> Jet_vec, RVec<int> Jet_hadronFlavour, RVec<float> Jet_btagDisc, std::string wp, float wp_val) {
    float val_tag_up   = 1.0;
    float val_notag_up = 1.0;
    float val_tag_dn   = 1.0;
    float val_notag_dn = 1.0;
    // Loop over all jets 
    for (size_t ijet = 0; ijet < Jet_vec.size(); ijet++) {
	// Check if jet is one of the target flavors
	auto it = std::find(hadron_flavours.begin(),hadron_flavours.end(),Jet_hadronFlavour[ijet]);
	if (it != hadron_flavours.end()) {
	    int fidx = it - hadron_flavours.begin();
            float pt = Jet_vec[ijet].Pt();
            float abseta = std::abs(Jet_vec[ijet].Eta());
	    // get efficiency/scale factors
	    int globalbin = eff_maps[fidx]->FindFixBin(pt,abseta);
            float efficiency = eff_maps[fidx]->GetEfficiency(globalbin);
            float sf_up = (float)sf_set->evaluate({"up_"+corr, wp.c_str(), Jet_hadronFlavour[ijet], abseta, pt});
            float sf_down = (float)sf_set->evaluate({"down_"+corr, wp.c_str(), Jet_hadronFlavour[ijet], abseta, pt});
        
	    // calculate contribution to event weight
	    if (Jet_btagDisc[ijet] > wp_val) {
                float valUpTag = (sf_up * efficiency) / (efficiency);
                float valDnTag = (sf_down * efficiency) / (efficiency);
                val_tag_up *= valUpTag;
                val_tag_dn *= valDnTag;
            }
            else {
                float valUpNoTag = (1. - (sf_up * efficiency)) / (1. - efficiency);
                float valDnNoTag = (1. - (sf_down * efficiency)) / (1. - efficiency);
                val_notag_up *= valUpNoTag;
                val_notag_dn *= valDnNoTag;
            }
	}
    }
    float weight_up = val_tag_up * val_notag_up;
    float weight_dn = val_tag_dn * val_notag_dn;
    return {weight_up, weight_dn};
}

