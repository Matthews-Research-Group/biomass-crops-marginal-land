#ifndef MISCANTHUS_SENESCENCE_H
#define MISCANTHUS_SENESCENCE_H

#include "../modules.h"
#include "../state_map.h"
#include <R.h>
/**
 * @class miscanthus_senescence
 * 
 * @brief Calculates leaf senescence based on: (1) ageing of leaves (2) maximum lai (3)frost.
 * maximum of the three senescnece is taken to represent current leaf senescence
 * This module is based on leaf senescene in APSIM (https://www.apsim.info/documentation/model-documentation/crop-module-documentation/plant/)
 * stem, root, and rhizome senescene are set to zero (hard coded )and yet to be implemented. 
 * for water stress see Figure 4 of https://link.springer.com/article/10.1007/s11258-013-0228-4.
 * A constant value of 0.5 for soil porosity is used, which should be updated from select soil (missing porosity)
 */

class miscanthus_senescence : public DerivModule {
	public:
	  miscanthus_senescence(const std::unordered_map<std::string, double>* input_parameters, std::unordered_map<std::string, double>* output_parameters) :
			// Define basic module properties by passing its name to its parent class
            DerivModule("miscanthus_senescence"),
			// Get pointers to input parameters
            Leaf_ip(get_ip(input_parameters, "Leaf")),
            Stem_ip(get_ip(input_parameters, "Stem")),
            Root_ip(get_ip(input_parameters, "Root")),
            Rhizome_ip(get_ip(input_parameters, "Rhizome")),
			      lai_max_ip(get_ip(input_parameters, "lai_max")),
			      lai_ip(get_ip(input_parameters, "lai")),
			      Sp_ip(get_ip(input_parameters, "Sp")),
			      TTc_ip(get_ip(input_parameters, "TTc")),
			      delta_TT_ip(get_ip(input_parameters, "delta_TT")),
            remobilization_fraction_leaf_to_rhizome_ip(get_ip(input_parameters, "remobilization_fraction_leaf_to_rhizome")),
            remobilization_fraction_stem_to_rhizome_ip(get_ip(input_parameters, "remobilization_fraction_stem_to_rhizome")),
            remobilization_fraction_root_to_rhizome_ip(get_ip(input_parameters, "remobilization_fraction_root_to_rhizome")),
            remobilization_fraction_rhizome_to_rhizome_ip(get_ip(input_parameters, "remobilization_fraction_rhizome_to_rhizome")),
            leaf_turnover_rate_ip(get_ip(input_parameters, "leaf_turnover_rate")),
            stem_turnover_rate_ip(get_ip(input_parameters, "stem_turnover_rate")),
            root_turnover_rate_ip(get_ip(input_parameters, "root_turnover_rate")),
            rhizome_turnover_rate_ip(get_ip(input_parameters, "rhizome_turnover_rate")),
            TTc_leafsenescence_threshold_ip(get_ip(input_parameters, "TTc_leafsenescence_threshold")),
            TTc_stemsenescence_threshold_ip(get_ip(input_parameters, "TTc_stemsenescence_threshold")),
            TTc_rootsenescence_threshold_ip(get_ip(input_parameters, "TTc_rootsenescence_threshold")),
            TTc_rhizomesenescence_threshold_ip(get_ip(input_parameters, "TTc_rhizomesenescence_threshold")),
            Tfrostlow_ip(get_ip(input_parameters, "Tfrostlow")),
            Tfrosthigh_ip(get_ip(input_parameters, "Tfrosthigh")),
            temp_ip(get_ip(input_parameters, "temp")),
            StomataWS_ip(get_ip(input_parameters, "StomataWS")),
            soil_field_capacity_ip(get_ip(input_parameters, "soil_field_capacity")),
            phi_waterstress_induced_leafsenescence_ip(get_ip(input_parameters,"phi_waterstress_induced_leafsenescence")),
            sene_factor_when_sws_eq_0_ip(get_ip(input_parameters,"sene_factor_when_sws_eq_0")),
            sene_factor_when_sws_eq_1_ip(get_ip(input_parameters,"sene_factor_when_sws_eq_1")),
            
			// Get pointers to output parameters
			Leaf_op(get_op(output_parameters, "Leaf")),
			LeafLitter_op(get_op(output_parameters, "LeafLitter")),
			Stem_op(get_op(output_parameters, "Stem")),
			StemLitter_op(get_op(output_parameters, "StemLitter")),
			Root_op(get_op(output_parameters, "Root")),
			RootLitter_op(get_op(output_parameters, "RootLitter")),
			Rhizome_op(get_op(output_parameters, "Rhizome")),
			RhizomeLitter_op(get_op(output_parameters, "RhizomeLitter")),
	    Grain_op(get_op(output_parameters, "Grain"))
		{}
		static std::vector<std::string> get_inputs();
		static std::vector<std::string> get_outputs();
	private:
		// Pointers to input parameters
        const double* Leaf_ip;
        const double* Stem_ip;
        const double* Root_ip;
        const double* Rhizome_ip;
        const double* lai_max_ip;
        const double* lai_ip;
        const double* Sp_ip;
        const double* TTc_ip;
        const double* delta_TT_ip;
        const double* remobilization_fraction_leaf_to_rhizome_ip;
        const double* remobilization_fraction_stem_to_rhizome_ip;
        const double* remobilization_fraction_root_to_rhizome_ip;
        const double* remobilization_fraction_rhizome_to_rhizome_ip;
        const double* leaf_turnover_rate_ip;
        const double* stem_turnover_rate_ip;
        const double* root_turnover_rate_ip;
        const double* TTc_leafsenescence_threshold_ip;
        const double* TTc_stemsenescence_threshold_ip;
        const double* TTc_rootsenescence_threshold_ip;
        const double* TTc_rhizomesenescence_threshold_ip;
        const double* rhizome_turnover_rate_ip;
        const double* Tfrostlow_ip;
        const double* Tfrosthigh_ip;
        const double* temp_ip;
        const double* soil_field_capacity_ip;
        const double* StomataWS_ip;
        const double* phi_waterstress_induced_leafsenescence_ip;
        const double* sene_factor_when_sws_eq_0_ip;
        const double* sene_factor_when_sws_eq_1_ip;

		// Pointers to output parameters
		double* Leaf_op;
		double* LeafLitter_op;
		double* Stem_op;
		double* StemLitter_op;
		double* Root_op;
		double* RootLitter_op;
		double* Rhizome_op;
		double* RhizomeLitter_op;
		double* Grain_op;
		// Main operation
		void do_operation() const;
};

std::vector<std::string> miscanthus_senescence::get_inputs() {
	return {
		"Leaf", "Stem", "Root", "Rhizome",
		"lai_max","lai","Sp","TTc", "delta_TT","remobilization_fraction_leaf_to_rhizome", 
		"remobilization_fraction_stem_to_rhizome", 
		"remobilization_fraction_root_to_rhizome", 
		"remobilization_fraction_rhizome_to_rhizome", 
		"leaf_turnover_rate",
		"stem_turnover_rate",
		"root_turnover_rate",
		"rhizome_turnover_rate",
		"TTc_leafsenescence_threshold",
		"TTc_stemsenescence_threshold",
		"TTc_rootsenescence_threshold",
		"TTc_rhizomesenescence_threshold",
		"Tfrostlow", "Tfrosthigh","temp",
		"soil_field_capacity", "StomataWS","phi_waterstress_induced_leafsenescence","sene_factor_when_sws_eq_0",
		 "sene_factor_when_sws_eq_1"
	};
}

std::vector<std::string> miscanthus_senescence::get_outputs() {
	return {
		"Leaf", "LeafLitter","Stem", "StemLitter",
		"Root", "RootLitter", "Rhizome", "RhizomeLitter","Grain"
	};
}

void miscanthus_senescence::do_operation() const {
  
    double Leaf = *Leaf_ip; //Current leaf biomass, Mg/ha
    double Stem = *Stem_ip; // current stem biomass, Mg/ha
    double Root = *Root_ip; // current root biomass, Mg/ha
    double Rhizome = *Rhizome_ip; // current rhizome biomass, Mg/ha
    double lai_max = *lai_max_ip; //maximum permissible leaf area index
    double lai = *lai_ip; //maximum permissible leaf area index
    double Sp = *Sp_ip; //maximum permissible leaf area index
    double TTc = *TTc_ip;
    double delta_TT = *delta_TT_ip;
    double remobilization_fraction_leaf_to_rhizome = *remobilization_fraction_leaf_to_rhizome_ip; // remobilization fraction from leaf to rhizome upon senescence
    double remobilization_fraction_stem_to_rhizome = *remobilization_fraction_stem_to_rhizome_ip;
    double remobilization_fraction_root_to_rhizome = *remobilization_fraction_root_to_rhizome_ip;
    double remobilization_fraction_rhizome_to_rhizome = *remobilization_fraction_rhizome_to_rhizome_ip;
    double leaf_turnover_rate = *leaf_turnover_rate_ip;
    double stem_turnover_rate = *stem_turnover_rate_ip;
    double root_turnover_rate = *root_turnover_rate_ip;
    double rhizome_turnover_rate = *rhizome_turnover_rate_ip;
    double TTc_leafsenescence_threshold = *TTc_leafsenescence_threshold_ip;
    double TTc_stemsenescence_threshold = *TTc_stemsenescence_threshold_ip;
    double TTc_rootsenescence_threshold = *TTc_rootsenescence_threshold_ip;
    double TTc_rhizomesenescence_threshold = *TTc_rhizomesenescence_threshold_ip;
    double Tfrostlow = *Tfrostlow_ip;
    double Tfrosthigh = *Tfrosthigh_ip;
    double temp = *temp_ip;
    double soil_field_capacity = *soil_field_capacity_ip;
    double StomataWS = *StomataWS_ip;
    double phi_waterstress_induced_leafsenescence = *phi_waterstress_induced_leafsenescence_ip;
    double sene_factor_when_sws_eq_0 = *sene_factor_when_sws_eq_0_ip;
    double sene_factor_when_sws_eq_1 = *sene_factor_when_sws_eq_1_ip;

	  double senescence_leaf_frost = std::max(0.0, std::min(1.0,(Tfrosthigh - temp) / (Tfrosthigh - Tfrostlow))) * Leaf*0.005; // Amount of leaf senesced due to frost, Mg/ha
	  double senescence_leaf_lightcompetition = (lai > lai_max) ? (lai - lai_max)/Sp : 0.0 ; // Amount of leaf senesced due to light competition, Mg/ha
	  double senescence_leaf_ageing = (TTc > TTc_leafsenescence_threshold) ? Leaf*leaf_turnover_rate*delta_TT : 0.0;
   
    double waterstress_factor_on_senescence = (sene_factor_when_sws_eq_1 - sene_factor_when_sws_eq_0) * StomataWS + sene_factor_when_sws_eq_0;
    
    
    double senescence_leaf;
    if (senescence_leaf_lightcompetition > senescence_leaf_frost){
      senescence_leaf = (senescence_leaf_lightcompetition > senescence_leaf_ageing) ? senescence_leaf_lightcompetition : senescence_leaf_ageing;
    }else{
        senescence_leaf = (senescence_leaf_frost > senescence_leaf_ageing) ? senescence_leaf_frost : senescence_leaf_ageing;
    } 
    
    double senescence_root = (TTc > TTc_rootsenescence_threshold) ? Root * root_turnover_rate * waterstress_factor_on_senescence :0.0 ;
    double senescence_stem = (TTc > TTc_stemsenescence_threshold) ? Stem * stem_turnover_rate: 0.0 ; // Amount of stem senesced in Mg/ha, update it if needed
    double senescence_rhizome = (TTc > TTc_rhizomesenescence_threshold) ? Rhizome*rhizome_turnover_rate * waterstress_factor_on_senescence : 0.0; // Amount of rhizome senesced in Mg/ha
    
    
    double dLeaf = -senescence_leaf; 
    double dLeafLitter = senescence_leaf * (1 - remobilization_fraction_leaf_to_rhizome); // part of senescced leaf that does not remobilize become leaf litter
    
    double dStem = -senescence_stem; // No stem senescence for now. Need to update it later
    double dStemLitter = senescence_stem * (1 - remobilization_fraction_stem_to_rhizome); // part of senescced stem that does not remobilize become stem litter
    
    double dRoot = -senescence_root; //no root senesnce for now. Need to update it later
    double dRootLitter = senescence_root* (1 - remobilization_fraction_root_to_rhizome); // part of senescced roots that does not remobilize become leaf litter
    
    double dRhizome = -senescence_rhizome + 
                      senescence_leaf * remobilization_fraction_leaf_to_rhizome +
                      senescence_stem * remobilization_fraction_stem_to_rhizome +
                      senescence_root * remobilization_fraction_root_to_rhizome +
                      senescence_rhizome * remobilization_fraction_rhizome_to_rhizome;
    
    double dRhizomeLitter = senescence_rhizome * (1-remobilization_fraction_rhizome_to_rhizome);
    
    double dGrain = 0.0 ; // Grain is zero. I am avoiding using kGrain here as I removed this from the list of inputs

    
    update(Leaf_op, dLeaf); // Mg/ha
    update(Stem_op, dStem); // Mg/ha
    update(Root_op, dRoot); // Mg/ha
    update(Rhizome_op, dRhizome); // Mg/ha
    update(Grain_op, dGrain); // Mg/ha
    
    update(LeafLitter_op, dLeafLitter); // Mg/ha
    update(StemLitter_op, dStemLitter); // Mg/ha
    update(RootLitter_op, dRootLitter); // Mg/ha
    update(RhizomeLitter_op, dRhizomeLitter); // Mg/ha
    
}

#endif
