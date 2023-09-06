#ifndef LOGISTIC_PARTITIONING_WITH_RHIZOME_REMOBILIZATION_H
#define LOGISTIC_PARTITIONING_WITH_RHIZOME_REMOBILIZATION_H

#include <cmath>
#include "../modules.h"

/**
 * \brief
 * References:
 * Osborne, T. et al. 2015. “JULES-Crop: A Parametrisation of Crops in the Joint UK Land Environment Simulator.” Geoscientific Model Development 8(4): 1139–55. https://doi.org/10.5194/gmd-8-1139-2015
 */

class logistic_partitioning_with_rhizome_remobilization : public SteadyModule {
	public:
	  logistic_partitioning_with_rhizome_remobilization(const std::unordered_map<std::string, double>* input_parameters, std::unordered_map<std::string, double>* output_parameters) :
			// Define basic module properties by passing its name to its parent class
			SteadyModule("logistic_partitioning_with_rhizome_remobilization"),
			// Get pointers to input parameters
            DVI_ip(get_ip(input_parameters,"DVI")),
            alphaRoot_ip(get_ip(input_parameters, "alphaRoot")),
            alphaStem_ip(get_ip(input_parameters, "alphaStem")),
            alphaLeaf_ip(get_ip(input_parameters, "alphaLeaf")),
            betaRoot_ip(get_ip(input_parameters, "betaRoot")),
            betaStem_ip(get_ip(input_parameters, "betaStem")),
            betaLeaf_ip(get_ip(input_parameters, "betaLeaf")),
            kRhizome_emr_ip(get_ip(input_parameters, "kRhizome_emr")),
            kStem_emr_ip(get_ip(input_parameters, "kStem_emr")),
            kLeaf_emr_ip(get_ip(input_parameters, "kLeaf_emr")),
            kRoot_emr_ip(get_ip(input_parameters, "kRoot_emr")),
           
           
    
			// Get pointers to output parameters
			kStem_op(get_op(output_parameters, "kStem")),
			kLeaf_op(get_op(output_parameters, "kLeaf")),
			kRoot_op(get_op(output_parameters, "kRoot")),
			kRhizome_op(get_op(output_parameters, "kRhizome")),
	    kGrain_op(get_op(output_parameters, "kGrain"))

    {}
		static std::vector<std::string> get_inputs();
		static std::vector<std::string> get_outputs();
	private:
		// Pointers to input parameters
        const double* DVI_ip;
		const double* alphaRoot_ip;
        const double* alphaStem_ip;
		const double* alphaLeaf_ip;
		const double* betaRoot_ip;
		const double* betaStem_ip;
		const double* betaLeaf_ip;
        const double* kRhizome_emr_ip;
        const double* kStem_emr_ip;
        const double* kLeaf_emr_ip;
        const double* kRoot_emr_ip;
		
		// Pointers to output parameters
		double* kStem_op;
		double* kLeaf_op;
		double* kRoot_op;
		double* kRhizome_op;
		double* kGrain_op;
		// Main operation
		void do_operation() const;
};

std::vector<std::string> logistic_partitioning_with_rhizome_remobilization::get_inputs() {
	return {
		"DVI", "alphaRoot", "alphaStem", "alphaLeaf",
		"betaRoot", "betaStem", "betaLeaf", "kRhizome_emr",
		"kStem_emr","kLeaf_emr","kRoot_emr"
	};
}

std::vector<std::string> logistic_partitioning_with_rhizome_remobilization::get_outputs() {
	return {
		"kStem", "kLeaf", "kRoot", "kRhizome", "kGrain"
	};
}

void logistic_partitioning_with_rhizome_remobilization::do_operation() const {
    
   const double DVI = *DVI_ip; // dimensionless; Development Index
	 double kStem, kLeaf, kRoot, kRhizome, kGrain, kDenom;
    
    // Determine partitioning coefficients using multinomial logistic equations from Osborne et al., 2015 JULES-crop doi:10.5194/gmd-8-1139-2015
   
// IF DVI <0 tables based partitioning coefficients are used
// IF DVI >0  then root, stem, leaf,  and rhizome ---all the partitioning coefficients are calculated using logistic function
    if (DVI<0) {
//      kRoot = *kRoot_emr_ip; //dimensionless
      kStem = *kStem_emr_ip; //dimensionless
      kLeaf = *kLeaf_emr_ip; //dimensionless
      kRhizome = *kRhizome_emr_ip; //dimensionless
      kRoot = 1.0 - kLeaf - kStem;
    } else {
      kDenom = exp(*alphaRoot_ip + *betaRoot_ip * DVI) + exp(*alphaLeaf_ip + *betaLeaf_ip * DVI) + exp(*alphaStem_ip + *betaStem_ip * DVI)  + 1.0; // denominator term for kRoot, kStem, kLeaf, and kRhizome, dimensionless
      kRoot = exp(*alphaRoot_ip + *betaRoot_ip * DVI) / kDenom; // dimensionless
      kStem = exp(*alphaStem_ip + *betaStem_ip * DVI) / kDenom; // dimensionless
      kLeaf = exp(*alphaLeaf_ip + *betaLeaf_ip * DVI) / kDenom; // dimensionless
      kRhizome = 1.0/ kDenom; // dimensionless
    }
    kGrain= 0; // kGrain is always zero. It is being used because it is needed by other modules
    // Update the output parameters
    update(kStem_op, kStem); //dimensionless
    update(kLeaf_op, kLeaf); //dimensionless
    update(kRoot_op, kRoot); //dimensionless
    update(kRhizome_op, kRhizome); //dimensionless
}

#endif
