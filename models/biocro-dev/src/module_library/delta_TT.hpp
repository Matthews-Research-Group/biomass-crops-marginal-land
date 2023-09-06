#ifndef DELTA_TT_H
#define DELTA_TT_H

#include "../modules.h"

/**
 * @class DELTA TT
 * 
 * @brief Calculates the rate of thermal time accumulation. Thermal time is measured in
 * degree-days, so the resulting derivative has units of degree-days per hour. This model
 * includes an optimum temperature and a high temperature cutoff, as described in many
 * places, e.g. Ruiz‐Vera et al., Plant, Cell & Environment 41, 2806–2820 (2018).
 */
class delta_TT : public SteadyModule
{
   public:
    delta_TT(
        const std::unordered_map<std::string, double>* input_parameters,
        std::unordered_map<std::string, double>* output_parameters)
        :  // Define basic module properties by passing its name to its parent class
           SteadyModule("delta_TT"),
          // Get pointers to input parameters
          temp_ip(get_ip(input_parameters, "temp")),
          tbase_ip(get_ip(input_parameters, "tbase")),
          topt_lower_ip(get_ip(input_parameters, "topt_lower")),
          topt_upper_ip(get_ip(input_parameters, "topt_upper")),
          tmax_ip(get_ip(input_parameters, "tmax")),
          // Get pointers to output parameters
          delta_TT_op(get_op(output_parameters, "delta_TT"))
    {
    }
    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();

   private:
    // Pointers to input parameters
    const double* temp_ip;
    const double* tbase_ip;
    const double* topt_lower_ip;
    const double* topt_upper_ip;
    const double* tmax_ip;
    // Pointers to output parameters
    double* delta_TT_op;
    // Main operation
    void do_operation() const;
};

std::vector<std::string> delta_TT::get_inputs()
{
    return {
        "temp",   // deg. C
        "tbase",  // deg. C
        "topt_lower",   // deg. C
        "topt_upper",   // deg. C
        "tmax"    // deg. C
    };
}

std::vector<std::string> delta_TT::get_outputs()
{
    return {
        "delta_TT"  // deg. C * day / hr
    };
}

void delta_TT::do_operation() const
{
    // Collect inputs
    const double temp = *temp_ip;
    const double tbase = *tbase_ip;
    const double topt_lower = *topt_lower_ip;
    const double topt_upper = *topt_upper_ip;
    const double tmax = *tmax_ip;

    // Find the rate of change of the growing degree days
    double gdd_rate;
    if (temp <= tbase) {
        gdd_rate = 0.0;
    } else if (temp <= topt_lower) {
        gdd_rate = temp - tbase;
    } else if (temp > topt_lower && temp < topt_upper)  {
        gdd_rate = topt_lower-tbase;
    } else if(temp >= topt_upper && temp < tmax)  {
        gdd_rate = (tmax - temp) * (topt_lower - tbase) / (tmax - topt_upper);
    } else {
        gdd_rate = 0.0;
    }

    // Normalize to a rate per hour
    gdd_rate /= 24.0;

    // Update the output parameter list
    update(delta_TT_op, gdd_rate);
}

#endif