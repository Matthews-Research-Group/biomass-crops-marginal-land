#ifndef INCIDENT_SHORTWAVE_PAR_FROM_CWRF_H
#define INCIDENT_SHORTWAVE_PAR_FROM_CWRF_H

#include "../modules.h"
#include "../constants.h"

/**
 * @class incident_shortwave_par_from_cwrf
 * 
 * @brief Uses direct and diffure PAR radiation from CWRF model.
 * 
 * 
 * Here the input quantity `solar` represents the total flux density of direct (beam) PAR photons  which comes out of CWRF model.
 * Diffused radiation (PAR) `solardiff` also comed out of CWRF in
 * units of micromol / m^2 / s. However, the output parameters will be converted into the
 * SI standard energy flux density units of J / m^2 / s or equivalently W / m^2.
 * 
 * The `par_energy_fraction_of_sunlight` quantity should represent the fraction of solar energy
 * that lies in the PAR band. This value is often near 0.5. This is used in calculating nir.
 */
class incident_shortwave_par_from_cwrf : public SteadyModule
{
   public:
     incident_shortwave_par_from_cwrf(
        const std::unordered_map<std::string, double>* input_parameters,
        std::unordered_map<std::string, double>* output_parameters)
        :  // Define basic module properties by passing its name to its parent class
          SteadyModule("incident_shortwave_par_from_cwrf"),
          // Get pointers to input parameters
          solar(get_input(input_parameters, "solar")),
          solardiff(get_input(input_parameters, "solardiff")),
          irradiance_direct_fraction(get_input(input_parameters, "irradiance_direct_fraction")),
          irradiance_diffuse_fraction(get_input(input_parameters, "irradiance_diffuse_fraction")),
          par_energy_fraction_of_sunlight(get_input(input_parameters, "par_energy_fraction_of_sunlight")),
          par_energy_content(get_input(input_parameters, "par_energy_content")),
          // Get pointers to output parameters
          par_incident_direct_op(get_op(output_parameters, "par_incident_direct")),
          par_incident_diffuse_op(get_op(output_parameters, "par_incident_diffuse")),
          nir_incident_direct_op(get_op(output_parameters, "nir_incident_direct")),
          nir_incident_diffuse_op(get_op(output_parameters, "nir_incident_diffuse"))
    {
    }
    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_description();

   private:
    // References to input parameters
    double const& solar;
    double const& solardiff;
    double const& irradiance_direct_fraction;
    double const& irradiance_diffuse_fraction;
    double const& par_energy_fraction_of_sunlight;
    double const& par_energy_content;
    // Pointers to output parameters
    double* par_incident_direct_op;
    double* par_incident_diffuse_op;
    double* nir_incident_direct_op;
    double* nir_incident_diffuse_op;
    // Main operation
    void do_operation() const;
};

std::vector<std::string>  incident_shortwave_par_from_cwrf::get_inputs()
{
    return {
        "solar",                            // micromol / (m^2 beam) / s [area perpendicular to beam]
        "solardiff",                        // micromol / (m^2 diffused) / s 
        "irradiance_direct_fraction",       // dimensionless
        "irradiance_diffuse_fraction",      // dimensionless
        "par_energy_fraction_of_sunlight",  // dimensionless
        "par_energy_content"                // J / micromol
    };
}

std::vector<std::string>  incident_shortwave_par_from_cwrf::get_outputs()
{
    return {
        "par_incident_direct",   // J / (m^2 beam) / s [area perpendicular to beam]
        "par_incident_diffuse",  // J / m^2 / s        [through any plane]
        "nir_incident_direct",   // J / (m^2 beam) / s [area perpendicular to beam]
        "nir_incident_diffuse"   // J / m^2 / s        [through any plane]
    };
}

void  incident_shortwave_par_from_cwrf::do_operation() const
{
    // Calculate the PAR flux density in the direct beam and as diffuse radiation
    const double par_incident_direct = solar  * par_energy_content;
    const double par_incident_diffuse = solardiff * par_energy_content;

    // Calculate the NIR flux density in the direct beam and as diffuse radiation
    // Here we just assume that a constant fraction of solar energy lies in the PAR band,
    //  with the remainder being in the NIR band.
    const double nir_incident_direct = par_incident_direct * (1.0 / par_energy_fraction_of_sunlight - 1.0);
    const double nir_incident_diffuse = par_incident_diffuse * (1.0 / par_energy_fraction_of_sunlight - 1.0);

    // Check for error conditions
    std::map<std::string, bool> errors_to_check = {
        {"par_energy_fraction_of_sunlight cannot be zero", par_energy_fraction_of_sunlight == 0}  // divide by zero
    };

    check_error_conditions(errors_to_check, get_name());

    // Update outputs
    update(par_incident_direct_op, par_incident_direct);
    update(par_incident_diffuse_op, par_incident_diffuse);
    update(nir_incident_direct_op, nir_incident_direct);
    update(nir_incident_diffuse_op, nir_incident_diffuse);
}

#endif
