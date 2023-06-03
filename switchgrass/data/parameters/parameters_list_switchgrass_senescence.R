parameters_list_switchgrass = with(list(), {
  datalines =
    "symbol value
  lat                 26.639          # latitude, needed in many modules
  soil_clod_size      0.04        # needed in two_layer_soil_profile module
  soil_reflectance    0.2         # needed in two_layer_soil_profile module
  soil_transmission   0.01        # needed in two_layer_soil_profile module
  specific_heat       1010        # needed in two_layer_soil_profile module
  stefan_boltzman     5.67e-8     # needed in two_layer_soil_profile module
  iSp                 1.9         # needed in parameter_calculator (lai /leaf in hectare/ton)      
  Sp_thermal_time_decay  0.00038    # neeeded in parameter_calculator, linear decline in iSp with thermal time hecatre/ton/thermaltime
  mrc1             0.0008333      # needed in no_leaf_resp_partitioning_growth_calculator, 0.02/24, Bunce 1995, https://doi.org/10.1006/anbo.1996.0061
  mrc2             0.000025       # needed in no_leaf_resp_partitioning_growth_calculator in 03*mrc1 #0.000015    #0.03
  nlayers            10           # needed in c4_canopy module , number of canopy layers
  chil                1           # needed in c4_canopy module , leaf angle distribution
  kd                 0.3          # needed in c4_canopy module , leaf angle distribution. An Introduction to Environmental Biophysics, 2nd Edition, Figure 15.4, pg 254   
  heightf            3            # needed in c4_canopy module
  leafwidth           0.04         # needed in c4_canopy module
  et_equation         0            # needed in c4_canopy module
  growth_respiration_fraction 0    # needed in multilayer_canopy_integrator
  vmax1 32                         # needed in c4_canopy module
  alpha1 0.026                      # needed in c4_canopy module
  kparm 0.7                        # needed in c4_canopy module
  theta 0.83                       # needed in c4_canopy module
  beta 0.93                        # needed in c4_canopy module
  Rd 0.8                           # needed in c4_canopy module
  Catm 400                         # needed in c4_canopy module
  b0 0.082                          # needed in c4_canopy module
  b1 2.15                             # needed in c4_canopy module
  water_stress_approach     1      # needed in c4_canopy module
  upperT 41                      # needed in c4_canopy module
  lowerT 6                        # needed in c4_canopy module
  lnfun    0                       # needed in c4_canopy 1 means use N relationship, 0 means use provided values of photosynthesis parameters
  nvmaxb1 0.6938                   # needed in c4_canopy module, not used because lnfun =0
  nvmaxb0 -16.25                   # needed in c4_canopy module, not used because lnfun =0
  nalphab1 0.000488                # needed in c4_canopy module, not used because lnfun =0
  nalphab0 0.02367                 # needed in c4_canopy module, not used because lnfun =0
  nRdb1 0.1247                     # needed in c4_canopy module, not used because lnfun =0
  nRdb0 -4.5917                    # needed in c4_canopy module, not used because lnfun =0
  nileafn 85                       # needed in c4_canopy module, not used because lnfun = 0
  nkln 0.5                         # needed in c4_canopy module
  nkpLN 0.17                       # needed in c4_canopy module
  nlnb0 -5                         # needed in c4_canopy module
  nlnb1 18                         # needed in c4_canopy module
  kpLN     0                       # needed in c4 canopy module for changing leaf N within canopy [currently set to 0]
  TTemr  500                      # neededin in thermaltime_development_rate_calculator module
  TTveg  1500                     # based on miscanthus
  TTrep       2500                # based on miscanthus
  alphaLeaf   7.752615    # needed in miscanthus_partitioning_coefficient_logistic
     betaLeaf   -8.011449     # needed in miscanthus_partitioning_coefficient_logistic
     alphaRoot   5.750015         # needed in miscanthus_partitioning_coefficient_logistic
     betaRoot    -5.428443      # needed in miscanthus_partitioning_coefficient_logistic
     alphaStem   5.120094      # needed in miscanthus_partitioning_coefficient_logistic
     betaStem    -2.070156  # needed in miscanthus_partitioning_coefficient_logistic
     kRhizome_emr -0.00004        # needed in miscanthus_partitioning_coefficient_logistic
     kLeaf_emr    0.2           # needed in miscanthus_partitioning_coefficient_logistic
     kStem_emr     0.4          # needed in miscanthus_partitioning_coefficient_logistic
     kRoot_emr    0.4       # needed in miscanthus_partitioning_coefficient_logistic
     kGrain_emr    0.0               # needed in miscanthus_partitioning_coefficient_logistic
     remobilization_fraction    0    # needed in thermal_time_senescence_logistic
     tbase      4                      # base temperatue for thermal time calculations
     topt_lower  28                    # lower range of optiaml temp for thermal time calculations
     topt_upper  38                   # upper range of optiaml temp for thermal time calculations
  tmax        45                  # maximum temperature for for thermal time calculations
  phi1       0.01                 # needed in two_layer_soil_profile module
  phi2       1.5                  # needed in two_layer_soil_profile module        #jaiswal 2017 #10       
  soil_depth  0.38                   # needed in two_layer_soil_profile module
  soil_type_indicator   10         # needed in soil_type_selector
  rsec       0.2                  # needed in two_layer_soil_profile module
  LeafN_0   2                     # needed in linear_vmax_from_leaf_n
  vmax_n_intercept  0             # needed linear_vmax_from_leaf_n and parameter_calculator
  timestep   1
  alphab1    0
  alpha1     0
  lnb0   0
  lnb1   0
  soil_depth1    0.0 # from Justins sorghum parameters (see 10/3/19 email -mlm)
  soil_depth2    0.185 # from Justins sorghum parameters (see 10/3/19 email -mlm)
  soil_depth3    0.37 # from Justins sorghum parameters (see 10/3/19 email -mlm)
  wsFun  2 #exponential stomata water stress
  hydrDist   0 #copied from sorghum parameter file
  rfl    0.2 #copied from sorghum parameter file
  rsdf   0.44 #copied from sorghum parameter filep 
  atmospheric_pressure 101325
  atmospheric_transmittance 0.85
  atmospheric_scattering  0.3
  par_energy_fraction_of_sunlight 0.5
  par_energy_content 0.235
  lai_max 5
  remobilization_fraction_leaf_to_rhizome 0.5
  phi_waterstress_induced_leafsenescence  124.31
  sene_factor_when_sws_eq_0 1.0
  sene_factor_when_sws_eq_1 0.0
  Tfrostlow -2
  Tfrosthigh 0
  remobilization_fraction_leaf_to_rhizome 0.5
  remobilization_fraction_stem_to_rhizome 0.5
  remobilization_fraction_root_to_rhizome 0.5
  remobilization_fraction_rhizome_to_rhizome 0.5
  leaf_turnover_rate 0.000210083 
  stem_turnover_rate 0.0
  root_turnover_rate 0.0
  rhizome_turnover_rate 0.0
  TTc_leafsenescence_threshold 0
  TTc_stemsenescence_threshold 10000
  TTc_rootsenescence_threshold 10000
  TTc_rhizomesenescence_threshold 10000
  "
    data_frame = read.table(textConnection(datalines), header=TRUE)
    values = as.list(data_frame$value)
    names(values) = data_frame$symbol
    values
})
