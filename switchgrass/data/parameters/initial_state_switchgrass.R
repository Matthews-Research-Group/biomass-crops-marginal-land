initial_state_switchgrass = with(list(), {
  datalines =
    "symbol value
    Rhizome 0.02                                #Mg/ha
    Leaf    0.01                                #Mg/ha, 80% of total seed mass per land area
    Stem    0.01                                #Mg/ha, 10% of total seed mass per land area
    Root    0.02                                #Mg/ha, 10% of total seed mass per land area
    Grain   0.0                                #Mg/ha # treating this as the Pod/reproductive (minus seed)
    soil_water_content 0.37                    # Needed in two_layer_soil_profile module
    LeafN 2
    TTc 0                                     # Accumulated thermal time at the start  
    LeafLitter 0                              # Leaf litter at the start
    RootLitter 0                              # Root Litter at the start
    RhizomeLitter 0                           # Rhizome litter at the start
    StemLitter 0                              # Stem Litter at the start
     rhizome_senescence_index 0               # Not used but needed  in the module 'partitioning_growth'
    ## needed for two layer soil profile
    cws1    0.185                             # equal to initial soil_water_content, from Justins sorghum parameters (see 10/3/19 email -mlm)
    cws2    0.185                              # equal to initial soil_water_content, from Justins sorghum parameters (see 10/3/19 email -mlm)
    DVI -1                                    # Used in development_index.h,miscanthus_partitioning_coefficient_logistic,thermaltime_development_rate_calculator, senescence_coefficient_logistic
    "
  data_frame = read.table(textConnection(datalines), header=TRUE)
  values = as.list(data_frame$value)
  names(values) = data_frame$symbol
  values
})
