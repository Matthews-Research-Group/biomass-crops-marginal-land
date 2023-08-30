il_objfunlogistic <- function (parms,parameter_names,growing_season_list,observed,years){

  run_cwrfsoilwater=TRUE 

  # Modules required for miscanthus
  load("../../data/parameters/miscanthus_giganteus_ss_logistic_modules.rdata")
  load("../../data/parameters/miscanthus_giganteus_deriv_logistic_modules.rdata")
  load("../../data/parameters/miscanthus_giganteus_initial_state.rdata")
  load("../../data/parameters/miscanthus_giganteus_logistic_parameters.rdata")

  if(run_cwrfsoilwater){
    miscanthus_giganteus_deriv_logistic_modules = miscanthus_giganteus_deriv_logistic_modules[-3] #remove two_layer_soil_profile
    miscanthus_giganteus_initial_state =
      miscanthus_giganteus_initial_state[names(miscanthus_giganteus_initial_state)!=c('soil_water_content')]
  }
  miscanthus_giganteus_initial_state$Rhizome = 2.0  #t/ha,https://encyclopedia.pub/entry/31968  
  #obs data years 2006-2008 from doi: 10.1111/j.1757-1707.2011.01153.x
  #plant year: 2002. We start simulation from the plant year for a continuous run
  #This will be consistent with our multi-site validation
  initial_state = miscanthus_giganteus_initial_state

  miscanthus_giganteus_logistic_parameters$TTemr  = 400
  miscanthus_giganteus_logistic_parameters$alpha1 = 0.045 #use a larger value based on Charles' paper

  solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

  #these predicted values will be in vectors
  #although they are initialized as a single value
  predicted = list()
  predicted$stem = 0
  predicted$root = 0
  predicted$rhizome = 0
  predicted$leaf = 0
  obs_doys = observed$doy #get the obs' DOYs
  
  #hard-coded obs' DOY index in BioCro's hourly vector
  time_ind = c(1,1273,2713,4201,5665)
 
  for (i in 1:length(years))
  {
    year_i = years[i]

    growing_season = growing_season_list[[i]]

    if(run_cwrfsoilwater){
    #read in CWRF soil water for the initial soil water content and soil type
    #hard-coded path here!
      soil_data = 
      read.csv(paste0("~/MFEW_CWRF_soil_water/miscanthus/site_",2,"/cwrf_soilwater_",year_i,".csv"))
      miscanthus_giganteus_logistic_parameters$soil_type_indicator = soil_data$soiltype[1]
      growing_season$soil_water_content = soil_data$swc[soil_data$doy>=growing_season$doy[1] 
                                                   & soil_data$doy<=tail(growing_season$doy,1)]
    }
    partial_gro_function <- BioCroMis::partial_gro_solver(initial_state =  initial_state,
                            parameters = miscanthus_giganteus_logistic_parameters,
                            varying_parameters = growing_season,
                            steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
                            derivative_module_names = miscanthus_giganteus_deriv_logistic_modules,
                            arg_names = parameter_names,
                            solver=solver,verbose = FALSE)

    current_res <- partial_gro_function(parms)
    last_rhizome = tail(current_res$Rhizome,1)
    if(is.na(last_rhizome)) stop('NA rhizome')
    initial_state$Rhizome  = last_rhizome * 0.67  #assume 33% loss 

    if(i>4){ #only sum from 2006-2008
      predicted$stem    = predicted$stem    + current_res$Stem[time_ind]
      predicted$root    = predicted$root    + current_res$Root[time_ind]
      predicted$rhizome = predicted$rhizome + current_res$Rhizome[time_ind]
      predicted$leaf    = predicted$leaf    + current_res$Leaf[time_ind]
    }
  }
  #3-year averages
  predicted$stem    = predicted$stem/3
  predicted$root    = predicted$root/3
  predicted$rhizome = predicted$rhizome/3
  predicted$leaf    = predicted$leaf/3
  
# First data point of root should be zero. Observation does not make difference between alive and dead roots
  observed$root[1] = 0

#be careful of the variable names of the observation
#if they don't match, the program won't stop!
#you may see warning message like "stack imbalance". Very annoying!
  stemE    = observed$stem - predicted$stem 
  rootE    = observed$root - predicted$root 
  rhizomeE = observed$rhizome - predicted$rhizome 
  leafE    = observed$leaf - predicted$leaf 

#this weights apply to the DOYs' of errors where the last DOY is assigned a low weight
  E = sum(stemE^2* observed$wt)  + sum(rootE^2* observed$wt) + sum(rhizomeE^2* observed$wt)  + sum(leafE^2* observed$wt)

  if(is.na(E)) {
    print(parms)
    E=999
  }
 #if kLeaf_emr + kStem_emr >1, discard it 
  if(parms[2] + parms[3] >1 ) E=999
 return(E)
}
