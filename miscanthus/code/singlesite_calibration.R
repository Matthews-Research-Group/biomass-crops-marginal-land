library(BioCroMis)
#Yufeng He: Dec 12, 2022

run_cwrfsoilwater=FALSE
# parameters and weather inputs

load("../data//parameters/miscanthus_giganteus_initial_state.rdata")
load("../data//parameters/miscanthus_giganteus_logistic_parameters.rdata")
load("../data//parameters/miscanthus_giganteus_ss_logistic_modules.rdata")
load("../data//parameters/miscanthus_giganteus_deriv_logistic_modules.rdata")

miscanthus_giganteus_logistic_parameters$TTemr  = 400

miscanthus_giganteus_logistic_parameters$alpha1 = 0.045 #0.04-0.05

parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr",
                            "alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")

##this chunk can make leaf&stomata water stress to a constant value
# miscanthus_giganteus_logistic_parameters$LeafWS <- 1
# miscanthus_giganteus_logistic_parameters$StomataWS <- 1
# ss_module_list = miscanthus_giganteus_ss_logistic_modules
# ss_modules     = miscanthus_giganteus_ss_logistic_modules
# ss_modules <- ss_modules[c(-which(ss_module_list=="stomata_water_stress_linear"),
#                            -which(ss_module_list=="leaf_water_stress_exponential"))]
# ss_modules <- ss_modules[c(-which(ss_module_list=="leaf_water_stress_exponential"))]
# miscanthus_giganteus_ss_logistic_modules = ss_modules
## 

# #these values are from the optimization results on BioCluster
# x=c(-0.000555,   0.211191,    0.671859,
#     1.011031,   -0.228122,    0.081572,  -38.765151,    0.160947,  -39.855814)

#this one is for turning OFF cwrf soil water
x=c(-0.000362,    0.269631,    0.609721,
    29.939436,  -39.155449,    0.199284,  -36.468715,    2.891147,   -5.831675)

miscanthus_giganteus_logistic_parameters[parameters_to_optimize] = x

if(run_cwrfsoilwater){
  miscanthus_giganteus_deriv_logistic_modules = miscanthus_giganteus_deriv_logistic_modules[-3] #remove two_layer_soil_profile
  miscanthus_giganteus_initial_state =
    miscanthus_giganteus_initial_state[names(miscanthus_giganteus_initial_state)!=c('soil_water_content')]
}

miscanthus_giganteus_initial_state$Rhizome = 24.1

solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
            adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

weather_all<-read.csv(
  "../data/weather_data_NASA_POWER/BioCro_input_NASA/site_2_2002_2018_LowerTransmittance.csv")
result = list()
years = 2006:2008
for (i in 1:length(years)){
  year_i = years[i]
  growing_season_weather = weather_all[weather_all$year == year_i,]
  growing_season_weather =  growing_season_weather[growing_season_weather$doy>=106 &
                                                    growing_season_weather$doy<=350,]
  
  if(run_cwrfsoilwater){
    #read in CWRF soil water for the initial soil water content and soil type
    soil_data = 
      read.csv(paste0("../data/weather_data_NASA_POWER//cwrf_soil_data//site_",2,"/cwrf_soilwater_",year_i,".csv"))
    miscanthus_giganteus_logistic_parameters$soil_type_indicator = soil_data$soiltype[1]
    growing_season_weather$soil_water_content = soil_data$swc[soil_data$doy>=growing_season_weather$doy[1] 
                                                              & soil_data$doy<=tail(growing_season_weather$doy,1)]
  }
  
  
  result[[i]] <- Gro_solver(initial_state =  miscanthus_giganteus_initial_state,
                       parameters = miscanthus_giganteus_logistic_parameters,
                       varying_parameters = growing_season_weather,
                       steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
                       derivative_module_names = miscanthus_giganteus_deriv_logistic_modules, solver=solver,verbose = FALSE)
  
  
  ##########################################################################################
  # Correcting for winter loss of atem based on 0.07 tons/ha per day
  non_frost_weather <- get_growing_season_climate(growing_season_weather, threshold_temperature = 0)
  
  for ( j in dim(non_frost_weather)[1]: (dim(growing_season_weather)[1])){
    result[[i]]$Stem[j] = result[[i]]$Stem[j] - (0.07/24)*(j- dim(non_frost_weather)[1])
  }
  ########################################################################################
}

source("plot_single_site.R")
#avg of 2006-2008
predicted = 0
for (i in 1:length(years)){
  predicted = predicted+result[[i]][,c("doy","Stem","Leaf","Rhizome","Root")]
}
predicted = predicted/length(years)
predicted = reshape2::melt(predicted,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root","Rhizome"))

observed_biomass <- read.csv("../data/biomass_observation/il_observation_without_index.csv")

names(observed_biomass)  = c("doy", "Rhizome",     "Leaf"     ,"Stem", "Root")
observed_biomass$Root[1] = 0 # First data point of root should be zero. Observation does not make difference between alive and dead roots
#observed <- observed_biomass[1:4, c("doy","Stem","Leaf","Rhizome","Root","lai")]
observed = observed_biomass[, c("doy","Stem","Leaf","Rhizome","Root")]
observed = reshape2::melt(observed,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root","Rhizome"))
names(observed)  = c("x","varname","y")
names(predicted) = c("x","varname","y")
xlabtitle = expression(paste("Day of Year"))
ylabtitle = expression(paste("Dry Biomass (Mg/ha)"))
BiomassPartitioning <- 
  Compare_Observed_and_Predicted(observeddata = observed,predicteddata=predicted  ,xlabtitle,ylabtitle)
plot(BiomassPartitioning)
# ggsave(filename="./BiomassPartitioning.png",dpi = 500, width =6.5, height =4.2)

doy = result[[1]]$doy
# 
# plot(doy,result[[1]]$DVI)
# lines(doy,rep(0,length(doy)),col='red')
# 
# plot(doy,result[[1]]$kRoot)
