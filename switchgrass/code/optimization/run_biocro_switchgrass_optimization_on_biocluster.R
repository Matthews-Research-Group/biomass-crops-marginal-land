library(BioCroMis)
library(DEoptim)

run_cwrfsoilwater=TRUE

# reading observation and climate data
observation_switchgrass_IL = read.csv("../../data/biomass_observation//observations_switchgrass_2006_2008.csv")
#ND location: 46.77,-100.92
observation_switchgrass_ND = read.csv("../../data/biomass_observation//lit_data.csv")

#exlude the last point in Dec.
observation_switchgrass_IL <- observation_switchgrass_IL[1:4,] 
print(observation_switchgrass_IL)

weather_urbana = read.csv('../../data/weather/NASA_data/BioCroInputs/site_1_lowerTransmittance.csv')
weather_ND     = read.csv('../../data/weather/NASA_data/BioCroInputs/site_99_lowerTransmittance.csv')

# Modules/initial conditions required for running biocro
source("../../data/parameters/ss_modules_list_switchgrass.R")
source("../../data/parameters/deriv_modules_list_with_switchgrass_senescence.R")
source("../../data/parameters/initial_state_switchgrass.R")
source("../../data/parameters/parameters_list_switchgrass_senescence.R")

# upload R script files to source objective functions
source("objectivefunction_switchgrass.R")


# Running optimization for the logistical phase
parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot","betaRoot",
                             "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
#parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
parameters_list_switchgrass$TTc_leafsenescence_threshold = 5 
parameters_list_switchgrass$soil_type_indicator = 4 #YH: this by default was 10!?  
parameters_list_switchgrass$iSp    = 1.1
#parameters_list_switchgrass$TTemr  = 300 

#parameters_list_switchgrass$tbase      = 10
#parameters_list_switchgrass$topt_upper = 31
#parameters_list_switchgrass$tmax       = 40 
parameters_list_switchgrass$upperT  = 28 

years  = 2006:2008
years2 = 2001:2002
if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass = 
  initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
  cwrf_soilwater = list()
  for(i in 1:length(years)){
   year = years[i]
#first site's soil water
   soil_data = read.csv(paste0('~/MFEW_CWRF_soil_water/switchgrass/site_1/cwrf_soilwater_',year,'.csv'))
   cwrf_soilwater[[i]] = soil_data[c('doy','swc')] 
  }

  cwrf_soilwater2 = list()
  for(i in 1:length(years2)){
   year2 = years2[i]
#second site's soil water
   soil_data = read.csv(paste0('~/MFEW_CWRF_soil_water/switchgrass/site_99/cwrf_soilwater_',year2,'.csv'))
   cwrf_soilwater2[[i]] = soil_data[c('doy','swc')] 
  }
}else{
  ##turn off water stress
  ss_modules_list_swithcgrass=ss_modules_list_swithcgrass[-c(3)]
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  parameters_list_switchgrass$LeafWS=1
}

lower_bound_parameters <- c(0 ,-20,0 ,-20,0 ,-20)
upper_bound_parameters <- c(20,0  ,20,0  ,20,0  )
#lower_bound_parameters <- c(-20 ,-20 ,-20 ,-20 ,-20 ,-20)
#upper_bound_parameters <- c(20  ,20  ,20  ,20  ,20  ,20  )
lower_bound_parameters <- c(lower_bound_parameters,c(-0.01   ,0.01,0.01,0.0  ,100))
upper_bound_parameters <- c(upper_bound_parameters,c(-0.00001,0.80,0.80,0.002,500))

wt <- list(leaf=1,stem=1,aboveground = 1, root = 1, rhizome = 1)

solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

set.seed(1234)

optsolver <- list()
initial_state_switchgrass$Rhizome = 10.0
initial_state_switchgrass$Root    = 0.1
for (i in 1:length(years)){
   ind = which(weather_urbana$year==years[i])
   growing_season = weather_urbana[ind,]
   #growing_season = get_growing_season_climate(growing_season,threshold_temperature = 0)
   growing_season = growing_season[growing_season$doy>=105 & growing_season$doy<=290,]
   if(run_cwrfsoilwater){
     tmp = cwrf_soilwater[[i]]
     growing_season$soil_water_content = tmp$swc[tmp$doy>=growing_season$doy[1] & tmp$doy<=tail(growing_season$doy,1)] 
   }else{
   #this currently does not work 
      initial_state_switchgrass$soil_water_content = soilmoisture[[i]]
   }
   optsolver[[i]] <- partial_gro_solver(initial_state =  initial_state_switchgrass,
                        parameters = parameters_list_switchgrass,
                        varying_parameters = growing_season,
                        steady_state_module_names = ss_modules_list_swithcgrass,
                        derivative_module_names = deriv_modules_list_switchgrass, parameters_to_optimize,solver=solver,verbose = FALSE)
}

optsolver2 <- list()
initial_state_switchgrass$Rhizome = 4.0
initial_state_switchgrass$Root    = 0.1
for (i in 1:length(years2)){
   ind = which(weather_ND$year==years2[i])
   growing_season = weather_ND[ind,]
   #growing_season = get_growing_season_climate(growing_season,threshold_temperature = 0)
   growing_season = growing_season[growing_season$doy>=130 & growing_season$doy<=260,]
   if(run_cwrfsoilwater){
     tmp = cwrf_soilwater2[[i]]
     growing_season$soil_water_content = tmp$swc[tmp$doy>=growing_season$doy[1] & tmp$doy<=tail(growing_season$doy,1)] 
   }else{
   
      initial_state_switchgrass$soil_water_content = soilmoisture[[i]]
   }
   optsolver2[[i]] <- partial_gro_solver(initial_state =  initial_state_switchgrass,
                        parameters = parameters_list_switchgrass,
                        varying_parameters = growing_season,
                        steady_state_module_names = ss_modules_list_swithcgrass,
                        derivative_module_names = deriv_modules_list_switchgrass, parameters_to_optimize,solver=solver,verbose = FALSE)
}

cost_func <- function(x){
  objectivefunction_switchgrass(x,optsolver,optsolver2,observation_switchgrass_IL,observation_switchgrass_ND,wt)
}

set.seed(1234)
# maximum number of iterations
max.iter <- 1000

# Call DEoptim function to run optimization
parVars <- c('objectivefunction_switchgrass','optsolver','optsolver2',
              'observation_switchgrass_IL', 'observation_switchgrass_ND','wt')

cl <- makeCluster(8)
clusterExport(cl, parVars,envir=environment())

optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
                      control=list(VTR=20,itermax=max.iter,parallelType=1,packages=c('BioCroMis'),cluster=cl))
#optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
#                      control=list(VTR=20,itermax=max.iter,parallelType=0,packages=c('BioCroMis')))

opt_result = data.frame(optim_result$par,MSE=optim_result$value)

saveRDS(opt_result,'opt_result_DEoptim_cwrfsoil_NDsite_r3.rds')
