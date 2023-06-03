library(BioCroMis)
library(DEoptim)

run_cwrfsoilwater=TRUE

# reading observation and climate data
observation_switchgrass_IL = read.csv("../../data/biomass_observation//observations_switchgrass_2006_2008.csv")

# File containing information to map result index and observation corresponding to date of measurements
#observation_result_mapping <- read.csv("../../data/biomass_observation//switchgrass_observation_result_mapping_for_optimization.csv")

observation_switchgrass_IL <- observation_switchgrass_IL[1:4,] 
#observation_result_mapping <- observation_result_mapping[1:4,]

weather_urbana = read.csv('../../data/weather/NASA_data/BioCroInputs/site_1_lowerTransmittance.csv')

# Modules/initial conditions required for running biocro
source("../../data/parameters/ss_modules_list_switchgrass.R")
source("../../data/parameters/deriv_modules_list_with_switchgrass_senescence.R")
source("../../data/parameters/initial_state_switchgrass.R")
source("../../data/parameters/parameters_list_switchgrass_senescence.R")

# upload R script files to source objective functions
source("objectivefunction_switchgrass.R")


# Running optimization for the logistical phase
#parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot","betaRoot",
#                             "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate")
parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
parameters_list_switchgrass$TTc_leafsenescence_threshold = 5 
parameters_list_switchgrass$soil_type_indicator = 4 #YH: this by default was 10!?  
parameters_list_switchgrass$iSp = 1.7
years = 2006:2008
if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass = 
  initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
  cwrf_soilwater = list()
  k=1
  for(year in years){
   soil_data = read.csv(paste0('~/MFEW/switchgrass/site_1/cwrf_soilwater_',year,'.csv'))
   cwrf_soilwater[[k]] = soil_data[c('doy','swc')] 
   k=k+1
  }
}else{
  ##turn off water stress
  ss_modules_list_swithcgrass=ss_modules_list_swithcgrass[-c(3)]
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  parameters_list_switchgrass$LeafWS=1
}

#lower_bound_parameters <- c(0 ,-40,0 ,-40,0 ,-40,    -0.01   ,0.01,0.01,0.0)
#upper_bound_parameters <- c(40,0  ,40,0  ,40,0  ,    -0.00001,0.20,0.20,0.001)
lower_bound_parameters <- c(-0.01   ,0.01,0.01,0.0  ,  100)
upper_bound_parameters <- c(-0.00001,0.50,0.50,0.001,  400)

wt <- list(aboveground = 1, root = 2, rhizome = 2)

solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

set.seed(1234)

optsolver <- list()
initial_state_switchgrass$Rhizome = 14.0
initial_state_switchgrass$Root    = 0.1
for (i in 1:length(years)){
   ind = which(weather_urbana$year==years[i])
   growing_season = weather_urbana[ind,]
   growing_season = get_growing_season_climate(growing_season,threshold_temperature = 0)
   if(run_cwrfsoilwater){
     tmp = cwrf_soilwater[[i]]
     growing_season$soil_water_content = tmp$swc[tmp$doy>=growing_season$doy[1] & tmp$doy<=tail(growing_season$doy,1)] 
   }else{
      initial_state_switchgrass$soil_water_content = soilmoisture[[4+i]]
   }
   optsolver[[i]] <- partial_gro_solver(initial_state =  initial_state_switchgrass,
                        parameters = parameters_list_switchgrass,
                        varying_parameters = growing_season,
                        steady_state_module_names = ss_modules_list_swithcgrass,
                        derivative_module_names = deriv_modules_list_switchgrass, parameters_to_optimize,solver=solver,verbose = FALSE)
}

cost_func <- function(x){
  objectivefunction_switchgrass(x,optsolver,observation_switchgrass_IL,wt)
}

set.seed(1234)
# maximum number of iterations
max.iter <- 500

# Call DEoptim function to run optimization
parVars <- c('objectivefunction_switchgrass','optsolver','observation_switchgrass_IL',
                                           'wt')
cl <- makeCluster(8)
clusterExport(cl, parVars,envir=environment())

optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
                      control=list(VTR=30,itermax=max.iter,parallelType=1,packages=c('BioCroMis'),cluster=cl))

opt_result = data.frame(optim_result$par,MSE=optim_result$value)

saveRDS(opt_result,'opt_result_DEoptim_cwrfsoil_r1.rds')

stopCluster(cl)
