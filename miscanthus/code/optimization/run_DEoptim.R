library(BioCroMis)
library(DEoptim)

run_cwrfsoilwater=TRUE

fn.observed = "il_observation_without_index.csv"
#Yufeng: the last DOY the obs has reduced STEM.
#since the model does not have STEM senescenece, it has a low weight on the last DOY
observed <- read.csv(fn.observed)
print(observed)

#climate file for running the model
fn.climate <- "site_2_2002_2018_LowerTransmittance.csv"
# growing season weather from IL 
weather_all <- read.csv(fn.climate)
partial_gro_list = list()

# Modules required for miscanthus
load("../../data/parameters/miscanthus_giganteus_ss_logistic_modules.rdata")
load("../../data/parameters/miscanthus_giganteus_deriv_logistic_modules.rdata")
load("../../data/parameters/miscanthus_giganteus_initial_state.rdata")
load("../../data/parameters/miscanthus_giganteus_logistic_parameters.rdata")

parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")
lower_bound_parameters <- c(-0.01   ,0.05 ,0.05)
upper_bound_parameters <- c(-0.00001,0.8  ,0.8)
lower_bound_parameters <- c(lower_bound_parameters,c(0 ,-10,0 ,-40,0 ,-40))
upper_bound_parameters <- c(upper_bound_parameters,c(10,0  ,40,0  ,40,0))

miscanthus_giganteus_logistic_parameters$TTemr  = 400
miscanthus_giganteus_logistic_parameters$alpha1 = 0.045 #use a larger value based on Charles' paper

if(run_cwrfsoilwater){
  miscanthus_giganteus_deriv_logistic_modules = miscanthus_giganteus_deriv_logistic_modules[-3] #remove two_layer_soil_profile
  miscanthus_giganteus_initial_state =
    miscanthus_giganteus_initial_state[names(miscanthus_giganteus_initial_state)!=c('soil_water_content')]
}

solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

miscanthus_giganteus_initial_state$Rhizome = 24.1  #t/ha,https://encyclopedia.pub/entry/31968  
#obs data years 2006-2008 from doi: 10.1111/j.1757-1707.2011.01153.x
#plant year: 2002. We start simulation from the plant year for a continuous run
#This will be consistent with our multi-site validation
years = 2006:2008 
for (i in 1:length(years)){
  year_i   = years[i] 
  growing_season_weather = weather_all[weather_all$year == year_i, ]
  
  growing_season <- growing_season_weather[with(growing_season_weather,doy >=106 & doy <=350),] 

  if(run_cwrfsoilwater){
    #read in CWRF soil water for the initial soil water content and soil type
    #hard-coded path here!
    soil_data = 
      read.csv(paste0("~/MFEW_CWRF_soil_water/miscanthus/site_",2,"/cwrf_soilwater_",year_i,".csv"))
    miscanthus_giganteus_logistic_parameters$soil_type_indicator = soil_data$soiltype[1]
    growing_season$soil_water_content = soil_data$swc[soil_data$doy>=growing_season$doy[1] 
                                                              & soil_data$doy<=tail(growing_season$doy,1)]
  }
  
  partial_gro_function <- partial_gro_solver(initial_state =  miscanthus_giganteus_initial_state,
                                             parameters = miscanthus_giganteus_logistic_parameters,
                                             varying_parameters = growing_season,
                                             steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
                                             derivative_module_names = miscanthus_giganteus_deriv_logistic_modules,
                                             arg_names = parameters_to_optimize,
                                             solver=solver,verbose = FALSE)
  partial_gro_list[[i]] = partial_gro_function
}

source("il_objfunlogistic.R")
cost_func <- function(x){
 	il_objfunlogistic(x,partial_gro_list,observed)
}
# maximum number of iterations
max.iter <- 500

set.seed(1234)
# Call DEoptim function to run optimization
parVars <- c('il_objfunlogistic','partial_gro_list','observed')

cl <- makeCluster(12)
clusterExport(cl, parVars,envir=environment())
optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
                      control=list(VTR=10,itermax=max.iter,parallelType=1,packages=c('BioCroMis'),parVar=parVars,cl=cl))

#optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
#                      control=list(VTR=10,itermax=max.iter,parallelType=0,packages=c('BioCroMis')))

opt_result <- data.frame(para=optim_result$par,MSE=optim_result$value)

saveRDS(opt_result,'opt_result_DEoptim_3year_run_r1.rds')
