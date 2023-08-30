library(BioCroMis)
library(DEoptim)

fn.observed = "il_observation_without_index.csv"
#Yufeng: the last DOY the obs has reduced STEM.
#since the model does not have STEM senescenece, it has a low weight on the last DOY
observed <- read.csv(fn.observed)
print(observed)

parameter_names <- c("kRhizome_emr","kLeaf_emr","kStem_emr","alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")
lower_bound_parameters <- c(-0.01   ,0.05 ,0.05)
upper_bound_parameters <- c(-0.00001,0.8  ,0.8)
lower_bound_parameters <- c(lower_bound_parameters,c(0 ,-10,0 ,-40,0 ,-40))
upper_bound_parameters <- c(upper_bound_parameters,c(10,0  ,40,0  ,40,0))

#climate file for running the model
fn.climate <- "site_2_2002_2018_LowerTransmittance.csv"
# growing season weather from IL 
weather_all <- read.csv(fn.climate)
growing_season_list   = list()

years = 2002:2008 
for (i in 1:length(years)){
  year_i   = years[i] 
  growing_season_weather = weather_all[weather_all$year == year_i, ]
  
  growing_season <- growing_season_weather[with(growing_season_weather,doy >=106 & doy <=350),] 

  growing_season_list[[i]] = growing_season
  
}

source("il_objfunlogistic_continue.R")
cost_func <- function(x){
 	il_objfunlogistic(x,parameter_names,growing_season_list,observed,years)
}
# maximum number of iterations
max.iter <- 500

set.seed(1234)
# Call DEoptim function to run optimization
parVars <- c('il_objfunlogistic','parameter_names','growing_season_list','observed','years')

cl <- makeCluster(12)
clusterExport(cl, parVars,envir=environment())
optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
                      control=list(VTR=10,itermax=max.iter,parallelType=1,packages=c('BioCroMis'),parVar=parVars,cl=cl))

#optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
#                      control=list(VTR=10,itermax=max.iter,parallelType=0,packages=c('BioCroMis')))

opt_result <- data.frame(para=optim_result$par,MSE=optim_result$value)

saveRDS(opt_result,'opt_result_DEoptim_continue_run_r1.rds')
