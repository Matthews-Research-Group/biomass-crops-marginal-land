library(BioCroMis) #change this to your BioCro library name
library(ggplot2)
library(epiR)
run_cwrfsoilwater=TRUE
#obs data was pre-processed. See "averaging_biomass.R" on how it was done
obs = read.csv('../data/biomass_observation/switchgrass_observation_averaged_plantyear_maturestand_r2.csv')
latlon_id <- read.csv('../data/weather/NASA_data/unique_latlon_and_id.csv')

rmse<-function(obs,pred){
  NAs  = which(is.na(obs) | is.na(pred))
  pred = pred[-NAs]
  obs  = obs[-NAs]
  rmse = sqrt(sum((pred-obs)^2)/length(obs))
}
ttc_function<-function(temp_array){
  tbase = 4
  topt_lower = 28
  topt_upper = 38
  tmax = 45
  TTc = 0
  for (i in 1:length(temp_array)){
    temp = temp_array[i]
    if (temp <= tbase){
      gdd_rate = 0.0
    }else if (temp <= topt_lower){
      gdd_rate = temp - tbase
    }else if (temp > topt_lower & temp < topt_upper){
      gdd_rate = topt_lower-tbase
    }else if (temp >= topt_upper &  temp < tmax){
      gdd_rate = (tmax - temp) * (topt_lower - tbase) / (tmax - topt_upper)
    }else{
      gdd_rate = 0.0
    }
    #  Normalize to a rate per hour
    gdd_rate = gdd_rate/ 24.0
    
    TTc = TTc + gdd_rate
  }
  return(TTc)
}
years = 2001:2015
TTc_urbana=data.frame(year = years,total_TTC=NA)
weather_urbana = read.csv('../data/weather/NASA_data/BioCroInputs/site_1_lowerTransmittance.csv')
for (i in 1:length(years)){
  weather_data = weather_urbana[which(weather_urbana$year==years[i]),]
  growing_season <- get_growing_season_climate(weather_data, threshold_temperature = 0)
  tmp = growing_season$temp
  TTc_urbana$total_TTC[i] = ttc_function(tmp)
}

source("../data/parameters/initial_state_switchgrass.R")
source("../data/parameters/ss_modules_list_switchgrass.R")
source("../data/parameters/deriv_modules_list_with_switchgrass_senescence.R")
source("../data/parameters/parameters_list_switchgrass_senescence.R")

parameters_list_switchgrass$TTc_leafsenescence_threshold = 5
parameters_list_switchgrass$iSp   = 1.1

parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf", "alphaRoot", "betaRoot",
                            "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
x=c( 12.172727,   -8.076308,   14.594526,  -12.076529,    9.122148,  -13.958066,
     -0.003602,    0.019288,    0.015474,    0.000361,  100.848453)

parameters_list_switchgrass[parameters_to_optimize] = x

parameters_list_switchgrass0 = parameters_list_switchgrass

# initial_state_switchgrass$Rhizome = 14
initial_state_switchgrass$Rhizome = 0.1
initial_state_switchgrass$Root    = 0.1
parameters_list_switchgrass$soil_type_indicator = 4
solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
           adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

#http://dx.doi.org/10.1080/07352680500316433
#https://doi.org/10.2135/cropsci2004.1391

if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass =
  initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
}

initial_state_switchgrass0 = initial_state_switchgrass

rhizome_loss =  0.33
# root_loss    =  0.33
post_season_loss_per_day = 0.0213 #t/ha/day
#ref paper: [Over-winter yield decline in Switchgrass and Miscanthus]

#This paper suggests loss of 35% to 45% due to machine harvest
#https://www.sciencedirect.com/science/article/pii/S0961953409000166
# 
predicted = c()
avg_temp  = c()
post_season_loss = c(1:dim(obs)[1])*NA
ttc_scaling_factor_all=c()
for (i in 1:dim(obs)[1]){
  print(i)
  #re-init for every location
  initial_state_switchgrass = initial_state_switchgrass0
  #
  parameters_list_switchgrass$lat=obs$lat[i]
  plant_year = obs$plantyear[i]
  harvest_year = obs$harvest_year[i]
  harvest_month = obs$harvest_month[i]
  harvest_day = obs$harvest_day[i]
  harvest_date = as.Date(paste(harvest_year,harvest_month,harvest_day,sep="-"))
  # harvest_doy = as.numeric(format(harvest_date,"%j"))
  if(harvest_year-plant_year+1<3){
    predicted[i] = NA
    next
  }
  #here I use lat/lon to match site id so that the id of the weather data could be different from 
  #the current obs data's id order
  #it's much safer this way instead of matching id index directly
  site.id <- latlon_id$id[latlon_id$lat==obs$lat[i] & latlon_id$lon==obs$lon[i]]
  biocroinput <- read.csv(paste0('../data/weather/NASA_data/BioCroInputs/site_',site.id,'_lowerTransmittance.csv')) 
  
  if(harvest_month<6){
    harvest_year_actual = harvest_year -1
  }else{
    harvest_year_actual = harvest_year
  }
  
  for (year_j in plant_year:harvest_year_actual){
    weather_data = biocroinput[which(biocroinput$year==year_j),]
    growing_season <- get_growing_season_climate(weather_data, threshold_temperature = 0)

    #read in CWRF soil water for the initial soil water content and soil type
    soil_data = read.csv(paste0("../data/weather/CWRF_soil_water/site_",site.id,"/cwrf_soilwater_",year_j,".csv"))
    # parameters_list_switchgrass$soil_type_indicator = soil_data$soiltype[1]
    if(run_cwrfsoilwater){
      growing_season$soil_water_content = soil_data$swc[soil_data$doy>=growing_season$doy[1] 
                                                        & soil_data$doy<=tail(growing_season$doy,1)]
    }else{
      initial_state_switchgrass$soil_water_content = soil_data$swc[soil_data$doy==growing_season$doy[1] & soil_data$hour==0]
    }
    
    total_ttc = ttc_function(growing_season$temp)
    ttc_scaling_factor = total_ttc / TTc_urbana$total_TTC[TTc_urbana$year==year_j]
    parameters_list_switchgrass$TTemr = parameters_list_switchgrass0$TTemr * ttc_scaling_factor
    parameters_list_switchgrass$TTveg = parameters_list_switchgrass0$TTveg * ttc_scaling_factor
    parameters_list_switchgrass$TTrep = parameters_list_switchgrass0$TTrep * ttc_scaling_factor
    ttc_scaling_factor_all = c(ttc_scaling_factor_all,ttc_scaling_factor)
    
    result <- Gro_solver(initial_state =  initial_state_switchgrass,
                         parameters = parameters_list_switchgrass,
                         varying_parameters = growing_season,
                         steady_state_module_names = ss_modules_list_swithcgrass,
                         derivative_module_names = deriv_modules_list_switchgrass, solver=solver,verbose = FALSE)
    n=dim(result)[1]
    # reduce rhizome and root at the end of each year
    initial_state_switchgrass$Rhizome = result$Rhizome[n]*(1-rhizome_loss)
    # initial_state_switchgrass$Root    = result$Root[n]*(1-root_loss)
    # stop()
  }
  #save mean temp
  avg_temp[i] = mean(growing_season$temp)
  #save the aboveground biomass for each location
  aboveground = result$Stem+result$Leaf
  predicted[i]  = tail(aboveground,1)
  season_last_day = as.Date(tail(growing_season$doy,1),origin=paste0(growing_season$year[1],"-01-01"))
  if(season_last_day < harvest_date){
    post_season_loss[i] = post_season_loss_per_day*(as.numeric(harvest_date - season_last_day))
    predicted[i] = predicted[i] - post_season_loss[i]
  }
}

obs_and_model = cbind(obs,predicted,avg_temp)
obs_and_model = as.data.frame(obs_and_model)
colnames(obs_and_model)[colnames(obs_and_model) =="MEAN"] = "observed"
colnames(obs_and_model)[colnames(obs_and_model) =="use_or_not"] = "treatment"
obs_and_model$treatment = as.character(obs_and_model$treatment)
obs_and_model$cultivars=ifelse(obs_and_model$iscaveinrock==1,"cave-in-rock","others")

obs_and_model_plot = obs_and_model
#[obs_and_model$treatment==1,]
obs_and_model_plot = obs_and_model_plot[!is.na(obs_and_model_plot$predicted) & !is.na(obs_and_model_plot$observed),]
cc=cor(obs_and_model_plot$observed,obs_and_model_plot$predicted)

obs_and_model$cultivars=ifelse(obs_and_model$iscaveinrock==1,"cave-in-rock","others")
#write.csv(obs_and_model,'obs_and_model.csv')

# #
# remove young stands
# ages = obs_and_model$harvest_year - obs_and_model$plantyear +1
# obs_and_model = obs_and_model[ages>3,]

correction_term=1
plot1 <- 
  ggplot(obs_and_model_plot, aes(x=observed,y=predicted*correction_term,color = treatment))+
  geom_point()+#colour = "red") + 
  geom_errorbarh(aes(xmin=observed-SE, xmax=observed+SE), height=.2,position=position_dodge(.9)) +
  # ggtitle('All')+
  ylab('predicted biomass (Mg/ha)')+
  xlab('observed biomass (Mg/ha)')+
  ylim(c(5,30)) + xlim(c(5,30))+
  geom_abline(intercept = 0 , slope = 1)+
  theme(text = element_text(size = 16)) 
print(plot1)

# #average by locations [lat, lon]
# unique_latlon = unique(obs_and_model[,c('lat','lon')])
# matrix_plot = data.frame(id = 1:nrow(unique_latlon),lat=NA,lon=NA,observed=NA,predicted=NA,SE=NA)
# for (i in 1:nrow(unique_latlon)){
#   lat_i = unique_latlon$lat[i]
#   lon_i = unique_latlon$lon[i]
#   matrix_plot$lat[i]=lat_i
#   matrix_plot$lon[i]=lon_i
#   obs_and_model_sub = obs_and_model[obs_and_model$lat==lat_i & obs_and_model$lon==lon_i,]
#   matrix_plot$observed[i]  = mean(obs_and_model_sub$observed,na.rm=TRUE)
#   matrix_plot$predicted[i] = mean(obs_and_model_sub$predicted,na.rm=TRUE)
#   matrix_plot$SE[i] = mean(obs_and_model_sub$SE,na.rm=TRUE)
# }
# 
# plot2 <- 
#   ggplot(matrix_plot, aes(x=observed,y=predicted*correction_term))+
#   geom_point()+#colour = "red") + 
#   geom_errorbarh(aes(xmin=observed-SE, xmax=observed+SE), height=.2,position=position_dodge(.9)) +
#   # ggtitle('All')+
#   ylab('predicted biomass (Mg/ha)')+
#   xlab('observed biomass (Mg/ha)')+
#   ylim(c(5,20)) + xlim(c(5,20))+
#   geom_abline(intercept = 0 , slope = 1)+
#   theme(text = element_text(size = 16)) 
# print(plot2)
# cc=cor(matrix_plot$observed,matrix_plot$predicted,use="na.or.complete")


c(mean(obs_and_model_plot$predicted),sd(obs_and_model_plot$predicted))
mean_rmse = rmse(obs_and_model_plot$observed,obs_and_model_plot$predicted) # 
mean_ccc =  epi.ccc(obs_and_model_plot$observed,obs_and_model_plot$predicted)$rho.c$est  #
