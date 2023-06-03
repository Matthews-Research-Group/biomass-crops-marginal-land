library(BioCroMis) #change this to your BioCro library name
library(ggplot2)
library(epiR)
run_cwrfsoilwater=TRUE
#obs data was pre-processed. See "averaging_biomass.R" on how it was done
obs = read.csv('../data/biomass_observation/switchgrass_observation_averaged_plantyear_maturestand.csv')
latlon_id <- read.csv('../data/weather/NASA_data/unique_latlon_and_id.csv')

rmse<-function(obs,pred){
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
parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
x=c(-0.001814,    0.010002,   0.437179,    0.000734,  229.688196)

parameters_list_switchgrass[parameters_to_optimize] = x

parameters_list_switchgrass0 = parameters_list_switchgrass

# initial_state_switchgrass$Rhizome = 14
initial_state_switchgrass$Rhizome = 0.1
initial_state_switchgrass$Root    = 0.1
parameters_list_switchgrass$iSp   = 1.7
if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass =
  initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
}

initial_state_switchgrass0 = initial_state_switchgrass

rhizome_loss =  0.1
root_loss    =  0.33
post_season_loss_per_day = 0.0213 
#t/ha/day: ref: [Over-winter yield decline in Switchgrass and Miscanthus]

predicted = c()
post_season_loss = c(1:dim(obs)[1])*NA
ttc_scaling_factor_all=c()
for (i in 1:dim(obs)[1]){
  print(i)
  initial_state_switchgrass = initial_state_switchgrass0
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
  site.id <- latlon_id$id[latlon_id$lat==obs$lat[i] & latlon_id$lon==obs$lon[i]]
  biocroinput <- read.csv(paste0('../data/weather/NASA_data/BioCroInputs/site_',site.id,'_lowerTransmittance.csv')) 
  
  if(harvest_month<6){
    year_j = harvest_year -1
  }else{year_j = harvest_year}
  for (year_j in plant_year:harvest_year){
    weather_data = biocroinput[which(biocroinput$year==year_j),]
    growing_season <- get_growing_season_climate(weather_data, threshold_temperature = 0)

    #read in CWRF soil water for the initial soil water content and soil type
    soil_data = read.csv(paste0("../data/weather/CWRF_soil_water/site_",site.id,"/cwrf_soilwater_",year_j,".csv"))
    parameters_list_switchgrass$soil_type_indicator = soil_data$soiltype[1]
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
                         derivative_module_names = deriv_modules_list_switchgrass, verbose = FALSE)
    n=dim(result)[1]
    # reduce rhizome and root at the end of each year
    initial_state_switchgrass$Rhizome = result$Rhizome[n]*(1-rhizome_loss)
    # initial_state_switchgrass$Root    = result$Root[n]*(1-root_loss)
  }
  #save the aboveground biomass for each location
  aboveground = result$Stem+result$Leaf #+result$LeafLitter
  predicted[i]  = tail(aboveground,1)
  season_last_day = as.Date(tail(growing_season$doy,1),origin=paste0(growing_season$year[1],"-01-01"))
  if(season_last_day < harvest_date){
    post_season_loss[i] = post_season_loss_per_day*(as.numeric(harvest_date - season_last_day))
    predicted[i] = predicted[i] - post_season_loss[i]
  }
}

obs_and_model = cbind(obs,predicted)
obs_and_model = as.data.frame(obs_and_model)
colnames(obs_and_model)[colnames(obs_and_model) =="MEAN"] = "observed"
obs_and_model$cultivars=ifelse(obs_and_model$iscaveinrock==1,"cave-in-rock","others")
correction_term=1
plot1 <- 
  ggplot(obs_and_model, aes(x=observed,y=predicted*correction_term,color = cultivars))+
  geom_point()+#colour = "red") + 
  geom_errorbarh(aes(xmin=observed-SE, xmax=observed+SE), height=.2,position=position_dodge(.9)) +
  # ggtitle('All')+
  ylab('predicted biomass (Mg/ha)')+
  xlab('observed biomass (Mg/ha)')+
  ylim(c(5,20)) + xlim(c(5,20))+
  geom_abline(intercept = 0 , slope = 1)+
  theme(text = element_text(size = 16)) 
print(plot1)

c(mean(obs_and_model$predicted),sd(obs_and_model$predicted))

obs_and_model = obs_and_model[!is.na(obs_and_model$predicted) & !is.na(obs_and_model$observed),]
cor(obs_and_model$observed,obs_and_model$predicted)

mean_rmse = rmse(obs_and_model$observed,obs_and_model$predicted) # 
mean_ccc =  epi.ccc(obs_and_model$observed,obs_and_model$predicted)$rho.c$est  #
