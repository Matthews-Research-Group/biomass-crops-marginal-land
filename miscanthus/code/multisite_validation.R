#Yufeng He: Dec 12, 2022

library(ncdf4)
library(BioCroMis)
library(lubridate)
library(ggplot2)
library(epiR)
# library(sf)
# library(raster)
library(ggmap)
ttc_function<-function(temp_array){
  tbase = 10
  topt_lower = 28
  topt_upper = 31
  tmax = 40
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
rmse<-function(obs,pred){
  rmse = sqrt(sum((pred-obs)^2)/length(obs))
}
source("output_for_different_harvest_dates.R")

nasa_path = "../data/weather_data_NASA_POWER/BioCro_input_NASA/"
################################ COMMENT 1 STARTS ###############################
# parameters inputs
#################################################################################
load("../data//parameters/miscanthus_giganteus_initial_state.rdata")
load("../data//parameters/miscanthus_giganteus_logistic_parameters.rdata")
load("../data//parameters/miscanthus_giganteus_ss_logistic_modules.rdata")
load("../data//parameters/miscanthus_giganteus_deriv_logistic_modules.rdata")
#################################COMMENT 1 ENDS #################################
# miscanthus_giganteus_logistic_parameters$kd = 0.7         #An Introduction to Environmental Biophysics
# miscanthus_giganteus_logistic_parameters$alpha1 = 0.05    #https://academic.oup.com/jxb/article/68/2/335/2932218#88082389
parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","alphaStem",
                            "betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")
x= c(-0.000351,    0.295862,    0.614672,
     25.991699,  -34.721051,    8.931936,  -33.450919,   22.825293,  -33.522782)
miscanthus_giganteus_logistic_parameters[parameters_to_optimize] = x

multisite <- read.csv("../data/biomass_observation/Miscanthus_Observation_20230529.csv")
multisite <- multisite[!is.na(multisite$biomass),] #remove NA observed biomass

# latlon=multisite[,c('lat','lon')]
# latlon=unique(latlon)
# # us <- readRDS('/Users/yufeng/Desktop/model_coupling_UMD/BioCro/Ob data/GADM_2.8_USA_adm1.rds')
# register_google("AIzaSyBRFL6mFkIl7EFNMBQdpiWKZfjjz-Ncdgg")
# CONUS_map <- get_map(location = c(lon = -96.5795, lat = 39.8283), # center map on CONUS
#                      color = "bw",
#                      source = "google",
#                      maptype = "terrain",
#                      # maptype = "toner-background",
#                      zoom = 4) # set the zoom to include all data
# ggmap(CONUS_map,
#       extent = "device",
#       xlab = "Longitude",
#       ylab = "Latitude") +
#   # Add core lat/long coordinates as points, colored by fraction organic matter
#   geom_point(data = latlon, aes_string(x = 'lon', y = 'lat'), color='blue',size = 1) +
#   xlab("Longitude") + # Change x axis title
#   ylab("Latitude") # Change y axis title

Num_obs <- dim(multisite)[1]

# new column for biocro soil type, exact date of measurement if known, and identification about certainty of date of measurement
# Adding columns for biocro soil type based on CWRF soil type
multisite$biocro_soiltype <- numeric(dim(multisite)[1])
multisite$measure_day_known <- numeric(dim(multisite)[1])
multisite$measure_day_n <- numeric(dim(multisite)[1])
multisite$measure_month_n <- numeric(dim(multisite)[1])
multisite$measure_doy <- numeric(dim(multisite)[1])

multisite$measure_day_known <- ifelse(multisite$measure_month == "November-April", 0,1)
multisite$measure_day_n <- ifelse(is.na(multisite$measure_day), 1, multisite$measure_day)

multisite$measure_month_n <- as.numeric(as.character(multisite$measure_month))

for (i in 1:Num_obs){
  multisite$biocro_soiltype[i] = switch(multisite$soi_type_CWRF[i],0,1,2,4,4,3,5,7,6,8,9,10)
  if(multisite$measure_day_known[i]==0){
    multisite$measure_doy[i] = "NA"
  } else {
    tmp1 <-  paste(
             sprintf("%04d",multisite$measure_year[i]),
             sprintf("%02d",multisite$measure_month_n[i]),
             sprintf("%02d",multisite$measure_day_n[i]),
             sep="")
     multisite$measure_doy[i] = tmp1
  }
}


multisite$predicted_stem_actual  = numeric(Num_obs)
multisite$predicted_stem_peak    = numeric(Num_obs)
multisite$predicted_stem_march   = numeric(Num_obs)

years = 2002:2018
TTc_urbana=data.frame(year = years,total_TTC=NA)
weather_urbana = read.csv(paste0(nasa_path,'site_',2,'_2002_2018.csv'))
for (i in 1:length(years)){
  tmp = weather_urbana$temp[weather_urbana$year==years[i]]
  TTc_urbana$total_TTC[i] = ttc_function(tmp)
}
miscanthus_giganteus_logistic_parameters0 = miscanthus_giganteus_logistic_parameters
ttc_scaling_factor_all=c()
for (i in 1:Num_obs){
    print(i)
    #years/duration of simulations
    planted.year <- multisite$planted_year[i]
    measure.year <- multisite$measure_year[i]
    number_of_years <- measure.year - planted.year
    
    #local parameters for simulations
    lat_obs = multisite$lat[i]
    lon_obs = multisite$lon[i]
    siteID = multisite$LocationID[i]
    weather_all = read.csv(paste0(nasa_path,'site_',siteID,'_2002_2018.csv'))
  
    #initial rhizome, this should fo inside age loop to account for multiple years of growth
    initialRhizome = 0.3
    
    for (age in 1:number_of_years){
      year_current = planted.year + age -1
      weather <- weather_all[weather_all$year==year_current,]
      growing_season_weather <- get_growing_season_climate(weather, threshold_temperature = 0)
      total_ttc = ttc_function(growing_season_weather$temp)
      ttc_scaling_factor = total_ttc / TTc_urbana$total_TTC[TTc_urbana$year==year_current]
      miscanthus_giganteus_logistic_parameters$TTemr = miscanthus_giganteus_logistic_parameters0$TTemr * ttc_scaling_factor
      miscanthus_giganteus_logistic_parameters$TTveg = miscanthus_giganteus_logistic_parameters0$TTveg * ttc_scaling_factor
      miscanthus_giganteus_logistic_parameters$TTrep = miscanthus_giganteus_logistic_parameters0$TTrep * ttc_scaling_factor
      ttc_scaling_factor_all = c(ttc_scaling_factor_all,ttc_scaling_factor)
      
      miscanthus_giganteus_initial_state$Rhizome = initialRhizome
      
      miscanthus_giganteus_logistic_parameters$soil_depth  = multisite$bedrock[i]
      miscanthus_giganteus_logistic_parameters$soil_depth2 = multisite$bedrock[i]/2
      miscanthus_giganteus_logistic_parameters$soil_depth3 = multisite$bedrock[i]
      
      miscanthus_giganteus_logistic_parameters$soil_type_indicator = multisite$biocro_soiltype[i]
  
      result <- Gro_solver(initial_state =  miscanthus_giganteus_initial_state,
                           parameters = miscanthus_giganteus_logistic_parameters,
                           varying_parameters = growing_season_weather,
                           steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
                           derivative_module_names = miscanthus_giganteus_deriv_logistic_modules, verbose = FALSE)
      
      initialRhizome = result$Rhizome[dim(result)[1]]*0.66 # 34% rhizome dies during winter
      print(paste("age=", age, "ttc_scaling_factor=",ttc_scaling_factor))
    }
    if(is.na(multisite$measure_month_n[i])){
      multisite$predicted_stem_actual[i] = result$Stem[dim(result)[1]] - 0.07*(60+365-tail(growing_season_weather$doy,1))
    }else{
      if(multisite$measure_month_n[i]>6){
        measure_doy = multisite$measure_doy[i]
        measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
        number_of_loss_days = measure_doy-tail(growing_season_weather$doy,1)+1
        multisite$predicted_stem_actual[i] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
      }else{
        measure_doy = multisite$measure_doy[i]
        measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
        number_of_loss_days = 365-tail(growing_season_weather$doy,1)+1+measure_doy
        multisite$predicted_stem_actual[i] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
      }
      if(number_of_loss_days<0) stop(c('i is',i))
      multisite$predicted_stem_actual[i]  = multisite$predicted_stem_actual[i] * multisite$correction_factor[i]
    }
    # harvest_option <-  output_for_different_harvest_dates(result,multisite$measure_day_n[i],
    #                                                        multisite$measure_month_n[i],
    #                                                        year_current,
    #                                                        multisite$measure_day_known[i])
    # multisite$predicted_stem_actual[i]  = harvest_option$actual_harvested_biomass * multisite$correction_factor[i]
    # multisite$predicted_stem_peak[i]    = harvest_option$peak_harvested_biomass* multisite$correction_factor[i]
    # multisite$predicted_stem_march[i]   = harvest_option$march15_harvested_biomass * multisite$correction_factor[i]
    # print(harvest_option$actual_harvested_biomass)
    # print(harvest_option$peak_harvested_biomass)
    # print(harvest_option$march15_harvested_biomass)
    print(multisite$predicted_stem_actual[i])
}

# This correction is required because some yields are becoming negative due to excess loss calculated based on daily loss after reaching peak. We may need to find a better way
# to calculate yield loss in winter using percentage instead of absolute yield loss per day.
multisite_clean <- multisite[multisite$predicted_stem_actual>5,]
multisite_clean <- multisite_clean[multisite_clean$biomass>5,] #observed
multisite_clean = multisite_clean[(multisite_clean$LocationID!=16),]
# library(epiR)
# library(Metrics)
uniqueID <- unique(multisite_clean$LocationID)


observed = numeric(length(uniqueID))
obs_se   = numeric(length(uniqueID))
predicted =  numeric(length(uniqueID))
ii=1
for (i in uniqueID){
  print(i)
  tmp <- multisite_clean[multisite_clean$LocationID==i,]
  tmp <- tmp[!is.na(tmp$biomass),]
  observed[ii] = mean(tmp$biomass)
  obs_se[ii]   = sd(tmp$biomass)/sqrt(length(tmp$biomass))
  predicted[ii] = mean(tmp$predicted_stem_actual)
  ii=ii+1
}
mean_rmse = rmse(observed,predicted) # 
mean_ccc =  epi.ccc(observed,predicted)$rho.c$est  #

cor(observed,predicted)

slope_parameters  = lm(predicted~observed+0) # 
meandataset = data.frame(obs=observed, pred = predicted, SE=obs_se, location=uniqueID)
actual_stem_comparison <- ggplot(data =  meandataset,aes(x = obs ,  y = pred))+ 
  geom_point() +
  # stat_smooth(method="lm",formula=y~0+x, level = 0.9) +
  geom_errorbarh(aes(xmin=obs-SE, xmax=obs+SE), height=.2,position=position_dodge(.9)) +
  xlab ("Observed Yield (Mg/ha)") + ylab("Predicted Yield(Mg/ha)") +
  ylim(0,35) + xlim(0,35)  +
  geom_abline() 
#  geom_text(x=28 , y = 22 , label = "1:1 line" ) + 
  # geom_text(x=10 , y = 30 , label = paste0("Predicted = ",round(slope_parameters$coefficients,2)," × Observed"), color="blue")

# write.csv(meandataset, file="./multisite_validation_data_FigS3b.csv")
plot(actual_stem_comparison)
# ggsave(actual_stem_comparison, filename = "multisite_miscanthus_validation.png", dpi = 500, width =6.5, height =4.2)
