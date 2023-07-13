#Yufeng He: Dec 12, 2022

library(ncdf4)
library(BioCroMis)
library(lubridate)
library(ggplot2)
library(epiR)
# library(sf)
# library(raster)
library(ggmap)

run_cwrfsoilwater=TRUE

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

# miscanthus_giganteus_logistic_parameters$LeafWS <- 0.5
# miscanthus_giganteus_logistic_parameters$StomataWS <- 0.5
# ss_module_list = miscanthus_giganteus_ss_logistic_modules
# ss_modules     = miscanthus_giganteus_ss_logistic_modules
# ss_modules <- ss_modules[c(-which(ss_module_list=="stomata_water_stress_linear"),
#                            -which(ss_module_list=="leaf_water_stress_exponential"))]
# # ss_modules <- ss_modules[c(-which(ss_module_list=="leaf_water_stress_exponential"))]
# miscanthus_giganteus_ss_logistic_modules = ss_modules

miscanthus_giganteus_logistic_parameters$alpha1 = 0.045    #https://academic.oup.com/jxb/article/68/2/335/2932218#88082389
miscanthus_giganteus_logistic_parameters$TTemr  = 400
parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","alphaStem",
                            "betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")
x=c(-0.000362,    0.269631,    0.609721,
    29.939436,  -39.155449,    0.199284,  -36.468715,    2.891147,   -5.831675)
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

multisite$predicted_stem_actual  = numeric(Num_obs)
multisite$predicted_stem_peak    = numeric(Num_obs)
multisite$predicted_stem_march   = numeric(Num_obs)

years = 2002:2018
TTc_urbana=data.frame(year = years,total_TTC=NA)
weather_urbana = read.csv(paste0(nasa_path,'site_',2,'_2002_2018_LowerTransmittance.csv'))
for (i in 1:length(years)){
  tmp = weather_urbana$temp[weather_urbana$year==years[i]]
  TTc_urbana$total_TTC[i] = ttc_function(tmp)
}
miscanthus_giganteus_logistic_parameters0 = miscanthus_giganteus_logistic_parameters
unique_IDs = unique(multisite$LocationID)
ttc_scaling_factor_all=c()

if(run_cwrfsoilwater){
  miscanthus_giganteus_deriv_logistic_modules = miscanthus_giganteus_deriv_logistic_modules[-3] #remove two_layer_soil_profile
  miscanthus_giganteus_initial_state =
    miscanthus_giganteus_initial_state[names(miscanthus_giganteus_initial_state)!=c('soil_water_content')]
}

rhizome_winter_loss = 0.34# 34% rhizome dies during winter

#create a column for predicted stem (with winter loss)
multisite$predicted_stem_actual = NA

#go through each site
for (i in 1:length(unique_IDs)){
    unique_ID = unique_IDs[i]
    weather_all = read.csv(paste0(nasa_path,'site_',unique_ID,'_2002_2018_LowerTransmittance.csv'))
    
    print(c("processing site",unique_ID))
    
    site_i = multisite[multisite$LocationID==unique_ID,]
    
    #years/duration of simulations
    planted_year <- site_i$planted_year
    measure_year <- site_i$measure_year
    stand_age    <- site_i$age     #this is actual age of measurements
    
    #maximum number of years
    number_of_years <- max(measure_year) - planted_year[1] + 1
    age_all = 1:number_of_years
    
    #biomass for each age at site_i
    biomass_site_i = (1:number_of_years)*NA
    count_index = 1  #first time the loop age reaches the first value of the measured age
    
    #initial rhizome, this should fo inside age loop to account for multiple years of growth
    initialRhizome = 0.3
    #always start with the first year stand and do a multi-year continuous simulation
    for (age in age_all){
      
      year_current = planted_year[1] + age -1
      weather <- weather_all[weather_all$year==year_current,]
      growing_season_weather <- get_growing_season_climate(weather, threshold_temperature = 0)
      
      #read in CWRF soil water for the initial soil water content and soil type
      soil_data = 
        read.csv(paste0("../data/weather_data_NASA_POWER/cwrf_soil_data_1m/site_",unique_ID,"/cwrf_soilwater_",year_current,".csv"))
      miscanthus_giganteus_logistic_parameters$soil_type_indicator = soil_data$soiltype[1]
      
      if(run_cwrfsoilwater){
        growing_season_weather$soil_water_content = soil_data$swc[soil_data$doy>=growing_season_weather$doy[1] 
                                                          & soil_data$doy<=tail(growing_season_weather$doy,1)]
        # growing_season_weather$soil_water_content = growing_season_weather$soil_water_content +0.2
      }
      
      total_ttc = ttc_function(growing_season_weather$temp)
      ttc_scaling_factor = total_ttc / TTc_urbana$total_TTC[TTc_urbana$year==year_current]
      miscanthus_giganteus_logistic_parameters$TTemr = miscanthus_giganteus_logistic_parameters0$TTemr * ttc_scaling_factor
      miscanthus_giganteus_logistic_parameters$TTveg = miscanthus_giganteus_logistic_parameters0$TTveg * ttc_scaling_factor
      miscanthus_giganteus_logistic_parameters$TTrep = miscanthus_giganteus_logistic_parameters0$TTrep * ttc_scaling_factor
      ttc_scaling_factor_all = c(ttc_scaling_factor_all,ttc_scaling_factor)
      
      miscanthus_giganteus_initial_state$Rhizome = initialRhizome
      
      miscanthus_giganteus_logistic_parameters$soil_depth  = site_i$bedrock[1]
      miscanthus_giganteus_logistic_parameters$soil_depth2 = site_i$bedrock[1]/2
      miscanthus_giganteus_logistic_parameters$soil_depth3 = site_i$bedrock[1]
      
      miscanthus_giganteus_logistic_parameters$soil_type_indicator = site_i$biocro_soiltype[1]
  
      result <- Gro_solver(initial_state =  miscanthus_giganteus_initial_state,
                           parameters = miscanthus_giganteus_logistic_parameters,
                           varying_parameters = growing_season_weather,
                           steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
                           derivative_module_names = miscanthus_giganteus_deriv_logistic_modules, verbose = FALSE)
      
      initialRhizome = result$Rhizome[dim(result)[1]]*(1-rhizome_winter_loss)
      print(paste("year=", year_current))
      print(paste("age=", age, "ttc_scaling_factor=",ttc_scaling_factor,"initialRhizome=",initialRhizome))
      print(paste("mean temp=", mean(growing_season_weather$temp)))
      print(paste("mean solar=", mean(growing_season_weather$solar)))
      print(paste("mean stomaWS=", mean(result$StomataWS)))
      print(paste("mean swc=", mean(growing_season_weather$soil_water_content)))
      #calculate winter loss of stem biomass
      #this is only valid when the loop age is matching an actual measured age
      if(age %in% stand_age){
        
        #derive from Fig. 10 Midwest https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12929
        #it makes no difference here since the oldest stand is just 7-year old
        #will use this in the multi-site validation script
          aging_loss = -0.7101*max(7,age) + 4.8145 #at least age 7. Negative values
        
        if(age<=3){
          #for young stands, we apply a constant loss of 33%
          biomass_site_i[age] = result$Stem[dim(result)[1]] * 0.67
        }else{
        #for older stands, we apply a loss of 0.07 t/ha per day
          if(is.na(site_i$measure_month_n[count_index])){ #if no measure month, assume the March 1st (DOY=60)
            site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*(60+365-tail(growing_season_weather$doy,1))
          }else{
            if(site_i$measure_month_n[count_index]>6){ #if measure month is later than June, it suggests harvesting in Winter
              measure_doy = site_i$measure_doy[count_index]
              measure_doy = as.character(measure_doy)
              measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
              number_of_loss_days = measure_doy-tail(growing_season_weather$doy,1)+1
              site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
            }else{              #if measure month is less than June, it suggests harvesting in Spring
              measure_doy = site_i$measure_doy[count_index]
              measure_doy = as.character(measure_doy)
              measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
              number_of_loss_days = 365-tail(growing_season_weather$doy,1)+1+measure_doy
              site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
            }
            if(number_of_loss_days<0) stop(c('age is',age))
            site_i$predicted_stem_actual[count_index]  = site_i$predicted_stem_actual[count_index]
            #* site_i$correction_factor[count_index]
          }
          biomass_site_i[age] = site_i$predicted_stem_actual[count_index] + aging_loss
        }
        count_index = count_index+1
      }else{
        biomass_site_i[age] = NaN
      }
      # if(age==2) stop()
      print(biomass_site_i[age])
    }
    #bind the biomass back to the full multisite matrix for plotting later
    #note that the multisite matrix does not need all ages' records
    multisite$predicted_stem_actual[multisite$LocationID==unique_ID] = biomass_site_i[age_all%in%stand_age]
}

# This correction is required because some yields are becoming negative due to excess loss calculated based on daily loss after reaching peak. We may need to find a better way
# to calculate yield loss in winter using percentage instead of absolute yield loss per day.
multisite_clean = multisite
# multisite_clean = multisite[multisite$predicted_stem_actual>5,]
# multisite_clean = multisite_clean[multisite_clean$biomass>5,] #observed
# multisite_clean = multisite_clean[(multisite_clean$LocationID!=16),]
uniqueID <- unique(multisite_clean$LocationID)
uniqueID = uniqueID[!is.na(uniqueID)]

observed = numeric(length(uniqueID))
obs_se   = numeric(length(uniqueID))
predicted =  numeric(length(uniqueID))
ii=1
for (i in uniqueID){
  print(i)
  tmp = multisite_clean[multisite_clean$LocationID==i,]
  tmp = tmp[!is.na(tmp$biomass),]
  observed[ii] = mean(tmp$biomass)
  obs_se[ii]   = sd(tmp$biomass)/sqrt(length(tmp$biomass))
  predicted[ii] = mean(tmp$predicted_stem_actual)
  ii=ii+1
}

mean_rmse = rmse(observed,predicted) # 
mean_ccc =  epi.ccc(observed,predicted)$rho.c$est  #

cc = cor(observed,predicted)

slope_parameters  = lm(predicted~observed+0) # 
meandataset = data.frame(obs=observed, pred = predicted, SE=obs_se, location=uniqueID)
actual_stem_comparison <- ggplot(data =  meandataset,aes(x = obs ,  y = pred))+ 
  geom_point() +
  stat_smooth(method="lm",formula=y~0+x, level = 0.9) +
  geom_errorbarh(aes(xmin=obs-SE, xmax=obs+SE), height=.2,position=position_dodge(.9)) +
  xlab ("Observed Yield (Mg/ha)") + ylab("Predicted Yield (Mg/ha)") +
  ylim(0,40) + xlim(0,40)  +
  geom_abline() +
  geom_text(aes(label=location),hjust=0,vjust=0)+
  geom_text(x=28 , y = 22 , label = "1:1 line" ) + 
  geom_text(x=10 , y = 30 , label = paste0("Predicted = ",round(slope_parameters$coefficients,2)," × Observed"), color="blue")

# write.csv(meandataset, file="./multisite_validation_data_FigS3b.csv")
plot(actual_stem_comparison)
# ggsave(actual_stem_comparison, filename = "multisite_miscanthus_validation.png", dpi = 500, width =6.5, height =4.2)
