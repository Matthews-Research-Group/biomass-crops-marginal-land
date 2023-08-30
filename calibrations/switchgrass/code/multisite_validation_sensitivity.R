library(BioCroMis) #change this to your BioCro library name
library(ggplot2)
library(epiR)
find_identical_rows<-function(A,x){
  nr = nrow(A)
  output=c()
  for (i in 1:nr){
    a <- A[i,]
    output = c(output,isTRUE(all.equal(a,x,check.attributes=FALSE)))
  }
  return(output)
}

run_cwrfsoilwater=TRUE
#obs data was pre-processed. See "averaging_biomass.R" on how it was done
obs = read.csv('../data/biomass_observation/switchgrass_observation_averaged_plantyear_maturestand_r2.csv')
latlon_id <- read.csv('../data/weather/NASA_data/unique_latlon_and_id.csv')
X = obs[,c('lat','lon','plantyear','use_or_not')]
unique_obs = unique(X)
row_ind_unique = list()
for (i in 1:nrow(unique_obs)){
  row_ind_unique[[i]]=which(find_identical_rows(X,unique_obs[i,])==TRUE)
}

rmse<-function(obs,pred){
  pred = pred[!is.na(obs)]
  obs  = obs[!is.na(obs)]
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
parameters_list_switchgrass$iSp   = 1.7

# parameters_list_switchgrass$tbase      = 4   #0-12
# parameters_list_switchgrass$topt_upper = 38  #20-35
# parameters_list_switchgrass$tmax       = 45  #36-45

parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf", "alphaRoot", "betaRoot",
                            "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
x=c(10.019875,   -5.734288,    5.418285,   -0.490853,   14.444107,  -17.618838,
    -0.001025,    0.143396,    0.615253,    0.000096,  253.538347)

parameters_list_switchgrass[parameters_to_optimize] = x

if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass =
  initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
}
initial_state_switchgrass0 = initial_state_switchgrass

upperT_vec = seq(25,40,by=2)
lowerT_vec = 6

obs_and_model_list=list()
para_grid = expand.grid(upperT_vec,lowerT_vec)
colnames(para_grid) = c('upperT','lowerT')

rhizome_loss =  0.33
# root_loss    =  0.33
post_season_loss_per_day = 0.0213 #t/ha/day
#ref paper: [Over-winter yield decline in Switchgrass and Miscanthus]

#[not used] This paper suggests loss of 35% to 45% due to machine harvest
#https://www.sciencedirect.com/science/article/pii/S0961953409000166
# 

for (i in 1:length(row_ind_unique)){
  row_ind = row_ind_unique[[i]]
  obs_sub = obs[row_ind,]
  lat_i = obs_sub$lat[1]
  lon_i = obs_sub$lon[1]
  site.id <- latlon_id$id[latlon_id$lat==lat_i & latlon_id$lon==lon_i]
  biocroinput = read.csv(paste0('../data/weather/NASA_data/BioCroInputs/site_',site.id,'_lowerTransmittance.csv')) 

  #re-init for every location
  initial_state_switchgrass = initial_state_switchgrass0
  #
  parameters_list_switchgrass$lat=lat_i
  plant_year     = obs_sub$plantyear[1]
  harvest_year   = max(obs_sub$harvest_year)
  harvest_months = obs_sub$harvest_month

  harvest_year_actual = harvest_year
  if(mean(harvest_months)<6) harvest_year_actual = harvest_year -1
  print(c("plant",plant_year,"harvest",harvest_year_actual))

  output = cbind(para_grid,rmse=NA,cc=NA)
  for (j in 1:nrow(para_grid)){
  
    parameters_list_switchgrass$upperT   = para_grid$upperT[j]
    parameters_list_switchgrass$lowerT   = para_grid$lowerT[j]
    
    parameters_list_switchgrass0 = parameters_list_switchgrass
    
    initial_state_switchgrass$Rhizome = 0.1
    initial_state_switchgrass$Root    = 0.1
    parameters_list_switchgrass$soil_type_indicator = 4
    solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
             adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)
    #http://dx.doi.org/10.1080/07352680500316433
    #https://doi.org/10.2135/cropsci2004.1391

    predict = c() 
    for (year_k in plant_year:harvest_year_actual){
      weather_data = biocroinput[which(biocroinput$year==year_k),]
      growing_season <- get_growing_season_climate(weather_data, threshold_temperature = 0)
      #read in CWRF soil water for the initial soil water content and soil type
      soil_data = read.csv(paste0("../data/weather/CWRF_soil_water/site_",site.id,"/cwrf_soilwater_",year_k,".csv"))
      if(run_cwrfsoilwater){
        growing_season$soil_water_content = soil_data$swc[soil_data$doy>=growing_season$doy[1] 
                                                        & soil_data$doy<=tail(growing_season$doy,1)]
      }else{
        initial_state_switchgrass$soil_water_content = soil_data$swc[soil_data$doy==growing_season$doy[1] & soil_data$hour==0]
      }
      
      total_ttc = ttc_function(growing_season$temp)
      ttc_scaling_factor = total_ttc / TTc_urbana$total_TTC[TTc_urbana$year==year_k]
      parameters_list_switchgrass$TTemr = parameters_list_switchgrass0$TTemr * ttc_scaling_factor
      parameters_list_switchgrass$TTveg = parameters_list_switchgrass0$TTveg * ttc_scaling_factor
      parameters_list_switchgrass$TTrep = parameters_list_switchgrass0$TTrep * ttc_scaling_factor

      result <- Gro_solver(initial_state =  initial_state_switchgrass,
                         parameters = parameters_list_switchgrass,
                         varying_parameters = growing_season,
                         steady_state_module_names = ss_modules_list_swithcgrass,
                         derivative_module_names = deriv_modules_list_switchgrass, solver=solver,verbose = FALSE)
      # reduce rhizome and root at the end of each year
      initial_state_switchgrass$Rhizome = tail(result$Rhizome,1)*(1-rhizome_loss)
      
      #if this year_k is in the obs data, we save the predicted results
      #otherwise, just go to the next year run
      if(mean(harvest_months)<6){
        #if havested in spring, we consider the previous year
        ind = which((obs_sub$harvest_year-1) == year_k)
      }else{
        ind = which(obs_sub$harvest_year == year_k)
      }
      if(length(ind)==1){
        harvest_month  = obs_sub$harvest_month[ind]
        harvest_day    = obs_sub$harvest_day[ind]
        harvest_year_j = obs_sub$harvest_year[ind]
        harvest_date   = as.Date(paste(harvest_year_j,harvest_month,harvest_day,sep="-"))

        aboveground  = result$Stem+result$Leaf
        predict_tmp  = tail(aboveground,1)
        season_last_day = as.Date(tail(growing_season$doy,1),origin=paste0(growing_season$year[1],"-01-01"))
        if(season_last_day < harvest_date){
          post_season_loss = post_season_loss_per_day*(as.numeric(harvest_date - season_last_day))
          predict_tmp = predict_tmp - post_season_loss
        }
        predict = c(predict,predict_tmp)
      }

    } #end for (year_k in plant_year:harvest_year) 
    output$rmse[j] = rmse(obs_sub$MEAN,predict) 
    output$cc[j]   = cor(obs_sub$MEAN,predict,use="complete.obs") 
  } #end for (j in 1:nrow(para_grid))
  print(output)
  stop()
}
