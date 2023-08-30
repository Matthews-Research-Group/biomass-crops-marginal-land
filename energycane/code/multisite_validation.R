library(stringr)
library(BioCroMis)
library(gridExtra)
library(epiR)
#Energy cane observations
rmse<-function(obs,pred){
  NAs  = which(is.na(obs) | is.na(pred))
  pred = pred[-NAs]
  obs  = obs[-NAs]
  rmse = sqrt(sum((pred-obs)^2)/length(obs))
}
observations <- read.csv("../data/observations/observation_averaged_r1.csv")
row_ind      <- readRDS('../data/observations/row_index_avg_r1.rds')
observations$biocro_stem = NA
observations$biocro_aboveground <- NA
observations$Tmean <- NA
observations$rainfall <- NA


#Parameters for running the energycane simulations
initial_state      <- readRDS("../data/parameters//initial_state.rds")
parameters_list    <- readRDS("../data/parameters//parameters_list.rds")
ss_module_list     <- readRDS("../data/parameters//ss_module_list.rds")
deriv_modules_list <- readRDS("../data/parameters//deriv_modules_list.rds")

parameters_list$TTemr = 900
parameters_list$TTc_leafsenescence_threshold = 50
parameters_list$alphaStem = 8.202316
parameters_list$betaStem = -0.28125
parameters_list$alphaLeaf = 7.565351
parameters_list$betaLeaf = -7.949926
parameters_list$alphaRoot = 2.372882
parameters_list$betaRoot = -39.87503
parameters_list$leaf_turnover_rate = 0.0005305176
parameters_list$tbase = 13
ss_module_list[13] = "stomata_water_stress_linear"
ss_module_list[14] = "leaf_water_stress_exponential"

#Number of rows
N = dim(observations)[1]
result.all <- list()
finalindex.all <- NULL
growing_period <- data.frame(matrix(ncol=2, nrow=N))
colnames(growing_period) <- c("start","end")

for (i in 1:N){
  #Find planting and harvest measurement years
  # plantingyear <- as.numeric(str_sub(observations$Date.of.Planting[i],-4,-1))
  # harvestyear <- as.numeric(str_sub(observations$Date.of.Measurement[i],-4,-1))
  rows = row_ind[[i]]
  weather_file_index = rows[1]
  # print(weather_file_index)
  #Find latitude and longitude
  latitude  = observations$Lat[i]
  longitude = observations$Lon[i]
  '
  Currently, I am just going to use default weather data using criteria ( T>0) to determine
  growing period. The weather data will correspond to harvest/measures year
  This approach can be later refined using actual date of planting and harvest
  '
  #Location of weather file
  growingseason.name <- paste("../data/weather/",weather_file_index,"/","growing_season",".csv",sep="")
  growing_year <- read.csv(growingseason.name)
  growing_season <- get_growing_season_climate(growing_year, threshold_temperature = 0)
  
  growing_period$start[i] <- growing_season$doy[1]
  growing_period$end[i] <- growing_season$doy[nrow(growing_season)]
  
  parameters_list$lat <- latitude
  
  params <- parameters_list
  ss_modules <- ss_module_list
  if(observations$treatment[i]=='Irrigated'){
    params$LeafWS <- 1
    params$StomataWS <- 1
    ss_modules <- ss_modules[c(-which(ss_module_list=="stomata_water_stress_linear"),
                               -which(ss_module_list=="leaf_water_stress_exponential"))]
  }
  
  result <- Gro_solver(initial_state =  initial_state,
                       parameters = params,
                       varying_parameters = growing_season,
                       steady_state_module_names = ss_modules,
                       derivative_module_names = deriv_modules_list, verbose = FALSE)
  
  result.all[[i]] <- result
  
  finalindex <- dim(result)[1]
  finalindex.all[i] <- finalindex
  
  observations$biocro_stem[i] <- result$Stem[finalindex]
  observations$biocro_aboveground[i] <- result$Stem[finalindex] + result$Leaf[finalindex] + 
                              result$LeafLitter[finalindex]
  observations$Tmean[i] = mean(growing_season$temp)
  observations$rainfall[i] = mean(growing_season$precip)
  if(i==1){
    test_growing_season <- cbind (ID=i,growing_season)
  } else {
    test_growing_season <- rbind(test_growing_season,cbind (ID=i,growing_season))
  }
}

####################################
correction_term = 0.91^(observations$Ratoon-1)

predicted <- ifelse(observations$YieldID=="Stem",observations$biocro_stem, observations$biocro_aboveground)
predicted = predicted * correction_term
validation_data_for_ggplot <- data.frame(observed=observations$MEAN,predicted=predicted,ID=observations$YieldID,
                                         citation_author = observations$citation_author, treatment=observations$treatment,
                                         SD=observations$SD,Ratoon = observations$Ratoon)
validation_data_for_ggplot$Ratoon=paste("Ratoon ",validation_data_for_ggplot$Ratoon, sep="")
# validation_data_for_ggplot <- validation_data_for_ggplot[validation_data_for_ggplot$citation_author == "Knoll",]
library(ggplot2)
energycane_validation <- 
  ggplot(validation_data_for_ggplot, aes(x=observed,y=predicted, color = Ratoon))+
  geom_point() + 
  geom_errorbarh(aes(xmin=observed-SD, xmax=observed+SD), height=.2,position=position_dodge(.9)) +
  xlab ("Observed Yield (Mg/ha)") + ylab("Predicted Yield(Mg/ha)") +
  ggtitle('All')+
  ylim(c(0,80)) + xlim(c(0,80))+
  geom_abline(intercept =0 , slope = 1)
plot(energycane_validation)
cor(validation_data_for_ggplot$observed,validation_data_for_ggplot$predicted)

mean_rmse = rmse(validation_data_for_ggplot$observed,validation_data_for_ggplot$predicted) #
mean_ccc =  epi.ccc(validation_data_for_ggplot$observed,validation_data_for_ggplot$predicted)$rho.c$est  #

stop()
validation.rainfed <- validation_data_for_ggplot[which(validation_data_for_ggplot$treatment=='Rainfed'),]
rainfed <- 
  ggplot(validation.rainfed, aes(x=observed,y=predicted, color = Ratoon))+
  geom_point() + 
  geom_errorbarh(aes(xmin=observed-SD, xmax=observed+SD), height=.2,position=position_dodge(.9)) +
  ggtitle('rainfed')+
  ylim(c(0,80)) + xlim(c(0,80))+
  geom_abline(intercept =0 , slope = 1)
# plot(rainfed)
cor(validation.rainfed$observed,validation.rainfed$predicted)

validation.irrigated <- validation_data_for_ggplot[which(validation_data_for_ggplot$treatment=='Irrigated'),]
irri <- 
  ggplot(validation.irrigated, aes(x=observed,y=predicted, color = Ratoon))+
  geom_point() + 
  geom_errorbarh(aes(xmin=observed-SD, xmax=observed+SD), height=.2,position=position_dodge(.9)) +
  ggtitle('Irrigated')+
  ylim(c(0,80)) + xlim(c(0,80))+
  geom_abline(intercept =0 , slope = 1)
# plot(irri)
cor(validation.irrigated$observed,validation.irrigated$predicted)

grid.arrange(energycane_validation, rainfed, irri, ncol=2, nrow =2)
