library(BioCroMis) #requires installation
source("set_model_parameters.R")
crop_type = "Miscanthus" #Three crops: Miscanthus, Switchgrass, Energycane
#prepare biocro inputs
biocro_module_list <- set_modules_parameters(crop_type)
#read in weather data
weather_all<-read.csv(
  "../calibrations/miscanthus/data/weather_data_NASA_POWER/BioCro_input_NASA/site_2_2002_2018_LowerTransmittance.csv")

growing_season_weather = weather_all[weather_all$year == 2008,]
#subset the full weather to the growing season
#the example here is for typical Miscanthus growth
#other crops may have very different growing seasons
growing_season_weather =  growing_season_weather[growing_season_weather$doy>=106 &
                                                    growing_season_weather$doy<=350,]

result <- Gro_solver(initial_state =  biocro_module_list[[1]],
                     parameters = biocro_module_list[[2]],
                     varying_parameters = growing_season_weather,
                     steady_state_module_names = biocro_module_list[[3]],
                     derivative_module_names = biocro_module_list[[4]], 
                     solver=biocro_module_list[[5]],verbose = FALSE)

plot(result$doy,result$Stem)
