# Biocro build from the branch "direct_radiation_input" of biocro-dev
library(BioCroMis)

run_cwrfsoilwater=TRUE
# parameters and weather inputs

observed_biomass <- read.csv('../data/biomass_observation/observations_switchgrass_2006_2008.csv')
observed_biomass2006 <- observed_biomass 
observed_biomass2006$date <- c("2006-04-15","2006-06-15","2006-08-15","2006-10-15","2006-12-15")

observed_biomass2007 <- observed_biomass 
observed_biomass2007$date <- c("2007-04-15","2007-06-15","2007-08-15","2007-10-15","2007-12-15")

observed_biomass2008 <- observed_biomass 
observed_biomass2008$date <- c("2008-04-15","2008-06-15","2008-08-15","2008-10-15","2008-12-15")

observed_biomass <- rbind(observed_biomass2006, observed_biomass2007, observed_biomass2008)
observed_biomass <- observed_biomass2008[, c("date", "Aboveground", "Rhizome", "Root")]  #
names(observed_biomass) <- c("date", "aboveground","rhizome","root")
observed <- reshape2::melt(observed_biomass, id.vars = c("date"),measure.vars = c("aboveground","rhizome","root"))
observed$date <- as.Date(observed$date)

source("../data/parameters/initial_state_switchgrass.R")
source("../data/parameters/ss_modules_list_switchgrass.R")
source("../data/parameters/deriv_modules_list_with_switchgrass_senescence.R")
source("../data/parameters/parameters_list_switchgrass_senescence.R")

if(TRUE){
parameters_list_switchgrass$TTc_leafsenescence_threshold = 5
parameters_list_switchgrass$iSp = 1.7
parameters_list_switchgrass$soil_type_indicator = 4

# parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf", "alphaRoot", "betaRoot",
#                             "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate")
parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")
x=c(-0.001814,    0.010002,   0.437179,    0.000734,  229.688196)

parameters_list_switchgrass[parameters_to_optimize] = x

initial_state_switchgrass$Rhizome = 14  #mature-stand's Rhizome
initial_state_switchgrass$Root    = 0.1
# #remove water stress
# initial_state_switchgrass = 
# initial_state_switchgrass[names(initial_state_switchgrass)!=c("soil_water_content")]
# parameters_list_switchgrass$LeafWS <- 0.5
# parameters_list_switchgrass$StomataWS <- 0.5
# ss_module_list = ss_modules_list_swithcgrass
# ss_modules     = ss_modules_list_swithcgrass
# ss_modules <- ss_modules[c(-which(ss_module_list=="stomata_water_stress_linear"),
#                            -which(ss_module_list=="leaf_water_stress_exponential"))]
# ss_modules_list_swithcgrass = ss_modules

parameters_list_switchgrass$iSp   = 1.7
}

if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass =
    initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
}

######multi year simulations####################
climate2002 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2002.csv")
climate2002 <- get_growing_season_climate(climate2002, threshold_temperature = 0)
climate2003 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2003.csv")
climate2003 <- get_growing_season_climate(climate2003, threshold_temperature = 0)
climate2004 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2004.csv")
climate2004 <- get_growing_season_climate(climate2004, threshold_temperature = 0)
climate2005 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2005.csv")
climate2005 <- get_growing_season_climate(climate2005, threshold_temperature = 0)
climate2006 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2006.csv")
climate2006 <- get_growing_season_climate(climate2006, threshold_temperature = 0)
climate2007 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2007.csv")
climate2007 <- get_growing_season_climate(climate2007, threshold_temperature = 0)
climate2008 <- read.csv("../data/weather/UIUC_weatherinputs_to_BioCro_2008.csv")
climate2008 <- get_growing_season_climate(climate2008, threshold_temperature = 0)
climate <- list(climate2002=climate2002,
                climate2003=climate2003,
                climate2004=climate2004,
                climate2005=climate2005,
                climate2006=climate2006,
                climate2007=climate2007,
                climate2008=climate2008)

# rhizome_loss = 0.1
# root_loss =  0.3
# initial_state_switchgrass$Rhizome = 0.1

years_all = 2002:2008
years     = 2006:2008
result_list=list()
for (i in 1:length(years)){
  ind = which(years_all==years[i])
  growing_season=climate[[ind]]
  soil_data = read.csv(paste0("../data/weather/CWRF_soil_water/site_",1,"/cwrf_soilwater_",years[i],".csv"))
  if(run_cwrfsoilwater){
    growing_season$soil_water_content = soil_data$swc[soil_data$doy>=growing_season$doy[1] 
                                                      & soil_data$doy<=tail(growing_season$doy,1)]
  }
  result <- Gro_solver(initial_state =  initial_state_switchgrass,
                       parameters = parameters_list_switchgrass,
                       varying_parameters = growing_season,
                       steady_state_module_names = ss_modules_list_swithcgrass,
                       derivative_module_names = deriv_modules_list_switchgrass, verbose = FALSE)
  result_list[[i]] = result

  #for continous interannual run, using the previous year's rhizome
  # initial_state_switchgrass$Rhizome = tail(result$Rhizome,1)
  
  if(years[i]==2006){
    result2006 <- data.frame(doy=result$doy,TT = result$TTc,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter,
                             rhizome = result$Rhizome, root = result$Root)
  }else if(years[i]==2007){
    result2007 <- data.frame(doy=result$doy,TT = result$TTc,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter,
                             rhizome = result$Rhizome, root = result$Root)
  }else if(years[i]==2008){
    result2008 <- data.frame(doy=result$doy,TT = result$TTc,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter,
                             rhizome = result$Rhizome, root = result$Root)
  }else if(years[i]==2002){
    result2002 <- data.frame(doy=result$doy,TT = result$TTc,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter,
                             rhizome = result$Rhizome, root = result$Root)
  }
  
}

result2006$date <- as.Date(result2006$doy,origin = "2006-01-01")
result2007$date <- as.Date(result2007$doy,origin = "2007-01-01")
result2008$date <- as.Date(result2008$doy,origin = "2008-01-01")

#repeat 2006's last day
result2006 = rbind(result2006,tail(result2006,24))
result2006$date[(nrow(result2006)-23):nrow(result2006)]=result2006$date[nrow(result2006)]+1

result <- rbind(result2006, result2007, result2008)
result$aboveground <-  result$stem + result$leaf #+ result$leaflitter
predicted <- result[,c("date","aboveground","rhizome","root")]
predicted <- reshape2::melt(predicted,id.vars = c("date"),measure.vars = c("aboveground","rhizome","root"))

#Only one year observed
library(ggplot2)
library(lubridate)

#mean of results
meanresult <- do.call(rbind, list(result2006, result2007,result2008))
meanresult$daymonth <- substr(as.character(meanresult$date),start = 6, stop = 10)

meanresult <- aggregate(meanresult,list(meanresult$daymonth),FUN=mean,na.rm=TRUE)
# meanresult = meanresult[2:(nrow(meanresult)-1),]

meanresult$aboveground <- meanresult$stem + meanresult$leaf #+ meanresult$leaflitter
meanresult$date = as.Date(paste0("2008-",meanresult$Group.1))

observed2 <- observed[month(observed$date)<12,]
observed2$doy = as.numeric(format(observed2$date,"%j"))
meanresult = meanresult[meanresult$date>=observed2$date[1] & meanresult$date <= tail(observed2$date,1),]

predicted <- meanresult[,c("date","doy","aboveground","rhizome","root")]
predicted <- reshape2::melt(predicted,id.vars = c("doy"),measure.vars = c("aboveground","rhizome","root"))
switchgrassplot <- ggplot(data=predicted) + geom_line(aes(x = doy, y = value,color = variable))
switchgrassplot <- switchgrassplot + geom_point(aes(x=doy,y=value, color=variable), data = observed2)
switchgrassplot <- switchgrassplot + xlab("Day of Year") + ylab ("dry biomass (Mg/ha)")
plot(switchgrassplot)



