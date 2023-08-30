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
parameters_list_switchgrass$iSp = 1.1
parameters_list_switchgrass$soil_type_indicator = 4

# parameters_list_switchgrass$tbase      = 10
# parameters_list_switchgrass$topt_upper = 31
# parameters_list_switchgrass$tmax       = 40
# parameters_list_switchgrass$upperT       = 28

parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf", "alphaRoot", "betaRoot",
                            "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")

x=c( 12.172727,   -8.076308,   14.594526,  -12.076529,    9.122148,  -13.958066,
     -0.003602,    0.019288,    0.015474,    0.000361,  100.848453)

parameters_list_switchgrass[parameters_to_optimize] = x

initial_state_switchgrass$Rhizome = 10  #mature-stand's Rhizome
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
}

if(run_cwrfsoilwater){
  parameters_list_switchgrass$wsFun=0
  deriv_modules_list_switchgrass = deriv_modules_list_switchgrass[-3]
  initial_state_switchgrass =
    initial_state_switchgrass[names(initial_state_switchgrass)!=c('soil_water_content')]
}

weather_urbana = read.csv('../data/weather/NASA_data/BioCroInputs/site_1_lowerTransmittance.csv')

# rhizome_loss = 0.1
# root_loss =  0.3

years = 2006:2008
result_list=list()
for (i in 1:length(years)){
  ind = which(weather_urbana$year==years[i])
  growing_season=weather_urbana[ind,]
  # growing_season = get_growing_season_climate(growing_season,threshold_temperature = 0)
  growing_season = growing_season[growing_season$doy>=105 & growing_season$doy<=290,]
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

result_bind <- rbind(result2006, result2007, result2008)
result_bind$aboveground <-  result_bind$stem + result_bind$leaf #+ result$leaflitter
predicted <- result_bind[,c("doy","aboveground","rhizome","root")]
predicted <- reshape2::melt(predicted,id.vars = c("doy"),measure.vars = c("aboveground","rhizome","root"))

#Only one year observed
library(ggplot2)
library(lubridate)

observed2 <- observed[month(observed$date)<12,]
observed2$doy = as.numeric(format(observed2$date,"%j"))
observed2$value[observed2$variable=='root'][1]=0.1 #make the first root 0.1

#mean of results
meanresult <- aggregate(predicted$value~predicted$doy+predicted$variable,data=predicted,FUN=mean,na.rm=TRUE)
meanresult_sd <- aggregate(predicted$value~predicted$doy+predicted$variable,data=predicted,FUN=sd,na.rm=TRUE)
meanresult = cbind(meanresult,meanresult_sd$`predicted$value`)
colnames(meanresult) = c('doy','variable','value','sd')

#ggplot
switchgrassplot <- ggplot(data=meanresult) +
             geom_line(aes(x = doy, y = value,color = variable))
switchgrassplot <- switchgrassplot + 
  geom_point(aes(x=doy,y=value, color=variable), data = observed2)
switchgrassplot <- switchgrassplot + xlab("Day of Year") + ylab ("dry biomass (Mg/ha)")
plot(switchgrassplot)

