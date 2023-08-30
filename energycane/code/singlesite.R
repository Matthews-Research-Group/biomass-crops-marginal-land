# Biocro build from the branch "direct_radiation_input" of biocro-dev
library(BioCroMis)

# parameters and weather inputs
observed_biomass <- read.csv('../data/observations/observation_for_optimization_correction7.csv')
growing_season <- read.csv("../data/weather/weather_NASA_calibration_LaBelle.csv")

#Parameters for running the energycane simulations
initial_state       <- readRDS("../data/parameters//initial_state.rds")
parameters_list     <- readRDS("../data/parameters//parameters_list.rds")
ss_modules_list     <- readRDS("../data/parameters//ss_module_list.rds")
deriv_modules_list  <- readRDS("../data/parameters//deriv_modules_list.rds")

# Here I am replacing initial parameter values by the optimized values
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
parameters_list$LeafWS = NULL
parameters_list$StomataWS = NULL
ss_modules_list[13] = "stomata_water_stress_linear"
ss_modules_list[14] = "leaf_water_stress_exponential"

result <- Gro_solver(initial_state =  initial_state,
                     parameters = parameters_list,
                     varying_parameters = growing_season,
                     steady_state_module_names = ss_modules_list,
                     derivative_module_names = deriv_modules_list, verbose = FALSE)


library(ggplot2)
source("Compare_Observed_and_Predicted.R")

observed <- observed_biomass[, c("doy","Stem","Leaf","Root")]
observed <- reshape2::melt(observed,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root"))

result =result[result$year==2012,]
predicted <- result[,c("doy","Stem","Leaf","Root")]
predicted <- reshape2::melt(predicted,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root"))
names(observed)<- c("x","varname","y")
names(predicted) <- c("x","varname","y")
xlabtitle=expression(paste("DOY"))
ylabtitle =expression(paste("Dry-biomass in Mg ha"^{-1}))
BiomassPartitioning <- Compare_Observed_and_Predicted(observeddata = observed,predicteddata=predicted  ,xlabtitle,ylabtitle)
plot(BiomassPartitioning)
