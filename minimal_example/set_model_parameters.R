set_modules_parameters<-function(crop_type){
   if(crop_type=="Miscanthus"){
     load("../calibrations/miscanthus//data//parameters/miscanthus_giganteus_initial_state.rdata")
     load("../calibrations/miscanthus//data//parameters/miscanthus_giganteus_logistic_parameters.rdata")
     load("../calibrations/miscanthus//data//parameters/miscanthus_giganteus_ss_logistic_modules.rdata")
     load("../calibrations/miscanthus//data//parameters/miscanthus_giganteus_deriv_logistic_modules.rdata")
     
     miscanthus_giganteus_logistic_parameters$TTemr  = 400
     
     miscanthus_giganteus_logistic_parameters$alpha1 = 0.045 #0.04-0.05
     
     parameters_to_optimize <- c("kRhizome_emr","kLeaf_emr","kStem_emr",
                                 "alphaStem","betaStem","alphaLeaf","betaLeaf","alphaRoot", "betaRoot")
     #this one is for turning OFF cwrf soil water
     x=c(-0.000362,    0.269631,    0.609721,
         29.939436,  -39.155449,    0.199284,  -36.468715,    2.891147,   -5.831675)
     miscanthus_giganteus_logistic_parameters[parameters_to_optimize] = x
     miscanthus_giganteus_initial_state$Rhizome = 24.1

     solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
                 adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)
     output_list = list(miscanthus_giganteus_initial_state,
                        miscanthus_giganteus_logistic_parameters,
                        miscanthus_giganteus_ss_logistic_modules,
			miscanthus_giganteus_deriv_logistic_modules,solver)

   }else if(crop_type=="Switchgrass"){
     source("../calibrations/switchgrass//data/parameters/initial_state_switchgrass.R")
     source("../calibrations/switchgrass//data/parameters/ss_modules_list_switchgrass.R")
     source("../calibrations/switchgrass//data/parameters/deriv_modules_list_with_switchgrass_senescence.R")
     source("../calibrations/switchgrass//data/parameters/parameters_list_switchgrass_senescence.R")
     parameters_list_switchgrass$TTc_leafsenescence_threshold = 5
     parameters_list_switchgrass$iSp = 1.1
     parameters_list_switchgrass$soil_type_indicator = 4

     parameters_to_optimize <- c("alphaStem","betaStem","alphaLeaf","betaLeaf", "alphaRoot", "betaRoot",
                            "kRhizome_emr","kLeaf_emr","kStem_emr","leaf_turnover_rate","TTemr")

     x=c( 12.172727,   -8.076308,   14.594526,  -12.076529,    9.122148,  -13.958066,
          -0.003602,    0.019288,    0.015474,    0.000361,  100.848453)

     parameters_list_switchgrass[parameters_to_optimize] = x
     
     initial_state_switchgrass$Rhizome = 10  #mature-stand's Rhizome
     initial_state_switchgrass$Root    = 0.1

     solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
                 adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)

     output_list = list(initial_state_switchgrass,
                        parameters_list_switchgrass,
                        ss_modules_list_swithcgrass,
                        deriv_modules_list_switchgrass,solver)

   }else if(crop_type=="Energycane"){
     #Parameters for running the energycane simulations
     initial_state       <- readRDS("../calibrations/energycane//data/parameters//initial_state.rds")
     parameters_list     <- readRDS("../calibrations/energycane//data/parameters//parameters_list.rds")
     ss_modules_list     <- readRDS("../calibrations/energycane//data/parameters//ss_module_list.rds")
     deriv_modules_list  <- readRDS("../calibrations/energycane//data/parameters//deriv_modules_list.rds")
     
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
     
     solver=list(type='Gro_euler', output_step_size=1.0, adaptive_rel_error_tol=1e-4, 
                 adaptive_abs_error_tole=1e-4, adaptive_max_steps=200)
     
     output_list = list(initial_state,
                        parameters_list,
                        ss_modules_list,
                        deriv_modules_list,solver)

   }else(stop("please input a valid crop type!")) 

   return(output_list)
}
