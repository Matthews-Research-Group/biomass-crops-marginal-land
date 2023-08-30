output_for_different_harvest_dates <- function(result,day,month,year,measure_day_known){

  
   #what is day of year on march 15
  march15_doy <-  yday(as.Date(paste(sprintf("%04d",year),"-",
                        sprintf("%02d",3),"-",
                        sprintf("%02d",15),
                        sep="")))
  
  # doy for apri 15
  april15_doy <-  yday(as.Date(paste(sprintf("%04d",year),"-",
                                     sprintf("%02d",4),"-",
                                     sprintf("%02d",15),
                                     sep="")))
  #what is actual day of year on the harves day.
  #if harvest day is not known assume that march 15 is the harvest day
  
  if(measure_day_known==1){
   actual_harvest_doy <-  yday(as.Date(paste(sprintf("%04d",year),"-",
                                       sprintf("%02d",month),"-",
                                       sprintf("%02d",day),
                                       sep="")))
  } else {
    actual_harvest_doy = march15_doy
  }
  
  #What is doy on the last day durign growing period or on the day of peak biomass
  peak_biomass_doy <- result$doy[dim(result)[1]-1]
  

# if measurement date is after april 15 and less than peak biomass then harvesting occured before reaching peak 
#  & there is no need to remove extra winter senescencd biomass
  if((actual_harvest_doy >= april15_doy) & (actual_harvest_doy <=  peak_biomass_doy)){
  actual_harvested_biomass <- result[(result$hour==0) & (result$doy ==actual_harvest_doy),]$Stem
  } else{
# otherwise we substract 0.07 tons/ha per day from peak biomass   
    #if harvest month is either in jan/feb/march then
    # loss days are calculated  by (365 - peak doy) + doy on actual harvest day
    if(actual_harvest_doy < april15_doy){
      number_of_loss_day =  (365 - peak_biomass_doy) + actual_harvest_doy
    } else{
      number_of_loss_day = actual_harvest_doy - peak_biomass_doy
    }
    actual_harvested_biomass <- result$Stem[dim(result)[1]] - 0.07*number_of_loss_day
  }
 peak_harvested_biomass <- result$Stem[dim(result)[1]]
 march15_harvested_biomass <- result$Stem[dim(result)[1]] - 0.07* ((365 - peak_biomass_doy) + march15_doy)
 
 output <- list(actual_harvested_biomass = actual_harvested_biomass, peak_harvested_biomass= peak_harvested_biomass,
                march15_harvested_biomass=march15_harvested_biomass )
 return(output)
}