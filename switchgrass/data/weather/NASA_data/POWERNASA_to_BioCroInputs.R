library(lubridate)
# This function coverts POWERNASA weather data to BioCro's inputs.
#NOTE: this assume we use Deepak's module incident_shortwave_par_from_cwrf, which
#takes in the solar radiation as Direct and Diffuse
#Therefore, we run the lightME function to separate the total radiation from NASA
#Also NOTE: the NASA's time stamp is NOT UTC! It's the local time. No need to convert time zone here

POWERNASA_to_BioCroInputs <- function(powernasa,latitude){
  
  source("lightME.R")
  
  solar <- powernasa$ALLSKY_SFC_PAR_TOT*4.25  #conversion from W/m2 to PPFD, See weach function. NotePOWERNASA is PAR and not total
  weatherdate <- as.Date(paste(powernasa$MO,"/",powernasa$DY,"/",powernasa$YEAR,sep=""),format=c("%m/%d/%Y"))
  doy <- yday(weatherdate)
  solardiff <- solar * lightME(lat=latitude,DOY=doy,t.d=powernasa$HR,t.sn=12,atm.P=1e5,alpha=0.6)$propIdiff
  rh=powernasa$RH2M*0.01
  solar_direct = solar - solardiff
  
  BioCroInputs <- data.frame(year=powernasa$YEAR,doy=doy,hour=powernasa$HR,temp=powernasa$T2M,rh=rh,
                             windspeed=powernasa$WS2M,precip=powernasa$PRECTOTCORR,solardiff=solardiff,
                             solar=solar_direct)
  return(BioCroInputs)
}
