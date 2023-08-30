source('POWERNASA_to_BioCroInputs.R')
weather_example = read.csv('weather_calibrationsite.csv') #this provides the period we need to truncate NASA's
NASA = read.csv('POWER_Weather_Hourly_20100101_20201231_Labelle.csv')
weatherdate <- as.Date(paste(NASA$MO,"/",NASA$DY,"/",NASA$YEAR,sep=""),format=c("%m/%d/%Y"))
doy <- yday(weatherdate)
ind_begin = which(NASA$YEAR==weather_example$year[1] & doy == weather_example$doy[1])[1]
ind_end = which(NASA$YEAR==tail(weather_example$year,1) & doy == tail(weather_example$doy,1))[24]
NASA_sub = NASA[ind_begin:ind_end,]

labelle_latitude = 26.76
  
biocro_weather = POWERNASA_to_BioCroInputs(NASA_sub,labelle_latitude)
write.csv(biocro_weather,'weather_NASA_calibration_LaBelle.csv')
