source('POWERNASA_to_BioCroInputs.R')
source('lightME.R')

latlon_unique <- read.csv('unique_latlon_and_id.csv')

for (i in 1:nrow(latlon_unique)) {

  NASA.weather <- read.csv(paste0("NASA_powerdata_2001_2015_site_",i,'.csv'),skip = 13 )
  
  lat = latlon_unique$lat[i]
  lon = latlon_unique$lon[i]
  site.id = latlon_unique$id[i]
  
  biocroinput <- POWERNASA_to_BioCroInputs(powernasa = NASA.weather, latitude = lat)
 
  write.csv(biocroinput,file = paste0("BioCroInputs/site_", site.id ,'_lowerTransmittance.csv'), row.names = FALSE)
   
}
