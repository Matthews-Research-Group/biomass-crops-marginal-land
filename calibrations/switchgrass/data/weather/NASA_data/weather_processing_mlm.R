source('POWERNASA_to_BioCroInputs.R')
source('lightME.R')

#latlon_unique <- read.csv('unique_latlon_and_id.csv')

latlon_unique <- data.frame(v1=99,lat=46.77,lon=-100.92)
colnames(latlon_unique)=c('id','lat','lon')

for (i in 1:nrow(latlon_unique)) {
  lat = latlon_unique$lat[i]
  lon = latlon_unique$lon[i]
  site.id = latlon_unique$id[i]
  
  NASA.weather <- read.csv(paste0("NASA_powerdata_2001_2015_site_",site.id,'.csv'),skip = 13 )

  biocroinput <- POWERNASA_to_BioCroInputs(powernasa = NASA.weather, latitude = lat)
 
  write.csv(biocroinput,file = paste0("BioCroInputs/site_", site.id ,'_lowerTransmittance.csv'), row.names = FALSE)
   
}
