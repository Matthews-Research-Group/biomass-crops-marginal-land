us <- readOGR(dsn = "./",layer = "us_conus")
cwrf_proj4="+proj=lcc +lat_1=30 +lat_2=60 +lat_0=37.49999 +lon_0=-95.5 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
us_lcc=spTransform(us,CRS(cwrf_proj4))
dsn='/Users/yufeng/Desktop/model_coupling_UMD/USshp/us_lcc_YH/'
writeOGR(us_lcc,dsn=dsn,layer = "us_lcc",driver="ESRI Shapefile")
