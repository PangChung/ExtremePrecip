args <- commandArgs(TRUE)
source("code/utility.R")
library(ncdf4)
library(parallel)
year <- 1956

for (arg in args) eval(parse(text = arg)) #input is the year

## prepare the data and files ##
filename1 = paste0("data/2t_day_era5-land_",year,".nc")
filename2 = paste0("data/temp_EUROPE/",year,"_daily_temp.nc")
file2save = paste0("data/temperature_",year,".RData")
data1 <- nc_open(filename1)
data2 <- nc_open(filename2)
#################################################
## create era5_geoinfo.RData file if not exist ##
#################################################

if(!file.exists("data/era5_geoinfo.RData")){
library(sf)
lon <- ncvar_get(data1, "lon")
lat <- ncvar_get(data1, "lat")

shape1 <- read_sf("data/shapefiles/US_river_basins.shp")
shape2 <- read_sf("data/shapefiles/mississippi_region.shp")
shape3 <- read_sf("data/shapefiles/danube_region.shp")

sf::sf_use_s2(FALSE)

location_id = locate2shape(loc = as.data.frame( expand.grid(lon - 180, lat) ),
				shape=shape2)
ind = which(!sapply(location_id,is.integer0))
loc_df = data.frame(id = ind, group.id =  unlist(location_id[ind]) )

location_id = locate2shape(loc = as.data.frame( expand.grid(lon - 180, lat) )[ind,],
				shape=shape1)
ind = which(!sapply(location_id,is.integer0))
loc_df = data.frame(id = loc_df$id[ind], group.id =  sapply(location_id[ind],max))
loc_df$lon = as.data.frame( expand.grid(lon - 180, lat) )[loc_df$id,1]
loc_df$lat = as.data.frame( expand.grid(lon - 180, lat) )[loc_df$id,2]
## the coordinates in the NC file: (loc_df$id %% length(lon), loc_df$id %/% length(lon) + 1)

## Danube region ##
lat_danube <- ncvar_get(data2,varid = "latitude")
lon_danube <- ncvar_get(data2,varid = "longitude")
location_id = locate2shape(loc = as.data.frame( expand.grid(lon_danube, lat_danube) ),
				shape=shape3)
ind = which(!sapply(location_id,is.integer0))
loc_df_temp = data.frame(id = ind, group.id =  19,lon=expand.grid(lon_danube, lat_danube)[ind,1],
	lat=expand.grid(lon_danube, lat_danube)[ind,2])
loc_df = rbind.data.frame(loc_df_temp,loc_df)
loc_df = loc_df[!(loc_df$group.id %in% c(4,5,7)),]
save(shape1,shape2,shape3,lon,lat,lat_danube,lon_danube,loc_df,file="data/era5_geoinfo.RData")

}

# plot(shape1)
# ind = loc_df$group.id== 19
# plot(x=loc_df$lon[ind],y=loc_df$lat[ind],type="p",col=loc_df$group.id[ind],cex=0.1)
# plot(shape2)


## extract the temperature data for each region ##
load("data/era5_geoinfo.RData")
region.id = unique(loc_df$group.id)
temperature = list()
for(i in 2:length(region.id)){
	ids = loc_df$id[loc_df$group.id == region.id[i]]
	ndays = length(extract_func(ids[1],vari="2t",data=data1,res=length(lon)))
	val = matrix(data=NA,nrow=ndays,ncol=length(ids))
	for(j in 1:length(ids)){
		val[,j] = extract_func(ids[j],vari="2t",data=data1,res=length(lon))
	}
	temperature[[i]] = val
}

## extract the temperature data for Danube region ##
ids = loc_df$id[loc_df$group.id == 19]
ndays = length(extract_func(ids[1],vari="t2m",data=data2,res=length(lon_danube)))
val = matrix(data=NA,nrow=ndays,ncol=length(ids))
for(j in 1:length(ids)){
	val[,j] = extract_func(ids[j],vari="t2m",data=data2,res=length(lon_danube))
}
temperature[[1]] = val
save(temperature,file=file2save)

# for(i in 1950:2020){
# filename = paste0("data/temperature_",i,".RData")
# data2 <- nc_open(paste0("../Temperatures/temp_EUROPE/",i,"_daily_temp.nc"))
# load(filename)
# ids = loc_df$id[loc_df$group.id == 19]
# ndays = length(extract_func(ids[1],vari="t2m",data=data2,res=length(lon_danube)))
# val = matrix(data=NA,nrow=ndays,ncol=length(ids))
# for(j in 1:length(ids)){
# 	val[,j] = extract_func(ids[j],vari="t2m",data=data2,res=length(lon_danube))
# }
# temperature[[1]] = val
# save(temperature,file=filename)
# }
