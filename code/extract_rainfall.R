## code to extract the precipitation data from raw files for each subregion ##
library(sf)
load("data/data.Rdata")
source("code/utility.R")

shape1 <- read_sf("data/shapefiles/US_river_basins.shp")
shape3 <- read_sf("data/shapefiles/danube_region.shp")
station.US <- station.df[[3]]
location_id <- locate2shape(data.frame(x=station.US$Y,y=station.US$X),shape1)
ind = which(!sapply(location_id,is.integer0))
station.US = station.US[ind,]
station.US$group.id = unlist(location_id[ind])

station.danube = station.df[[1]]
station.danube$group.id = 19

station = rbind(station.danube,station.US)
region.name <- shape1$names[unique(station$group.id)]
region.id <- unique(station$group.id)
region.name[1] <- "Danube"

### extract the precipitation data from raw files ###
precip <- list(Danube=NULL,Miss1=NULL,Miss2=NULL,Miss3=NULL,
					Miss4=NULL,Miss5=NULL,Miss6=NULL,Miss7=NULL)
START.date = as.Date("1965-01-01"); END.date = as.Date("2020-12-31")
n = as.integer(END.date - START.date + 1)
for(i in 1:length(region.id)){
	station.sub = subset(station, group.id==region.id[i])
	for(j in 1:nrow(station.sub)){
	file = ifelse(region.name[i]=="Danube",paste0("data/Danube/",station.sub$ID[j],".csv"),
			paste0("data/Mississippi/",station.sub$ID[j],".csv"))
	data = read.csv2(file,header=TRUE,sep = ",")
	date = as.integer(as.Date(data$date) - START.date + 1)
	print(range(date))
	ind = date >= 1 & date <= n 
	precip[[i]][[j]] = rep(NA, n)
	precip[[i]][[j]][date[ind]] = data$value[ind]
	}
	print(i)
}

save(precip,START.date,END.date,station,region.id,region.name,file="data/precip.RData")

## compare the dataset ##
load("~/Documents/R/rainfall_project/data/data.Rdata",e<-new.env())
date.old = e$date.ts
date = seq(START.date,END.date,1)
date.inter = which(date %in% date.old) 
station.sub = subset(station,group.id == 19)
ind = 9
dat = precip[[1]][[ind]][date.inter]/10 ## the two datasets have different units
dat.old = e$precip.ts.df[[1]][,ind+3]
dat[dat>1000 | dat < 0] = NA
idx.11 = !is.na(dat) & !is.na(dat.old)
idx.10 = !is.na(dat) & is.na(dat.old)
idx.01 = is.na(dat) & !is.na(dat.old)
idx.00 = is.na(dat) & is.na(dat.old)
sum(idx.11)
sum(idx.10)
sum(idx.01)
sum(idx.00)


