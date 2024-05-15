## code to generate the temperature covaraite after extracting the temperature data using code extrect_temperature.R ##
library(lubridate)
library(dplyr)
source("code/utility.R")
load("data/era5_geoinfo.RData")
files <- list.files(path="data/",pattern="temperature_\\d+",full.names=TRUE)
date <- seq(as.Date("1950-01-01"), as.Date("2020-12-31"),1)

temperature <- list(Danube=NULL,Miss1=NULL,Miss2=NULL,Miss3=NULL,
					Miss4=NULL,Miss5=NULL,Miss6=NULL,Miss7=NULL)
for(f in files){
	load(f,e<-new.env())
	temperature <- lapply(1:8,FUN= function(id) c(temperature[[id]],apply(e$temperature[[id]],1,mean) - 273.15 ) )
	print(f)
}

temperature.covariate = temperature
temperature.covariate <- lapply(temperature,function(x) sapply(1:length(x),consective.mean,val=x,n=30))
date.df = data.frame(date=date,year=getYear(date),season=getSeason(date)) 


save(temperature,temperature.covariate,date.df,loc_df,file="data/temperature.RData")



