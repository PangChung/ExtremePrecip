library(RCurl)
library(reshape2)
load("data.Rdata")
source("get_ghcn.R") #load the function (add filepath if necessary)
FTP<-"ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/"
wd <- getwd()
region <- c("Danube","Lena","Mississippi","Murry","Yangtze")
for(i in c(2,4,5)){
cwd <- paste0(wd,"/data/",region[i],"/")
dir.create(cwd)
setwd(cwd)
get_ghcn(FTP,station.df[[i]]$ID) #will put them in your working directory
source(paste0(wd,"/parse_GHCN_data.R")) #load the function (add filepath if necessary)
fnames<-list.files(path=".",pattern="*.dly")
for(fname in fnames){
	df<-parse_GHCN_data(fname)
	fname<-stringi::stri_extract(fname,regex="\\w+\\d+[^\\.]")
	#Write out the dataframe as a csv file
	write.csv(df,paste0(fname,".csv"),row.names=FALSE,na="")
}
setwd(wd)

}
