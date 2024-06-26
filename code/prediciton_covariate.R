## code to generate the temperature covaraite for the prediction from the climate models ##
rm(list=ls())
args <- commandArgs(TRUE)
load("data/era5_geoinfo.RData")
load("data/precip.RData")
library(parallel)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(lubridate)
library(sf)
source("code/utility.R")
sf_use_s2(FALSE)
ls()
for (arg in args) eval(parse(text = arg))

## extract the prediciton covariates for Mississippi region ##

load("data/temperature/temperatures-mississippi.RData")
str(xyt[models=="AWI" & periods=="historical"])
str(xyt[models=="MIROC" & periods=="historical"])
str(xyt[models=="NorESM" & periods=="historical"])
str(xyt[models=="MPI" & periods=="historical"])
str(xyt[models=="CNRM" & periods=="historical"])

load("data/temperature/temperatures-danube.RData")
str(xyt[models=="AWI" & periods=="historical"])
str(xyt[models=="MIROC" & periods=="historical"])
str(xyt[models=="NorESM" & periods=="historical"])
str(xyt[models=="MPI" & periods=="historical"])
str(xyt[models=="CNRM" & periods=="historical"])

## AWI, NorESM, MIROC has the highest spatial resolution ##

load("data/temperature/temperatures-mississippi.RData")
idx.models = c("AWI","NorESM","MIROC","MPI","CNRM")
locate.idx.list <- idx.grid.list <- list()
for(idx in 1:5){
    i = which(models==idx.models[idx] & periods=="historical") ## select the model
    n.x = length(xyt[[i]]$lon); n.y = length(xyt[[i]]$lat)   
    idx.grid = as.matrix(expand.grid(x=1:n.x,y=1:n.y))
    loc.xy = data.frame(lon=xyt[[i]]$lon[idx.grid[,1]],lat=xyt[[i]]$lat[idx.grid[,2]])
    locate.idx = locate2shape(loc=loc.xy,shape=shape1)
    locate.idx = locate.idx[locate.idx[,2] %in% region.id,]
    print(lapply(region.id[-1],function(j){sum(locate.idx[,2]==j)}))
    locate.idx.list[[idx]] = locate.idx
    idx.grid.list[[idx]] = list(miss=idx.grid,danube=NULL)
}

load("data/temperature/temperatures-danube.RData")
idx.models = c("AWI","NorESM","MIROC","MPI","CNRM")
for(idx in 1:5){
    i = which(models==idx.models[idx] & periods=="historical") ## select the model
    n.x = length(xyt[[i]]$lon); n.y = length(xyt[[i]]$lat)   
    idx.grid = as.matrix(expand.grid(x=1:n.x,y=1:n.y))
    loc.xy = data.frame(lon=xyt[[i]]$lon[idx.grid[,1]],lat=xyt[[i]]$lat[idx.grid[,2]])
    locate.idx = locate2shape(loc=loc.xy,shape=shape3)
    locate.idx[,2] = region.id[1]
    print(nrow(locate.idx))
    locate.idx.list[[idx]] = rbind(locate.idx.list[[idx]],locate.idx)
    idx.grid.list[[idx]]$danube = idx.grid
}

load("data/temperature/temperatures-danube.RData")
temperature.hist <- temperature.245 <- temperature.585 <- list()
for(idx in 1:5){
    i = which(models==idx.models[idx] & periods=="historical")
    idx.temperature = idx.grid.list[[idx]]$danube[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==19,1],]
    temperature.hist[[idx]] = list()
    temperature.hist[[idx]][[1]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))

    i = which(models==idx.models[idx] & periods=="ssp245")
    idx.temperature = idx.grid.list[[idx]]$danube[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==19,1],]
    temperature.245[[idx]] = list()
    temperature.245[[idx]][[1]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))

    i = which(models==idx.models[idx] & periods=="ssp585")
    idx.temperature = idx.grid.list[[idx]]$danube[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==19,1],]
    temperature.585[[idx]] = list()
    temperature.585[[idx]][[1]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))
}

load("data/temperature/temperatures-mississippi.RData")
for(idx.region in 2:length(region.id)){
    for(idx in 1:5){
        i = which(models==idx.models[idx] & periods=="historical")
        idx.temperature = idx.grid.list[[idx]]$miss[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==region.id[idx.region],1],]
        temperature.hist[[idx]][[idx.region]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))

        i = which(models==idx.models[idx] & periods=="ssp245")
        idx.temperature = idx.grid.list[[idx]]$miss[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==region.id[idx.region],1],]
        temperature.245[[idx]][[idx.region]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))

        i = which(models==idx.models[idx] & periods=="ssp585")
        idx.temperature = idx.grid.list[[idx]]$miss[locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==region.id[idx.region],1],]
        temperature.585[[idx]][[idx.region]] = rowMeans(matrix(unlist(lapply(1:dim(temperatures[[i]])[3],function(j){temperatures[[i]][,,j][idx.temperature]})),ncol=nrow(idx.temperature),byrow=TRUE))
    }
}

temperature.hist.avg = lapply(temperature.hist,function(x){lapply(x,function(x.i){sapply(1:length(x.i),consective.mean,val=x.i,n=30)})})

temperature.245.avg = lapply(temperature.245,function(x){lapply(x,function(x.i){sapply(1:length(x.i),consective.mean,val=x.i,n=30)})})

temperature.585.avg = lapply(temperature.585,function(x){lapply(x,function(x.i){sapply(1:length(x.i),consective.mean,val=x.i,n=30)})})

date.hist = date(xyt[[1]]$time)
date.245 = date(xyt[[2]]$time)
date.585 = date(xyt[[3]]$time)

save(temperature.245.avg,temperature.hist.avg,temperature.585.avg,locate.idx.list,idx.grid.list,idx.models,temperature.hist,temperature.245,temperature.585,date.hist,date.245,date.585,file="data/temperature_pred.RData")

plot(date.hist,temperature.hist[[1]][[1]],type="l")
plot(date.245,temperature.245[[1]][[1]],type="l")
plot(date.585,temperature.585[[1]][[1]],type="l")

table(locate.idx.list[[1]][,2])
table(locate.idx.list[[2]][,2])
table(locate.idx.list[[3]][,2])
table(locate.idx.list[[4]][,2])
table(locate.idx.list[[5]][,2])

