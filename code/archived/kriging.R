## kriging for the temperature data ##
rm(list=ls())
args<-commandArgs(TRUE) 
idx.region=1;idx.model=1
for (arg in args) eval(parse(text = arg))
library(parallel)
library(fields)
library(mgcv)
library(mvnfast)
library(lubridate)
library(sf)
library(ggplot2)
library(RColorBrewer)

###  load the data ###
load("data/precip.RData")
load("data/temperatures-danube.RData",danube <-new.env())
load("data/temperatures-mississippi.RData",miss <- new.env())
load("data/map_danube.Rdata",danube_map <- new.env())
load("data/map_mississippi.Rdata",miss_map <- new.env())
source("code/utility.R")

shape1 <- read_sf("data/shapefiles/US_river_basins.shp")
shape2 <- read_sf("data/shapefiles/danube_region.shp")

lat = danube$xyt[[1]]$lat 
lon = danube$xyt[[1]]$lon
lat.lon <- expand.grid(lat,lon)
lat.lon <- data.frame(x=lat.lon[,2],y=lat.lon[,1])

id.lat.lon <- locate2shape(lat.lon,shape1)
coord.elev.df = as.data.frame(danube_map$coord.elev.df)
coord.elev.df = coord.elev.df[!is.na(coord.elev.df$danubez7),]

p1 <- ggplot() + geom_sf(data=shape2,mapping=aes(),fill="red")
p1 <- p1 + geom_raster(data=coord.elev.df,aes(x=lon,y=lat,fill=danubez7)) + scale_fill_gradient(low="white",high="black") 
p1 <- p1 + geom_point(data=lat.lon,aes(x,y)) + labs(x="Longitude",y="Latitude",fill="Elevation (m)")
p1

### prepare the data
coord.elev.df = as.list(as.data.frame(coord.elev.df))
N.lon = length(xyt[[i]]$lon)
N.lat = length(xyt[[i]]$lat)
station.loc <-expand.grid(xyt[[i]]$lon,xyt[[i]]$lat)
names(station.loc) <- c("lon","lat");station.loc <- as.data.frame(station.loc)
alt <- coord.elev.df[[1]][sapply(1:nrow(station.loc),function(id){which.min((coord.elev.df$lon-station.loc[id,1])^2 + (coord.elev.df$lat-station.loc[id,2])^2)})]
xyt[[i]]$alt = alt

### prepare the data for the GAM model ###
ind = !is.na(xyt[[i]]$alt)
D = sum(ind);N = length(xyt[[i]]$time)
lon = station.loc[ind,1]
lat = station.loc[ind,2]
alt = xyt[[i]]$alt[ind]
dayy = yday(xyt[[i]]$time)
year = xyt[[i]]$year; year = year - min(year) + 1
y = c(apply(which(matrix(ind,N.lon,N.lat),arr.ind=T),1,function(id){temperatures[[i]][id[1],id[2],]}))
data.df = data.frame(y=y,lon=rep(lon,each=N),lat = rep(lat,each=N),alt=rep(alt,each=N),dayy= rep(dayy,D),year=rep(year,D))
## fit the GAM model to extract the mean ##
fit <- gam(formula=y~ti(lat,lon)+ti(alt)+ti(dayy)+ti(year),data=data.df,family=gaussian(link="identity"),fit=TRUE)
data.fit <- matrix(residuals(fit),nrow=N,ncol=D,byrow=FALSE)

## compute the distance 
coords <- cbind(lon,lat)
dist.mat <- rdist.earth(coords,miles=FALSE)
pairs <- t(combn(D,2))
dist.vec <- dist.mat[pairs]
## fit a Gaussian dependence model for the residuals 
D.ind = dist.vec < 500
pairs <- pairs[D.ind,]
dist.vec <- dist.vec[D.ind]
fit.composite <- function(par){
    par = exp(par)
    cov.vec <- par[1]*exp(-(dist.vec/par[2])^par[3])
    fun <- function(ind){
      A0 <- apply(data.fit[,pairs[ind,]],1,function(x){sum(is.na(x))==0})
      val <- rep(NA,N)
      if(sum(A0)!=0){
        val[A0] <- dmvn(data.fit[A0,pairs[ind,]],mu=c(0,0),sigma=matrix(c(par[1],cov.vec[ind],cov.vec[ind],par[1]),2,2),log=TRUE,ncores=1)
      }
      return(val)
    }
    val <- -sum(unlist(mclapply(1:nrow(pairs),FUN = fun,mc.cores = 20,mc.preschedule = TRUE)),na.rm=T)
}
  
  fit.cov <- optim(log(c(26,440,0.42)),fit.composite,method = "Nelder-Mead",control = list(maxit=200,trace=TRUE))
  para.mat <- exp(fit.cov$par)
  
  nx.coord = length(unique(coord.elev.df$lon))
  ny.coord = length(coord.elev.df$lon)/nx.coord
  
  loc2pred = expand.grid(ceiling(100*c(0:floor(nx.coord/100))+1),ceiling(100*c(0:floor(ny.coord/100))+1) )
  loc2pred = lapply(coord.elev.df,function(x){x[(loc2pred[,2] - 1) * nx.coord + loc2pred[,1]]})
  loc2pred = lapply(loc2pred,function(x){x[!is.na(loc2pred[[1]])]})
  
  data2pred = data.frame(lat = rep(loc2pred$lat,each=N),lon=rep(loc2pred$lon,each=N),alt=rep(loc2pred[[1]],each=N),dayy=rep(dayy,length(loc2pred[[1]])),year=rep(year,length(loc2pred[[1]])))
  n.cores=1000
  L = nrow(data2pred) %/% n.cores
  mc.seq = c(seq(0,nrow(data2pred),L),nrow(data2pred))
  data.pred.list = sapply(1:(length(mc.seq)-1),function(ind){data2pred[(mc.seq[ind]+1):mc.seq[ind+1],]},simplify = FALSE)
  
  message("start fit prediction")
  data.pred = mclapply(data.pred.list,function(dat){as.numeric(predict(fit,dat,type="response"))},mc.cores=20,mc.preschedule = FALSE)
  data.pred = matrix(unlist(data.pred),nrow=N,ncol=length(loc2pred[[1]]),byrow=FALSE)

  N.pred = length(loc2pred[[1]])
  param=para.mat
  Sigma22 <- dist.mat;Sigma22 <- param[1]*exp(-(Sigma22/param[2])^param[3])
  Sigma12 <- rdist.earth(x1=cbind(loc2pred$lon,loc2pred$lat),x2=coords,miles=FALSE); Sigma12 <- param[1]*exp(-(Sigma12/param[2])^param[3])
  ind.list <- sapply(1:nrow(data.fit),function(id){ind=!is.na(data.fit[id,]);list(ind=ind,data=data.fit[id,ind])},simplify = F)
  
  Kriging <- function(ind){
    fun<-function(id){
      dat = id$data; id.ind<-id$ind
      if(length(dat) < 2){return(rep(0,N.pred))}      
      if(length(dat)>100){id.temp = order(apply(Sigma12[,id.ind],2,sum))[1:100];dat=dat[id.temp]}
      else{id.temp <- 1:sum(id.ind)}
      val = Sigma12[,id.ind][,id.temp] %*% solve(Sigma22[id.ind,id.ind][id.temp,id.temp]) %*% dat
      return(c(val))
    }
    kri <- mcmapply(fun,id=ind.list[ind],SIMPLIFY = F,mc.cores=28,mc.preschedule = TRUE)
    kri <- matrix(unlist(kri),ncol=N.pred,byrow=T)
    return(kri)
  }
  
  n.cores=100
  L = length(ind.list) %/% n.cores
  mc.seq = c(seq(0,length(ind.list),L),length(ind.list))
  id.list <- sapply(1:(length(mc.seq)-1),function(id){(mc.seq[id]+1):mc.seq[id+1]})
  message("start kriging")
  pr <- lapply(id.list,Kriging)
  pr <- do.call(rbind,pr)
  
  pred.tep.mean.list <- data.pred
  pred.tep.list <- pr + data.pred
  save.image(file=paste0("data/data_temp_",i,"_",regions[region.ind],".Rdata"),compress = "xz")

#### extract and summerize the results ####
  rm(list=ls())
  data.files <- list.files(path = "data/",pattern = "data_temp_.+mississippi.Rdata",full.names = T)
  pred.tep.list <- pred.tep.mean.list <- xyt <- list() 
  for(j in 1:length(data.files)){
    load(data.files[j],e<-new.env())
    pred.tep.list[[e$i]] <- e$pred.tep.list
    pred.tep.mean.list[[e$i]] <- e$pred.tep.mean.list
    xyt[[e$i]] <- e$xyt[[e$i]]
  }
  temperatures <- e$temperatures
  coord.elev.df <- e$coord.elev.df
  loc2pred <- e$loc2pred
  periods <- e$periods
  models <- e$models
  save(pred.tep.list,pred.tep.mean.list,xyt,temperatures,coord.elev.df,loc2pred,periods,models,file="data/data_temp_mississippi.Rdata")
  
  rm(list=ls())
  data.files <- list.files(path = "data/",pattern = "data_temp_.+danube.Rdata",full.names = T)
  pred.tep.list <- pred.tep.mean.list <- xyt <- list() 
  for(j in 1:length(data.files)){
    load(data.files[j],e<-new.env())
    pred.tep.list[[e$i]] <- e$pred.tep.list
    pred.tep.mean.list[[e$i]] <- e$pred.tep.mean.list
    xyt[[e$i]] <- e$xyt[[e$i]]
  }
  temperatures <- e$temperatures
  coord.elev.df <- e$coord.elev.df
  loc2pred <- e$loc2pred
  periods <- e$periods
  models <- e$models
  save(pred.tep.list,pred.tep.mean.list,xyt,temperatures,coord.elev.df,loc2pred,periods,models,file="data/data_temp_danube.Rdata")
  
  ###########################################
  ### realign the temperature predictions ###
  ###########################################
  rm(list=ls())
  load("data/data_pot.Rdata",e<-new.env())
  load("data/data_temp_mississippi.Rdata")
  load("data/temperatures-mississippi.RData")
  region.ind=3
  library(lubridate)
  getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
    return(s.id)
  }
  data.pred.list <- list()
  for(i in 1:15){
    if(!is.null(pred.tep.list[[i]]) && !date(rev(xyt[[i]]$time)[1]) > date("2100-12-31") ){
    if(periods[i]=="historical"){
      id.pred <- date(xyt[[i]]$time) %in% date(e$date.ts)
      id.obs <- as.Date(e$date.ts) %in% as.Date(xyt[[i]]$time)
    }else{
      id.pred <- year(xyt[[i]]$time) %in% 2015:2020
      id.obs <- year(e$date.ts) %in% 2010:2015
    }
    data.pred <- apply(pred.tep.list[[i]],1,mean)
    data.obs <- apply(e$pred.tep.list[[region.ind]],1,mean)
    value.diff1 <- aggregate(data.obs[id.obs],by=list(season=unlist(parallel::mclapply(e$date.ts[id.obs],getSeason,mc.cores=20))),mean)
    value.diff2 <- aggregate(data.pred[id.pred],by=list(season=unlist(parallel::mclapply(xyt[[i]]$time[id.pred],getSeason,mc.cores=20))),mean)
    print(value.diff1$season)
    print(value.diff2$season)
    value.diff  = value.diff1$x - value.diff2$x
    xyt[[i]]$season <- as.factor(unlist(parallel::mclapply(xyt[[i]]$time,getSeason,mc.cores=20)))
    data.pred.list[[i]] <- data.pred + value.diff[xyt[[i]]$season]
    }
  }
  date.ts  <- e$date.ts
  pred.obs.list <- e$pred.tep.list
  save(data.pred.list,pred.obs.list,xyt,date.ts,periods,models,file="data/aligned_temperature_mississippi.Rdata")

  rm(list=ls())
  load("data/data_pot.Rdata",e<-new.env())
  load("data/data_temp_danube.Rdata")
  region.ind=1
  library(lubridate)
  getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
    return(s.id)
  }
  data.pred.list <- list()
  for(i in 1:15){
    if(!is.null(pred.tep.list[[i]]) && !date(rev(xyt[[i]]$time)[1]) > date("2100-12-31") ){
      if(periods[i]=="historical"){
        id.pred <- date(xyt[[i]]$time) %in% date(e$date.ts)
        id.obs <- as.Date(e$date.ts) %in% as.Date(xyt[[i]]$time)
      }else{
        id.pred <- year(xyt[[i]]$time) %in% 2015:2020
        id.obs <- year(e$date.ts) %in% 2010:2015
      }
      data.pred <- apply(pred.tep.list[[i]],1,mean)
      sum(is.na(data.pred))
      data.obs <- apply(e$pred.tep.list[[region.ind]],1,mean)
      value.diff1 <- aggregate(data.obs[id.obs],by=list(season=unlist(parallel::mclapply(e$date.ts[id.obs],getSeason,mc.cores=20))),mean)
      value.diff2 <- aggregate(data.pred[id.pred],by=list(season=unlist(parallel::mclapply(xyt[[i]]$time[id.pred],getSeason,mc.cores=20))),mean)
      print(value.diff1$season)
      print(value.diff2$season)
      value.diff  = value.diff1$x - value.diff2$x
      xyt[[i]]$season <- as.factor(unlist(parallel::mclapply(xyt[[i]]$time,getSeason,mc.cores=20)))
      data.pred.list[[i]] <- data.pred + value.diff[xyt[[i]]$season]
    }
  }
  date.ts  <- e$date.ts
  pred.obs.list <- e$pred.tep.list
  save(data.pred.list,pred.obs.list,xyt,date.ts,periods,models,file="data/aligned_temperature_danube.Rdata")


  load("data/aligned_temperature_danube.Rdata")
  model.selected <- c("AWI","MIROC","NorESM")
  inds <- which(models %in% model.selected)
  for(i in inds){
    print(paste(models[i],":",diff(xyt[[i]]$lat)[1],diff(xyt[[i]]$lon)[1]))
  }
  load("data/aligned_temperature_mississippi.Rdata")
  model.selected <- c("AWI","MIROC","NorESM")
  inds <- which(models %in% model.selected)
  for(i in inds){
    print(paste(models[i],":",diff(xyt[[i]]$lat)[1],diff(xyt[[i]]$lon)[1]))
  }
  

  rm(list=ls())
	getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
    return(s.id)
  }

  load("data/aligned_temperature_danube.Rdata")
  for(i in 1:15){
	xyt[[i]]$season <- as.factor(unlist(parallel::mclapply(xyt[[i]]$time,getSeason,mc.cores=4)))
  }
  save(data.pred.list,pred.obs.list,xyt,date.ts,periods,models,file="data/aligned_temperature_danube.Rdata")

  load("data/aligned_temperature_mississippi.Rdata")
  for(i in 1:15){
	xyt[[i]]$season <- as.factor(unlist(parallel::mclapply(xyt[[i]]$time,getSeason,mc.cores=4)))
  }
  save(data.pred.list,pred.obs.list,xyt,date.ts,periods,models,file="data/aligned_temperature_mississippi.Rdata")

  #############################################################################
  #### calculate the dependence range with the moving averages ###

  solve.h.BR <- function(par,temp,logval=FALSE){
    alpha = 2/(1+exp(-par[1]))
    lambda = exp(par[2] + par[3]*temp)
    if(logval){
      h = log(lambda) + (1/alpha) * (log(2*qnorm(1.95/2)^2))
    }else{
      h=lambda * (2*qnorm(1.95/2)^2)^(1/alpha)
    }
    return(h)
  }
  
  emp.h.BR <- function(obs,u){
    ind.obs1 <- obs[,1] > u
    ind.obs2 <- obs[,2] > u
    val <- sum(ind.obs1 & ind.obs2)/sum(ind.obs1 | ind.obs2)
    return(val)
  }
  
  library(lubridate)
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","lena","mississippi","murray","yangtze")
  lab.regions <- c("Danube","Lena","Mississippi","Murray","Yangtze")
  
  consective.mean <- function(pos,vals,n=30){
    if(pos <= n){
      return(vals[pos])
    }else{
      return(mean(vals[(pos-n+1):pos]))
    }
  }
  rFun <- function(x,u,xi=est.shape.gpd){
    val = mean((x/u)^xi)^{1/xi}
    return(val)
  }
  #select time point
  ind.fit.list <- list();norm.constant.list <- list()
  norm.ind=1
  for(region.ind in c(1,3)){
    ind.fit <- matrix(NA,ncol=4,nrow=length(season.list))
    norm.constant <- matrix(NA,ncol=4,nrow=2)
    load(paste0("data/data_pot_10_",regions[region.ind],".Rdata"))
    if(norm.ind==1){
      est.shape.gpd <- fitted(results.gpd)[1,2]
    }else{
      est.shape.gpd <- 1}
    for(season.ind in 1:4){
      n.skip = 30
      ind.season = sapply(season.list, function(x){x[1]==season[season.ind]})
      obs = split(U,row(U))
      obs = subset(obs,ind.season)
      no.obs = sapply(obs,function(x){sum(!is.na(x))})
      obs[no.obs>0] = lapply(obs[no.obs>0],evd::qgpd,shape=1,loc=1,scale = 1)
      tep = temperature.st.moving;
      tep = tep[ind.season];reg.t = tep
      r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],est.shape.gpd)}else{NA}})))
      thres = quantile(r.obs[no.obs > 5],0.8,na.rm=TRUE)
      ind.exc = no.obs > 5 & r.obs > thres 
      reg.t = reg.t[ind.exc]
      norm.constant[,season.ind] = c(mean(reg.t),sd(reg.t))
      
      ind.season[ind.season] <- ind.exc
      ind.fit[,season.ind] <- ind.season
      
      print(paste(season[season.ind],regions[region.ind],sum(ind.season)))
    }
    ind.fit.list[[region.ind]] = ind.fit
    norm.constant.list[[region.ind]] = norm.constant
  }
  
  getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
    return(s.id)
  }
  model.selected <- c("AWI","MIROC","NorESM")
  data.mean.tep.245.day <- list()
  data.mean.tep.585.day <- list()
  data.mean.obs.hist.day <- list()
  data.mean.tep.hist.day <- list()
  
  data.mean.tep.245 <- list()
  data.mean.tep.585 <- list()
  data.mean.obs.hist <- list()
  data.mean.tep.hist <- list()
  
  load("data/result_pot_10_moving.Rdata",e.fit<-new.env())
  for(region.ind in c(1,3)){
    load(paste0("data/aligned_temperature_",regions[region.ind],".Rdata"))
    load(paste0("data/data_pot_10_",regions[region.ind],".Rdata"))
    season.list <- as.factor(sapply(date.ts,getSeason))
    print(table(periods))
    data.mean.tep.245[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.tep.585[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.obs.hist[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.tep.hist[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.tep.245.day[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.tep.585.day[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.obs.hist.day[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    data.mean.tep.hist.day[[which(c(1,3)==region.ind)]] <- list(list(),list(),list())
    for(model.name in model.selected){
      id.hist <- which(periods == "historical" & models == model.name)
      id.245 <- which(periods ==  "ssp245" & models == model.name)
      id.585 <- which(periods ==  "ssp585" & models == model.name)
      mean.tep.obs.hist.list <- apply(pred.obs.list[[region.ind]],1,mean)
      mean.tep.obs.hist.list <- sapply(1:length(mean.tep.obs.hist.list),consective.mean,vals=mean.tep.obs.hist.list,n=n.skip,simplify = TRUE) 
      ind.temp <- season.list;ind.temp <- ordered(ind.temp,levels=season)
      data.mean.obs.hist.day[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] = data.frame(date=date.ts,tep=mean.tep.obs.hist.list)
      mean.tep.obs.hist.list <- (mean.tep.obs.hist.list - norm.constant.list[[region.ind]][1,ind.temp])/norm.constant.list[[region.ind]][2,ind.temp]  
      mean.tep.obs.hist.list <- aggregate(mean.tep.obs.hist.list,by=list(year=year(date.ts),season=season.list),mean)
      if(is.null(data.pred.list[[id.245]])){
        print("haha")
      }else{
        print(paste0("working",":",model.name))
        mean.tep.hist.list <- sapply(1:length(data.pred.list[[id.hist]]),consective.mean,vals=data.pred.list[[id.hist]],n=n.skip,simplify = TRUE) 
        data.mean.tep.hist.day[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] <- data.frame(date=xyt[[id.hist]]$time,tep=mean.tep.hist.list)[!(xyt[[id.hist]]$year %in% c(1850:1964)),] 
        ind.temp <- xyt[[id.hist]]$season;ind.temp <- ordered(ind.temp,levels=season)
        
        mean.tep.hist.list <- (mean.tep.hist.list - norm.constant.list[[region.ind]][1,ind.temp])/norm.constant.list[[region.ind]][2,ind.temp]  
        mean.tep.hist.list <- aggregate(mean.tep.hist.list,by=list(year=year(xyt[[id.hist]]$time),season=xyt[[id.hist]]$season),mean)
        
        mean.tep.245.list <- sapply(1:length(data.pred.list[[id.245]]),consective.mean,vals=data.pred.list[[id.245]],n=n.skip,simplify = TRUE) 
        data.mean.tep.245.day[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] <- data.frame(date=xyt[[id.245]]$time,tep=mean.tep.245.list)
        ind.temp <- xyt[[id.245]]$season;ind.temp <- ordered(ind.temp,levels=season)
        mean.tep.245.list <- (mean.tep.245.list - norm.constant.list[[region.ind]][1,ind.temp])/norm.constant.list[[region.ind]][2,ind.temp]  
        mean.tep.245.list <- aggregate(mean.tep.245.list,by=list(year=year(xyt[[id.245]]$time),season=xyt[[id.245]]$season),mean)
        
        mean.tep.585.list <- sapply(1:length(data.pred.list[[id.585]]),consective.mean,vals=data.pred.list[[id.585]],n=n.skip,simplify = TRUE) 
        data.mean.tep.585.day[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] <- data.frame(date=xyt[[id.585]]$time,tep=mean.tep.585.list)
        ind.temp <- xyt[[id.585]]$season;ind.temp <- ordered(ind.temp,levels=season)
        mean.tep.585.list <- (mean.tep.585.list - norm.constant.list[[region.ind]][1,ind.temp])/norm.constant.list[[region.ind]][2,ind.temp]  
        mean.tep.585.list <- aggregate(mean.tep.585.list,by=list(year=year(xyt[[id.585]]$time),season=xyt[[id.585]]$season),mean)
        
        data.mean.obs.hist[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] = mean.tep.obs.hist.list
        data.mean.tep.hist[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] = mean.tep.hist.list[!(mean.tep.hist.list$year %in% c(1850:1964)),] 
        data.mean.tep.245[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] <- mean.tep.245.list
        data.mean.tep.585[[which(c(1,3)==region.ind)]][[which(model.selected==model.name)]] <- mean.tep.585.list
      
    }
  }}
  
  names(data.mean.tep.245) <- names(data.mean.obs.hist) <- names(data.mean.tep.585) <- names(data.mean.tep.hist) <- regions[c(1,3)]
  for(r in 1:2){
      names(data.mean.tep.245[[r]]) <- names(data.mean.obs.hist[[r]]) <- names(data.mean.tep.585[[r]]) <- names(data.mean.tep.hist[[r]]) <- model.selected
  }
  
  ## compute the averages ##
  load("data/temperature_covaraites.Rdata")
  for(i in 1:2){
    data.mean.tep.245[[i]]$AVG <- data.mean.tep.245[[i]][[1]]
    data.mean.tep.245[[i]]$AVG$x <- 1/3*(data.mean.tep.245[[i]][[1]]$x + data.mean.tep.245[[i]][[2]]$x + data.mean.tep.245[[i]][[3]]$x)
    
    data.mean.obs.hist[[i]]$AVG <- data.mean.obs.hist[[i]][[1]]
    data.mean.obs.hist[[i]]$AVG$x <- 1/3*(data.mean.obs.hist[[i]][[1]]$x + data.mean.obs.hist[[i]][[2]]$x + data.mean.obs.hist[[i]][[3]]$x)
    
    data.mean.tep.585[[i]]$AVG <- data.mean.tep.585[[i]][[1]]
    data.mean.tep.585[[i]]$AVG$x <- 1/3*(data.mean.tep.585[[i]][[1]]$x + data.mean.tep.585[[i]][[2]]$x + data.mean.tep.585[[i]][[3]]$x)
    
    data.mean.tep.hist[[i]]$AVG <- data.mean.tep.hist[[i]][[1]]
    data.mean.tep.hist[[i]]$AVG$x <- 1/3*(data.mean.tep.hist[[i]][[1]]$x + data.mean.tep.hist[[i]][[2]]$x + data.mean.tep.hist[[i]][[3]]$x)
    
    data = data.mean.tep.245.day[[i]][[1]];data$tep = NA
    data$tep[1:length(data.mean.tep.245.day[[i]][[3]]$date)] = data.mean.tep.245.day[[i]][[3]]$tep
    data.mean.tep.245.day[[i]][[3]] <- data 
    
    data = data.mean.tep.585.day[[i]][[1]];data$tep = NA
    data$tep[1:length(data.mean.tep.585.day[[i]][[3]]$date)] = data.mean.tep.585.day[[i]][[3]]$tep
    data.mean.tep.585.day[[i]][[3]] <- data 
    
    data = data.mean.tep.hist.day[[i]][[1]];data$tep = NA
    data$tep[1:length(data.mean.tep.hist.day[[i]][[3]]$date)] = data.mean.tep.hist.day[[i]][[3]]$tep
    data.mean.tep.hist.day[[i]][[3]] <- data 
    
    data.mean.tep.hist.day[[i]][[4]] <- data.mean.tep.hist.day[[i]][[1]]
    data.mean.tep.hist.day[[i]][[4]]$tep <- 1/3*(data.mean.tep.hist.day[[i]][[1]]$tep + data.mean.tep.hist.day[[i]][[2]]$tep + data.mean.tep.hist.day[[i]][[3]]$tep)
    
    data.mean.obs.hist.day[[i]][[4]] <- data.mean.obs.hist.day[[i]][[1]]
    data.mean.obs.hist.day[[i]][[4]]$tep <- 1/3*(data.mean.obs.hist.day[[i]][[1]]$tep + data.mean.obs.hist.day[[i]][[2]]$tep + data.mean.obs.hist.day[[i]][[3]]$tep)
    
    data.mean.tep.245.day[[i]][[4]] <- data.mean.tep.245.day[[i]][[1]]
    data.mean.tep.245.day[[i]][[4]]$tep <- 1/3*(data.mean.tep.245.day[[i]][[1]]$tep + data.mean.tep.245.day[[i]][[2]]$tep + data.mean.tep.245.day[[i]][[3]]$tep)
    
    data.mean.tep.585.day[[i]][[4]] <- data.mean.tep.585.day[[i]][[1]]
    data.mean.tep.585.day[[i]][[4]]$tep <- 1/3*(data.mean.tep.585.day[[i]][[1]]$tep + data.mean.tep.585.day[[i]][[2]]$tep + data.mean.tep.585.day[[i]][[3]]$tep)
    
  }
  
  save(ind.fit.list,data.mean.tep.245,data.mean.tep.585,data.mean.tep.hist,data.mean.obs.hist,
       data.mean.tep.245.day,data.mean.tep.585.day,data.mean.tep.hist.day,data.mean.obs.hist.day,file="data1/temperature_covaraites.Rdata")
  
  ## use ggplot2 ##
  #rm(list=ls())
  load("data/result_pot_10_moving.Rdata")
  library(ggplot2)
  library(cowplot)
  model.selected <- c("AWI","MIROC","NorESM")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","lena","mississippi","murray","yangtze")
  lab.regions <- c("Danube","Lena","Mississippi","Murray","Yangtze")
  
solve.h.BR <- function(par,temp,logval=FALSE){
    alpha = par[1]
    lambda = par[2] + par[3]*temp
    if(logval){
      h = lambda + (1/alpha) * (log(2*qnorm(1.95/2)^2))
    }else{
      h=exp(lambda)*(2*qnorm(1.95/2)^2)^(1/alpha)
    }
    return(h)
  }
  
emp.h.BR <- function(obs,u){
    ind.obs1 <- obs[,1] > u
    ind.obs2 <- obs[,2] > u
    val <- sum(ind.obs1 & ind.obs2)/sum(ind.obs1 | ind.obs2)
    return(val)
  }
  
ext_range <- function(x,season,r){
    val <- matrix(NA,ncol=3,nrow=length(x))
    season.types = c("Winter" ,"Spring" ,"Summer" ,"Fall")
    ind = sapply(season,function(x){which(season.types==x)})
    fun <- function(id){
      solve.h.BR(par=estimates.list[[r]][ind[id],],temp=x[id],logval=TRUE)
    }
    val = mapply(fun,id=1:length(x),SIMPLIFY = TRUE)
    return(val)
  }
  
for(r in 1:2){
    for(i in 1:4){
      data.mean.obs.hist[[r]][[i]]$range <- ext_range(data.mean.obs.hist[[r]][[i]]$x,data.mean.obs.hist[[r]][[i]]$season,r)   
      data.mean.tep.hist[[r]][[i]]$range <- ext_range(data.mean.tep.hist[[r]][[i]]$x,data.mean.tep.hist[[r]][[i]]$season,r)   
      data.mean.tep.245[[r]][[i]]$range <- ext_range(data.mean.tep.245[[r]][[i]]$x,data.mean.tep.245[[r]][[i]]$season,r)   
      data.mean.tep.585[[r]][[i]]$range <- ext_range(data.mean.tep.585[[r]][[i]]$x,data.mean.tep.585[[r]][[i]]$season,r)   
  }}
  
  save(data.mean.obs.hist.day,data.mean.tep.hist.day,data.mean.tep.245.day,data.mean.tep.585.day,data.mean.obs.hist,data.mean.tep.245,data.mean.tep.585,data.mean.tep.hist,file="data/data_mean_tep_10_moving.Rdata")
  
## plot the tail-correlation range ##
  model.selected <- c("AWI","MIROC","NorESM","AVG")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  load("data/temperature_covaraites.Rdata")
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  count = 1;p.list <- list()
  for(r in 1:2){
    for(s in 1:4){
        data <- rbind(data.mean.obs.hist[[r]][[s]],data.mean.tep.245[[r]][[s]],data.mean.tep.585[[r]][[s]])
        data$type = c(rep("Obs",nrow(data.mean.obs.hist[[r]][[s]])),
                      rep("SSP 2-4.5",nrow(data.mean.tep.245[[r]][[s]])),
                      rep("SSP 5-8.5",nrow(data.mean.tep.585[[r]][[s]])))
        data$Group = as.factor(paste(data$type,data$season))
        data$type = as.factor(data$type)
        p <- ggplot(data,aes(x=year,y=range)) + geom_line(aes(color=season,group=Group,linetype=type),size=0.6,alpha=1) 
        p <- p + scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
        p <- p + scale_x_continuous(breaks = seq(1965,2100,20),labels = seq(1965,2100,20),expand = c(0, 0), limits = c(1965,2100)) + xlab("Year") + ylab ("Logarithmic range")
        p <- p + labs(linetype="Group",color="Season") + guides(linetype=guide_legend(text=c("Obs","SSP 2-4.5","SSP 5-8.5"))) 
        p <- p + ggtitle(paste0(lab.regions[r],"; ",model.selected[s]))
        p <- p + theme(axis.text = element_text(size=10), 
                       axis.title.x = element_text(size=14), 
                       axis.title.y = element_text(size=14),
                       plot.title = element_text(hjust = 0.5),
                       panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                   color = "black",            # Color of the border
                                                   size = 1),
                       panel.background = element_rect(fill = "transparent"))
        #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
        p.list[[count]] <- p;count = count + 1
    }
  }
  
  pdf("figures/extreme_range_10_moving_avg.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:length(p.list)){
    show(p.list[i])
  }
  dev.off()
  
  
### plot the covariates daily ### 
################################
  model.selected <- c("AWI","MIROC","NorESM","AVG")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  count = 1;p.list <- list()
  for(r in 1:2){
    for(s in 1:4){
      data <- rbind(data.mean.obs.hist.day[[r]][[s]],data.mean.tep.245.day[[r]][[s]],data.mean.tep.585.day[[r]][[s]])
      data$type = c(rep("Obs",nrow(data.mean.obs.hist.day[[r]][[s]])),
                    rep("SSP 2-4.5",nrow(data.mean.tep.245.day[[r]][[s]])),
                    rep("SSP 5-8.5",nrow(data.mean.tep.585.day[[r]][[s]])))
      p <- ggplot(data,aes(x=date,y=tep,group=type,color=type)) + geom_line(alpha=0.7) 
      p <- p  + xlab("Year") + ylab ("Temperature (°C)")
      p <- p + labs(color='Group') 
      p <- p + scale_color_manual(values=hcl.colors(3, "Berlin")) 
      p <- p + ggtitle(paste0(lab.regions[r],"; ",model.selected[s]))
      p <- p + theme(axis.text = element_text(size=10), 
                     axis.title.x = element_text(size=14), 
                     axis.title.y = element_text(size=14),
                     plot.title = element_text(hjust = 0.5),
                     panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                 color = "black",            # Color of the border
                                                 size = 1),
                     panel.background = element_rect(fill = "transparent"))
      #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
      p.list[[count]] <- p;count = count + 1
    }
  }
  
  pdf("figures/temperature_covariate.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:8){
    show(p.list[i])
  }
  dev.off()
  
  ## plot the return level ##
  library(ggplot2)
  library(lubridate)
  library(mgcv)
  library(evgam)
  library(evd)

  getSeason <- function(d){
    y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date(paste0(y,"-3-20"),  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date(paste0(y,"-6-20"),  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date(paste0(y,"-9-20"),  format = "%Y-%m-%d") # Fall Equinox
    
    s.id<- ifelse (d > WS | d <= SE, "Winter",
                   ifelse (d > SE & d <= SS, "Spring",
                           ifelse (d > SS & d <= FE, "Summer", "Fall")))
    return(s.id)
  }

  getYear <- function(d){
	y = year(d)
    WS <- as.Date(paste0(y,"-12-20"), format = "%Y-%m-%d") # Winter Solstice
    YEND <- as.Date(paste0(y+1,"-1-1"),  format = "%Y-%m-%d")
    if(d > WS && d < YEND) return(y+1) else return(y)
 }
  
  load("data/data_pot.Rdata")
  load("data1/temperature_covaraites.Rdata")
  model.selected <- c("AWI","MIROC","NorESM","AVG")
  season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
  regions <- c("danube","mississippi")
  lab.regions <- c("Danube","Mississippi")
  y.thres=10
  count = 1
  p.list <- list()
  for(r in 1:2){
    load(paste0("data/data_pot_10_",regions[r],".Rdata"))
    for(s in 1:4){
    data = precip.ts.df[[r*2-1]][,-c(1:3)]
    data[data > 1500] = NA
    loc.ind <- which.max(colMeans(data,na.rm = T))
    ### prepare the data frame to predict ###
    data <- rbind(data.mean.obs.hist.day[[r]][[s]],data.mean.tep.245.day[[r]][[s]],data.mean.tep.585.day[[r]][[s]])
    data$type = c(rep("Obs",nrow(data.mean.obs.hist.day[[r]][[s]])),
                  rep("SSP 2-4.5",nrow(data.mean.tep.245.day[[r]][[s]])),
                  rep("SSP 5-8.5",nrow(data.mean.tep.585.day[[r]][[s]])))
    days <- data$date
    d.ind = yday(days)
    y.ind = year(days);y.ind = y.ind - y.ind[1] + 1
    alt <- station.df[[r*2-1]]$elev[loc.ind]/1000
    lon <- station.df[[r*2-1]]$Y[loc.ind]
    lat <- station.df[[r*2-1]]$X[loc.ind]
    main = paste0(lab.regions[r],"; ",model.selected[s],"; ","Station ",loc.ind," (", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
    D = length(lat);Dt = length(days)
    data.pred <- data.frame(temp = rep(data$tep,times=D),
                            day = rep(d.ind,times=D),
                            year = rep(y.ind,times=D),
                            alt = rep(alt,each=Dt),
                            lon = rep(lon,each=Dt),
                            lat = rep(lat,each=Dt),
                            col=rep(1:D,each=Dt),
                            row=rep(1:Dt,times=D))
    ind = !is.na(data.pred$temp);data.pred <- data.pred[ind,]
    mean.pred  <- exp(predict.gam(results.gam,newdata = data.pred))
    sig2.pred <- results.gam$sig2
    shape.pred = 1/sig2.pred;scale.pred <- mean.pred/shape.pred;
    quantile.pred <- qgamma(0.9,shape = shape.pred,scale=scale.pred)
    data.pred$est.quantile = quantile.pred
    prob.exceed.pred <- predict.gam(results.bin,newdata=data.pred,type="response")
    gpd.pred <- predict(results.gpd,newdata=data.pred)
    scale.gpd.pred = exp(gpd.pred[,1]);shape.gpd.pred = gpd.pred[1,2]
    
    return.level = (1-1/(100*365)) 
    prob.gpd <- (return.level - (1-prob.exceed.pred))/prob.exceed.pred
    
    return.value <- y.thres + quantile.pred + qgpd(prob.gpd,loc=0,scale=scale.gpd.pred,shape=shape.gpd.pred)
    by.list = list(season=getSeason(days)[ind],year=sapply(days,getYear)[ind],type=data$type[ind])
    data <-aggregate(return.value,by=by.list,FUN=mean)
	data$season <- as.factor(data$season)
	data$type <- as.factor(data$type)
	data$x[data$type == "Obs" & data$season == "Winter" & data$year == 2016] = NA
    p <- ggplot(data,aes(x=year,y=x,color=season,linetype=type)) + geom_line(alpha=1,size=0.6) 
    p <- p + scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
    p <- p +  xlab("Year") + ylab ("Return level (mm)") 
    p <- p + labs(color="Season",linetype="Group") 
    #p <- p + scale_color_manual(values=hcl.colors(4, "Berlin")) 
    p <- p + ggtitle(main)
    p <- p + theme(axis.text = element_text(size=10), 
                   axis.title.x = element_text(size=14), 
                  axis.title.y = element_text(size=14),
                   plot.title = element_text(hjust = 0.5),
                   panel.border = element_rect(fill = "transparent", # Needed to add the border
                                               color = "black",            # Color of the border
                                               size = 1),
                   panel.background = element_rect(fill = "transparent")) 
    
    #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
    p.list[[count]] <- p;count = count + 1
    print(count)
    }
    }

  pdf("figures/return_level_margins.pdf",width = 6,height = 3,onefile = TRUE)
  for(i in 1:8){
    show(p.list[i])
  }
  dev.off()

### heatmap for exploratory analysis ###
##########################################      
  library(ggplot2)
  library(lubridate)
  load("data/data_pot.Rdata",e.temp<-new.env())
  precip.ts.df <- e.temp$precip.ts.df
  station.df <- e.temp$station.df
  date.ts <- e.temp$date.ts
  elev.danube <- station.df[[1]]$elev
  D = length(elev.danube)
  Dt = length(date.ts)
  data.danube <- data.frame(precip = drop(as.vector(as.matrix(precip.ts.df[[1]][,-c(1:3)][,order(elev.danube)]))),
                            date=rep(date.ts,D),
                            year=rep(year(date.ts),D),
                            yday=rep(yday(date.ts),D),
                            id = rep(1:D,each=Dt))
  
  elev.miss <- station.df[[3]]$elev
  D = length(elev.miss)
  Dt = length(date.ts)
  data.miss <- data.frame(precip = drop(as.vector(as.matrix(precip.ts.df[[3]][,-c(1:3)][,order(elev.miss)]))),
                            date=rep(date.ts,D),
                            year=rep(year(date.ts),D),
                            yday=rep(yday(date.ts),D),
                            id = rep(1:D,each=Dt))
  data.miss <- data.miss[is.na(data.miss$precip) | (data.miss$precip < 1500 & data.miss$precip >=0),]
  data.danube <- data.danube[is.na(data.danube$precip) | (data.danube$precip < 1500 & data.danube$precip >=0), ]
  rm(e.temp,precip.ts.df)  
  data <- list(data.danube,data.miss)
  elev <- list(elev.danube[order(elev.danube)],elev.miss[order(elev.miss)])
  rm(data.miss,data.danube)
  mean.func <- function(x,y=TRUE){
    n <- sum(!is.na(x))
    if(y & n >= 52){
      return(mean(x,na.rm=TRUE))
    }
    if(!y & n >= 10){
      return(mean(x,na.rm=TRUE))
    }
    return(NA)
  }
  
  lab.regions <- c("Danube","Mississippi")

  ## prepare for the data.frame ##
  library(viridis)
  p1.list <- p2.list <- list()
  for(r in 1:2){
  mean.year <- aggregate(data[[r]]$precip,by=list(data[[r]]$year,data[[r]]$id),FUN=mean.func,y=TRUE)
  mean.day <- aggregate(data[[r]]$precip,by=list(data[[r]]$yday,data[[r]]$id),FUN=mean.func,y=FALSE)  
  mean.day$elev <- rep(elev[[r]],each=366)
  names(mean.day) <- c("day","ID","precip","elev")
  mean.year$elev <- rep(elev[[r]],each=51)
  names(mean.year) <- c("year","ID","precip","elev")
  color_breaks = c(0:10,max(mean.year$precip,na.rm=T))
  main = paste0("Yearly average: ",lab.regions[r])
  p1 <- ggplot(mean.year) + geom_tile(aes(x=year,y=ID,fill=precip)) 
  p1 <- p1 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,length(elev[[r]]),length.out=5),labels=as.character(elev[[r]][seq(1,length(elev[[r]]),length.out=5)],limits=c(0,length(elev[[r]]))))
  p1 <- p1 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev = TRUE),limits=c(0,max(mean.year$precip,na.rm=T)),
                values=scales::rescale(color_breaks, to = c(0,1), 
                      from = c(0,max(mean.year$precip,na.rm=T)))) 
  p1 <- p1 + ggtitle(main) + xlab("Year") + 
        theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = "transparent", # Needed to add the border
                                      color = "transparent",            # Color of the border
                                      size = 0.5))
  color_breaks = c(0:10,max(mean.day$precip,na.rm=T))
  main = paste0("Daily average: ",lab.regions[r])
  p2 <- ggplot(mean.day) + geom_tile(aes(x=day,y=ID,fill=precip)) 
  p2 <- p2 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,length(elev[[r]]),length.out=5),labels=as.character(elev[[r]][seq(1,length(elev[[r]]),length.out=5)],limits=c(0,length(elev[[r]]))))
  p2 <- p2 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev=TRUE),limits=c(0,max(mean.day$precip,na.rm=T)),
                  values=scales::rescale(color_breaks, to = c(0,1), 
                                    from = c(0, max(mean.day$precip,na.rm=T)))) 
  p2 <- p2 + ggtitle(main) + xlab("Day") + 
    theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = "transparent", # Needed to add the border
                                      color = "transparent",            # Color of the border
                                      size = 0.5))
  p1.list[[r]] <- p1; p2.list[[r]] <- p2
  }

  
  pdf("figures/heatmaps.pdf",width = 6,height = 5,onefile = TRUE)
  for(i in 1:2){
    show(p1.list[[i]])
    show(p2.list[[i]])
  }
  dev.off()  
  
  save.image("temp.RData")
