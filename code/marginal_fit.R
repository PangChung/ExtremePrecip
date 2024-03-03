rm(list=ls())
args <- commandArgs(TRUE)
load("data/precip.RData")
load("data/temperature.RData")
library(parallel)
library(mgcv)
library(ggplot2)
library(lubridate)
library(evgam)
library(evd)
ls()
idx = 4;model.id = 1
for (arg in args) eval(parse(text = arg))
## prepare the dataframe for marginal fit ##
# for(idx in 1:length(region.id)){
print(idx)
## prepare the time covariate: year and date
ind.data = date.df$date >= START.date & date.df$date <= END.date
y.ind = date.df$year[ind.data];y.ind = y.ind - y.ind[1] + 1
d.ind = yday(date.df$date[ind.data])
ind.station = station$group.id==region.id[idx]
alt <- station$elev[ind.station]/1000 ## elevation of the station
lon <- station$Y[ind.station]
lat <- station$X[ind.station]

## compute the consective temperature averages 
tep.covariate <- temperature.covariate[[idx]][ind.data]
D = sum(ind.station);Dt = sum(ind.data) # dimensions
y = unlist(precip[[idx]])
y.thres <- 10;y = y-y.thres ## remove the values that are below y.thres
print(quantile(y[!is.na(y) & y>0 & y<2000],0.99))
#}
## generate the data frame for marginal fitting ## 
data.df <- data.frame(y=y,temp = rep(tep.covariate,times=D),
                      day = rep(d.ind,times=D),
                      year = rep(y.ind,times=D),
                      alt = rep(alt,each=Dt),
                      lon = rep(lon,each=Dt),
                      lat = rep(lat,each=Dt),
                      col=rep(1:D,each=Dt),
                      row=rep(1:Dt,times=D))[!is.na(y) & y > 0 & y < 1000,]
data.df = data.df[complete.cases(data.df),] ## select the complete dataframe

## start fitting the marginal model ##
## WARNING! Long time to run ##
t0 <- proc.time()
message("start Gamma fitting")
thres.prob = 0.8
if(model.id==1){
    formula = y ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
    results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
    est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
    est.shape = 1/est.sig2;est.scale <- est.mean/est.shape
    data.df$est.quantile <- qgamma(thres.prob,shape = est.shape,scale=est.scale)
    data.df$est.shape = est.shape;data.df$est.scale = est.scale
    est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
}else{
    formula = list(log(y) ~ temp + s(alt,k=10) + s(day,k=10) + ti(lon,lat,k=10) ,~s(day,k=10)+s(alt,k=10))
    results.gam = gam(formula,family=gaulss(link=list("log","logb")),data=data.df)
    est.gam = fitted(results.gam)
    data.df$est.mean = est.gam[,1];data.df$est.sd = 1/est.gam[,2]
    data.df$est.quantile = exp(qnorm(thres.prob,mean=data.df$est.mean,sd=data.df$est.sd))
    data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
    est.prob <- pnorm(log(data.df$y),mean=data.df$est.mean,sd=data.df$est.sd)
}
data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
data.df$est.prob <- est.prob 
print(summary(est.prob[data.df$y.bin] - thres.prob))
message("start binominal fitting")
formula.bin = y.bin ~ temp + s(day,k=10) + s(alt,k=10) + ti(lon,lat,k=10)
results.bin <- gam(formula.bin,family = binomial(link="logit"),data=data.df)
est.prob.exceed <- fitted(results.bin) ## fitted exceeding probability
data.df$est.prob.exceed <- est.prob.exceed

message("start GPD fitting")
data.df$y.gpd <- data.df$y - data.df$est.quantile
data.df.gpd <- data.df[data.df$y.gpd>0,]
formula.gpd = list(y.gpd ~ log(est.quantile) + temp + s(day,k=10) + s(alt,k=10) + ti(lon,lat,k=10),~1)
results.gpd <- evgam(formula.gpd,data=data.df.gpd,family="gpd")
est.scale.gpd = exp(fitted(results.gpd)[,1]);est.shape.gpd = fitted(results.gpd)[1,2]
data.df.gpd$est.scale.gpd = est.scale.gpd;data.df.gpd$est.shape.gpd = est.shape.gpd

## transform the data to pseudo uniform scores  ##
est.prob[data.df$y.bin] <- 1 - est.prob.exceed[data.df$y.bin] + 
est.prob.exceed[data.df$y.bin]*pgpd(data.df.gpd$y.gpd,loc = 0,scale = est.scale.gpd,shape = est.shape.gpd)


U <- matrix(NA,nrow=Dt,ncol=D)
U[cbind(data.df$row,data.df$col)] <- est.prob/(1+1e-10) ## avoid computational issues

t1 <- proc.time() - t0

save(results.gam,results.gpd,results.bin,data.df,data.df.gpd,U,file=paste0("data/marginal_fit_",idx,"_model_",model.id,".RData"))




