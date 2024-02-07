rm(list=ls())
t0 <- proc.time()
args <- commandArgs(TRUE)
library(parallel)
library(lubridate)
library(mvPotST)
library(evd)
library(mgcv)
library(ggplot2)
library(evgam)
source("code/utility.R")
load("data/precip.RData")
load("data/transformed_coordinates.RData")
load("data/temperature.RData")

idx.region = 1;bootstrap.ind = 400
init = c(-4,4,0);fixed=c(F,F,F)
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")

for (arg in args) eval(parse(text = arg))

ncores = 10 #detectCores()
init.seed = as.integer((as.integer(Sys.time())/bootstrap.ind + sample.int(10^5,1))%%10^5)
set.seed(init.seed)
## prepare the time covariate: year and date
ind.data = date.df$date >= START.date & date.df$date <= END.date
y.ind = date.df$year[ind.data];y.ind=y.ind - y.ind[1] + 1
d.ind = yday(date.df$date[ind.data])  
date.ind = date.df$date[ind.data]
tep.covariate <- temperature.covariate[[idx.region]][ind.data]

## sample the data ##
idx_numbers = sort(rep(1:1000,length.out = sum(ind.data)))
#idx_samples = (1:300)[-bootstrap.ind]
idx_samples = sort(sample(1:1000,1000,replace=TRUE))
ind.sample = sort(unlist(lapply(idx_samples,function(x){which(idx_numbers==x)})))

ind.station = station$group.id==region.id[idx.region]
alt <- station$elev[ind.station]/1000 ## elevation of the station
lon <- station$Y[ind.station]
lat <- station$X[ind.station]

## compute the consective temperature averages ## 
D = sum(ind.station);Dt = length(ind.sample) # dimensions
y = unlist(lapply(precip[[idx.region]],function(x){x[ind.sample]}))
y.thres <- 10;y = y-y.thres ## remove the values that are below y.thres



## generate the data frame for marginal fitting ## 
data.df <- data.frame(y=y,
                    temp = rep(tep.covariate[ind.sample],each=D),
                    date = rep(date.ind[ind.sample],each=D),
                    day = rep(d.ind[ind.sample],each=D),
                    year = rep(y.ind[ind.sample],each=D),
                    alt = rep(alt,time=Dt),
                    lon = rep(lon,time=Dt),
                    lat = rep(lat,time=Dt),
                    col=rep(1:D,time=Dt),
                    row=rep((1:sum(ind.data))[ind.sample],each=D))[!is.na(y) & y > 0 & y < 1000,]

data.df = data.df[complete.cases(data.df),] ## select the complete dataframe


message("start Gamma fitting")
thres.prob = 0.8
formula = y ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
est.shape = 1/est.sig2;est.scale <- est.mean/est.shape
data.df$est.quantile <- qgamma(thres.prob,shape = est.shape,scale=est.scale)
data.df$est.shape = est.shape;data.df$est.scale = est.scale
est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
data.df$est.prob <- est.prob 

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

U <- matrix(NA,nrow=sum(ind.data),ncol=D)
U[cbind(data.df$row,data.df$col)] <- est.prob/(1+1e-10) ## avoid computational issues

## depdence fit ##
result.list <- list(list(),list())
for(norm.ind in 1:2){
     for(season.idx in 1:4){
        file = paste0("data/fit_pot_ST_",season.idx,"_",idx.region,"_",norm.ind,".Rdata") 
        if(file.exists(file)){
            load(file,e<-new.env())
            init = e$result$par
            rm(e)
            #stop("file already exists")3
        }
        ## choose the r risk functional...##
        if(norm.ind==1){
            est.shape.gpd <- data.df.gpd$est.shape.gpd[1]
        }else{
            est.shape.gpd <- 1
        }
        # Define locations 
        loc = loc.trans.list[[idx.region]]/1000
        date.df2 = date.df[ind.data,][ind.sample,]
        idx.season = date.df2$season == season[season.idx]
        ## load the observations and covariates ## 
        obs = split(U,row(U))
        obs = subset(obs,idx.season)
        no.obs = sapply(obs,function(x){sum(!is.na(x))})
        obs[no.obs>1] = lapply(obs[no.obs>1],evd::qgpd,shape=1,loc=1,scale = 1)
        reg.t = temperature[[idx.region]][ind.data][ind.sample][idx.season]
        reg.t = (reg.t - mean(reg.t))/sd(reg.t)

        r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],u=1,
        est.shape.gpd)}else{NA}})))

        thres = quantile(r.obs[no.obs > 5],0.95,na.rm=TRUE)

        ## select the exceedances
        idx.exc = no.obs > 5 & r.obs > thres 
        stopifnot( sum(idx.exc) > 0 & any(!is.na(reg.t)) )

        reg.t = reg.t[idx.exc]
        exceedances <- obs[idx.exc]

        result.list[[norm.ind]][[season.idx]] = fit.gradientScoreBR(obs=exceedances,loc=loc,init=init,fixed = fixed,vario = vario,u = thres,method="BFGS",ST = TRUE,nCores = ncores,weightFun = weightFun,dWeightFun = dWeightFun)
    }
}

file2save = paste0("data/fit_bootstrap_",bootstrap.ind,"_",idx.region,".RData")
save(result.list,file=file2save)


