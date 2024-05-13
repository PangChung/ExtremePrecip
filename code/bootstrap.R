rm(list=ls())
t0 <- proc.time()
args <- commandArgs(TRUE)
source("code/utility.R")
load("data/precip.RData")
load("data/transformed_coordinates.RData")
load("data/temperature.RData")
idx.region = 1;bootstrap.ind = 1;njobs=1
init = c(-2,log(100),0);fixed=c(F,F,F)
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
computer="local"
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/ExtremePrecip";.libPaths("/home/z3536974/R/x86_64-pc-linux-gnu-library/4.2")},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/ExtremePrecip"}
)
library(parallel)
library(lubridate)
library(mvPotST)
library(evd)
library(mgcv)
library(evgam)

file.marginal = paste0(DataPath,"/data/marginal_fit_1_",bootstrap.ind,"_",idx.region,".RData")
file2save = paste0(DataPath,"/data/fit_bootstrap_1_",bootstrap.ind,"_",idx.region,".RData")
if(file.exists(file2save)) {stop("fit is already done")}
ncores = floor(detectCores()/njobs)
if(idx.region==2 | idx.region==7){ncores=ncores*2}
init.seed = as.integer((as.integer(Sys.time()) + sample.int(10^5,1))%%10^5)
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
if(bootstrap.ind!=0){
    idx_samples = sort(sample(1:1000,1000,replace=TRUE))
} else {
    idx_samples = 1:1000
}
ind.sample = sort(unlist(lapply(idx_samples,function(x){which(idx_numbers==x)})))
ind.station = station$group.id==region.id[idx.region]
alt <- station$elev[ind.station]/1000 ## elevation of the station
lon <- station$Y[ind.station]
lat <- station$X[ind.station]

## compute the consective temperature averages ## 
D = sum(ind.station);Dt = length(ind.sample) # dimensions
y = unlist(lapply(precip[[idx.region]],function(x){x[ind.sample]}))
y.thres <- 0;y = y-y.thres ## remove the values that are below y.thres

## generate the data frame for marginal fitting ## 
data.df <- data.frame(y=y,
                    temp = rep(tep.covariate[ind.sample],time=D),
                    date = rep(date.ind[ind.sample],time=D),
                    day = rep(d.ind[ind.sample],time=D),
                    year = rep(y.ind[ind.sample],time=D),
                    alt = rep(alt,each=Dt),
                    lon = rep(lon,each=Dt),
                    lat = rep(lat,each=Dt),
                    col=rep(1:D,each=Dt),
                    row=rep((1:sum(ind.data))[ind.sample],time=D))[!is.na(y) & y > 0 & y < 1000,]
data.df = data.df[complete.cases(data.df),] ## select the complete dataframe

if(file.exists(file.marginal)){
    load(file.marginal)
}else{
    message("start Gamma fitting")
    thres.prob = 0.9
    formula = y ~ temp + s(day,k=10) + te(lon,lat,k=10) + s(alt,k=10)
    results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
    est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
    est.shape = 1/est.sig2;est.scale <- est.mean/est.shape
    data.df$est.quantile <- qgamma(thres.prob,shape = est.shape,scale=est.scale)
    data.df$est.shape = est.shape;data.df$est.scale = est.scale
    est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
    data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
    data.df$est.prob <- est.prob 

    message("start binominal fitting")
    formula.bin = y.bin ~ temp + s(day,k=10) + s(alt,k=10) + te(lon,lat,k=10)
    results.bin <- gam(formula.bin,family = binomial(link="logit"),data=data.df)
    est.prob.exceed <- fitted(results.bin) ## fitted exceeding probability
    data.df$est.prob.exceed <- est.prob.exceed

    message("start GPD fitting")
    data.df$y.gpd <- data.df$y - data.df$est.quantile
    data.df.gpd <- data.df[data.df$y.gpd>0,]
    formula.gpd = list(y.gpd ~ log(est.quantile) + temp + s(day,k=10) + s(alt,k=10) + te(lon,lat,k=10),~1)
    results.gpd <- evgam(formula.gpd,data=data.df.gpd,family="gpd")
    est.scale.gpd = exp(fitted(results.gpd)[,1]);est.shape.gpd = fitted(results.gpd)[1,2]
    data.df.gpd$est.scale.gpd = est.scale.gpd;data.df.gpd$est.shape.gpd = est.shape.gpd

    ## transform the data to pseudo uniform scores  ##
    est.prob[data.df$y.bin] <- 1 - est.prob.exceed[data.df$y.bin] + 
    est.prob.exceed[data.df$y.bin]*pgpd(data.df.gpd$y.gpd,loc = 0,scale = est.scale.gpd,shape = est.shape.gpd)

    U <- matrix(NA,nrow=sum(ind.data),ncol=D)
    U[cbind(data.df$row,data.df$col)] <- est.prob/(1+1e-10) ## avoid computational issues
    
    if(bootstrap.ind==0) save(results.bin,results.gam,results.gpd,data.df.gpd,U,data.df,file=file.marginal)
}

# depdence fit ##
result.list <- list(list(),list())
for(count in 1:8){
        norm.ind = (count-1) %/% 4 + 1
        season.idx = (count - 1) %% 4 + 1
        file = paste0(DataPath,"/data/fit_bootstrap_1_0_",idx.region,".RData") 
        # file2save = paste0("/srv/scratch/z3536974/data/fit_bootstrap_",bootstrap.ind,"_",idx.region,"_",count,".RData")
        if(file.exists(file)){
            load(file,e<-new.env())
            init = e$result.list[[norm.ind]][[season.idx]]$par
            rm(e)
        }
        ## choose the r risk functional...##
        if(norm.ind==1){
            est.shape.gpd <- data.df.gpd$est.shape.gpd[1]
        }else{
            est.shape.gpd <- 1
        }
        # Define locations 
        loc = loc.trans.list[[idx.region]]/1000
        loc = apply(loc,2,function(x){x - mean(x,na.rm=TRUE)})
        date.df2 = date.df[ind.data,][ind.sample,]
        idx.season = (date.df2$season == season[season.idx])
        ## load the observations and covariates ## 
        obs = split(U,row(U))[ind.sample]
        obs = subset(obs,idx.season)
        no.obs = sapply(obs,function(x){sum(!is.na(x))})
        obs[no.obs>1] = lapply(obs[no.obs>1],evd::qgpd,shape=1,loc=1,scale = 1)
        reg.t = temperature.covariate[[idx.region]][ind.data][ind.sample][idx.season]

        r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],u=1,est.shape.gpd)}else{NA}})))
        thres = quantile(r.obs[no.obs > 10],0.9,na.rm=TRUE)

        ## select the exceedances
        idx.exc = no.obs > 10 & r.obs > thres 
        stopifnot( sum(idx.exc) > 0 & any(!is.na(reg.t)) )
        print(sum(idx.exc))
        reg.t = reg.t[idx.exc]
        exceedances <- obs[idx.exc]

        result.list[[norm.ind]][[season.idx]] = fit.gradientScoreBR(obs=exceedances,loc=loc,init=init,fixed = fixed,vario = vario,u = thres,method="Nelder-Mead",ST = TRUE,nCores = ncores,weightFun = weightFun,dWeightFun = dWeightFun)
        
}

t1 = proc.time() - t0
print(t1)
save(t1,result.list,file=file2save)

