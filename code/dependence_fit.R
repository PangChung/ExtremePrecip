rm(list=ls())
args <- commandArgs(TRUE)
library(parallel)
library(ggplot2)
library(lubridate)
library(mvPotST)
library(evd)
#library(fields)
source("code/utility.R")
load("data/precip.RData")
load("data/transformed_coordinates.RData")
load("data/temperature.RData")
#load("data/temperature_pred.RData")

idx.region = 2;ST = TRUE
init = c(-1.5,4,0);fixed=c(F,F,F)
bootstrap=FALSE;
bootstrap.ind = 1;season.ind=1
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
norm.ind = 1;seaon.ind = 1
for (arg in args) eval(parse(text = arg))

## load the data from marginal fit ##
load(paste0("data/marginal_fit_",idx.region,".RData"))

## file where the data should be stored ##
file = paste0("data/fit_pot_ST_",season.ind,"_",region.name[idx.region],"_",norm.ind,".Rdata") 
if(file.exists(file)){stop("file already exists")} 
ncores = 4 #detectCores()

## choose the r risk functional...##
if(norm.ind==1){
    est.shape.gpd <- data.df.gpd$est.shape.gpd[1]
}else{
    est.shape.gpd <- 1
}

# Define locations 
loc = loc.trans.list[[idx.region]]
idx.season = date.df$season == season[season.ind]
idx.season = idx.season[date.df$date >= START.date & date.df$date <= END.date]

## load the observations and covariates ## 
obs = split(U,row(U))
obs = subset(obs,idx.season)
no.obs = sapply(obs,function(x){sum(!is.na(x))})
obs[no.obs>1] = lapply(obs[no.obs>1],evd::qgpd,shape=1,loc=1,scale = 1)
reg.t = temperature[[idx.region]][idx.season]

r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],u=1,
est.shape.gpd)}else{NA}})))

thres = quantile(r.obs[no.obs > 5],0.95,na.rm=TRUE)

## select the exceedances
idx.exc = no.obs > 5 & r.obs > thres 
stopifnot( sum(idx.exc) > 0 )

reg.t = reg.t[idx.exc]
reg.t = (reg.t - mean(reg.t))/sd(reg.t) 
exceedances <- obs[idx.exc]

if(bootstrap){
    while(!file.exists(file)){Sys.sleep(60)}
    load(file,e<-new.env())
    
    param = e$result$par
    init = param
    rm(e)
    init.seed = as.integer((as.integer(Sys.time())/bootstrap.ind + sample.int(10^5,1))%%10^5)
    set.seed(init.seed)
    boot.index = sample(1:length(exceedances),length(exceedances),replace = TRUE)
    
    exceedances = exceedances[boot.index]
    reg.t = reg.t[boot.index]
    file = paste0("data/fit_pot_ST_bootstrap_",init.seed,"_",season.ind,"_",regions[idx.region],".Rdata")
}

t0 <- proc.time()
result = fit.gradientScoreBR(obs=exceedances,loc=loc,init=init + runif(3)*0.1,fixed = fixed,vario = vario,u = thres,ST = ST,nCores = ncores,weightFun = weightFun,dWeightFun = dWeightFun)
t1 <- proc.time() - t0
print(t1)

save(idx.exc,result,exceedances,reg.t,est.shape.gpd,thres,file=file)



