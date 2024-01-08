rm(list=ls())
args <- commandArgs(TRUE)
library(parallel)
library(lubridate)
library(mvPotST)
library(evd)
source("code/utility.R")
load("data/precip.RData")
load("data/transformed_coordinates.RData")
load("data/temperature.RData")

idx.region = 2;ST = TRUE
init = c(-1.5,4,0);fixed=c(F,F,F)
bootstrap=FALSE;
bootstrap.idx = 1;season.idx=1
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
norm.idx = 1
for (arg in args) eval(parse(text = arg))

## file where the data should be stored ##
file = paste0("data/fit_pot_ST_",season.idx,"_",region.name[idx.region],"_",norm.ind,".Rdata") 
if(file.exists(file)){
    load(file,e<-new.env())
    init = e$result$par
    rm(e)
    #stop("file already exists")
} 
ncores = 5 #detectCores()
## load the data from marginal fit ##
load(paste0("data/marginal_fit_",idx.region,"_model_1.RData"))
## choose the r risk functional...##
if(norm.ind==1){
    est.shape.gpd <- data.df.gpd$est.shape.gpd[1]
}else{
    est.shape.gpd <- 1
}
# Define locations 
loc = loc.trans.list[[idx.region]]
idx.season = date.df$season == season[season.idx]
idx.season = idx.season[date.df$date >= START.date & date.df$date <= END.date]

## load the observations and covariates ## 
obs = split(U,row(U))
obs = subset(obs,idx.season)
no.obs = sapply(obs,function(x){sum(!is.na(x))})
obs[no.obs>1] = lapply(obs[no.obs>1],evd::qgpd,shape=1,loc=1,scale = 1)
reg.t = temperature[[idx.region]][idx.season]
reg.t = (reg.t - mean(reg.t))/sd(reg.t)

r.obs <- suppressWarnings(unlist(lapply(obs,function(x){if(sum(!is.na(x))!=0){rFun(x[!is.na(x)],u=1,
est.shape.gpd)}else{NA}})))

thres = quantile(r.obs[no.obs > 5],0.95,na.rm=TRUE)

## select the exceedances
idx.exc = no.obs > 5 & r.obs > thres 
stopifnot( sum(idx.exc) > 0 )

reg.t = reg.t[idx.exc]
reg.t = (reg.t - mean(reg.t))/sd(reg.t) 
exceedances <- obs[idx.exc]

t0 <- proc.time()
result = fit.gradientScoreBR(obs=exceedances,loc=loc,init=init + runif(3)*0.1,fixed = fixed,vario = vario,u = thres,ST = ST,nCores = ncores,weightFun = weightFun,dWeightFun = dWeightFun)
t1 <- proc.time() - t0
print(t1)

save(idx.exc,result,exceedances,reg.t,est.shape.gpd,thres,file=file)



