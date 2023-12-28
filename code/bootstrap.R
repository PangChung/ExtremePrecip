rm(list=ls())
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

idx.region = 2
norm.ind = 1;seaon.ind = 1
init = c(-1.5,4,0);fixed=c(F,F,F)
bootstrap.ind = 1;season.ind=1
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
for (arg in args) eval(parse(text = arg))

load(paste0("data/marginal_fit_",idx.region,".RData"),e<-new.env())
data.df <- e$data.df
rm(e)


message("start Gamma fitting")
formula = y ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
est.shape = 1/est.sig2;est.scale <- est.mean/est.shape
data.df$est.quantile <- qgamma(0.9,shape = est.shape,scale=est.scale)
data.df$est.shape = est.shape;data.df$est.scale = est.scale

message("start binominal fitting")
data.df$y.bin <- as.numeric(data.df$y > data.df$est.quantile)
formula.bin = y.bin ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10)
results.bin <- gam(formula.bin,family = binomial(link="logit"),data=data.df)
est.prob.exceed <- fitted(results.bin) ## fitted exceeding probability
data.df$est.prob.exceed <- est.prob.exceed

message("start GPD fitting")
data.df$y.gpd <- data.df$y - data.df$est.quantile
data.df.gpd <- data.df[data.df$y.gpd>0,]
formula.gpd = list(y.gpd ~ temp + s(day,k=10) + ti(lon,lat,k=10) + s(alt,k=10),~1)
results.gpd <- evgam(formula.gpd,data=data.df.gpd,family="gpd")
est.scale.gpd = exp(fitted(results.gpd)[,1]);est.shape.gpd = fitted(results.gpd)[1,2]
data.df.gpd$est.scale.gpd = est.scale.gpd;data.df.gpd$est.shape.gpd = est.shape.gpd

## transform the data to pseudo uniform scores  ##
est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
head(est.prob[data.df$y.bin])
est.prob[data.df$y.bin] <- 1 - est.prob.exceed[data.df$y.bin] + 
est.prob.exceed[data.df$y.bin]*pgpd(data.df.gpd$y.gpd,loc = 0,scale = est.scale.gpd,shape = est.shape.gpd)
data.df$est.prob <- est.prob

U <- matrix(NA,nrow=Dt,ncol=D)
U[cbind(data.df$row,data.df$col)] <- est.prob/(1+1e-10) ## avoid computational issues

save(results.gam,results.gpd,results.bin,data.df,data.df.gpd,U,file=paste0("data/marginal_fit_",idx,".RData"))

