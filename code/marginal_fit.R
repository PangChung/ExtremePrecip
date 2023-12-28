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

hist(data.df$y,50)
## start fitting the marginal model ##
## WARNING! Long time to run ##
t0 <- proc.time()
message("start Gamma fitting")
thres.prob = 0.8
if(model.id==1){
    formula = y ~ temp + s(day,k=5) + ti(lon,lat,k=5) + s(alt,k=5)
    results.gam = gam(formula,family=Gamma(link="log"),data=data.df)
    est.sig2 <- results.gam$sig2;est.mean <- results.gam$fitted.values
    est.shape = 1/est.sig2;est.scale <- est.mean/est.shape
    data.df$est.quantile <- qgamma(thres.prob,shape = est.shape,scale=est.scale)
    data.df$est.shape = est.shape;data.df$est.scale = est.scale
    est.prob <- pgamma(data.df$y,shape=est.shape,scale=est.scale)
}else{
    formula = list(log(y) ~ temp + s(alt,k=5) + s(day,k=5) + ti(lon,lat,k=5) ,~s(day,k=5)+s(alt,k=5))
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
formula.bin = y.bin ~ temp + s(day,k=5) + s(alt,k=5) + ti(lon,lat,k=5) + ti(alt,)
results.bin <- gam(formula.bin,family = binomial(link="logit"),data=data.df)
est.prob.exceed <- fitted(results.bin) ## fitted exceeding probability
data.df$est.prob.exceed <- est.prob.exceed

message("start GPD fitting")
data.df$y.gpd <- data.df$y - data.df$est.quantile
data.df.gpd <- data.df[data.df$y.gpd>0,]
formula.gpd = list(y.gpd ~ temp + s(day,k=5) + s(alt,k=5) + ti(lon,lat,k=5),~1)
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

#}

## plot the qqplot for random location ##
#for(idx in 1:length(region.id)){
load(paste0("data/marginal_fit_",idx,"_model_",model.id,".RData"))
set.seed(1000)
idx.list = sample(1:sum(station$group.id==region.id[idx]),8,replace = F,prob=apply(U,2,function(x){sum(!is.na(x))}))
png(file = paste0("figures/qqplot_marginal_",idx,"_model_",model.id,".png"),height=4*3,width=4*3,units="cm",res=300, pointsize=6)
par(mfrow=c(3,3),mar=c(5,5,3,1),mgp=c(2.5,1,0),cex.lab=2,cex.axis=1.5,cex.main=2)
theoretical.quantiles <- U[,idx.list]
empirical.quantiles <- 1:1000/(1+1000)
theoretical.quantiles <- split(theoretical.quantiles,col(theoretical.quantiles))
theoretical.quantiles <- sapply(theoretical.quantiles,function(x) quantile(x,prob=empirical.quantiles,na.rm=T),simplify = F)
qqplot(unlist(theoretical.quantiles),empirical.quantiles,cex=0.5,pch=20,xlim=c(0,1), ylim=c(0,1),
    xlab="Theoretical quantiles",ylab="Empirical quantiles",main="All stations")
abline(0,1,col=2,lwd=2)
for(ind in 1:length(idx.list)){
qqplot(theoretical.quantiles[[ind]],empirical.quantiles,
       cex=0.5,pch=20,
       xlab="Theoretical quantiles",ylab="Empirical quantiles",main=paste("Station",idx.list[ind]))
abline(0,1,col=2,lwd=2)
}        
dev.off() 
#}


## the pesudo-uniform scores based on this marginal fit 
## is stored in the list U.data
# load("data/marginal_fit.RData")
# ## plot the marginal return level ##
#   model.selected <- c("AWI","MIROC","NorESM","AVG") #selected climate models with their averages
#   season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
#   regions <- c("danube","mississippi")
#   lab.regions <- c("Danube","Mississippi")
#   y.thres=10
#   count = 1
#   p.list <- list()
#   for(r in 1:2){
#     load(paste0("data/marginal_model_for_",regions[r],".RData"))
#     for(s in 1:4){
#     data = precip.ts.df[[r*2-1]][,-c(1:3)]
#     data[data > 1500] = NA
#     loc.ind <- which.max(colMeans(data,na.rm = T))
#     ### prepare the data frame to predict ###
#     data <- rbind(data.mean.obs.hist.day[[r]][[s]],data.mean.tep.245.day[[r]][[s]],
# 			data.mean.tep.585.day[[r]][[s]])
#     data$type = c(rep("Obs",nrow(data.mean.obs.hist.day[[r]][[s]])),
#                   rep("SSP 2-4.5",nrow(data.mean.tep.245.day[[r]][[s]])),
#                   rep("SSP 5-8.5",nrow(data.mean.tep.585.day[[r]][[s]])))
#     days <- data$date
#     d.ind = yday(days)
#     y.ind = year(days);y.ind = y.ind - y.ind[1] + 1
#     alt <- station.df[[r*2-1]]$elev[loc.ind]/1000
#     lon <- station.df[[r*2-1]]$Y[loc.ind]
#     lat <- station.df[[r*2-1]]$X[loc.ind]
#     main = paste0(lab.regions[r],"; ",model.selected[s],"; ","Station ",loc.ind," (", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
#     D = length(lat);Dt = length(days)
#     data.pred <- data.frame(temp = rep(data$tep,times=D),
#                             day = rep(d.ind,times=D),
#                             year = rep(y.ind,times=D),
#                             alt = rep(alt,each=Dt),
#                             lon = rep(lon,each=Dt),
#                             lat = rep(lat,each=Dt),
#                             col=rep(1:D,each=Dt),
#                             row=rep(1:Dt,times=D))
#     ind = !is.na(data.pred$temp);data.pred <- data.pred[ind,]
#     mean.pred  <- exp(predict.gam(results.gam,newdata = data.pred))
#     sig2.pred <- results.gam$sig2
#     shape.pred = 1/sig2.pred;scale.pred <- mean.pred/shape.pred;
#     quantile.pred <- qgamma(0.9,shape = shape.pred,scale=scale.pred)
#     data.pred$est.quantile = quantile.pred
#     prob.exceed.pred <- predict.gam(results.bin,newdata=data.pred,type="response")
#     gpd.pred <- predict(results.gpd,newdata=data.pred)
#     scale.gpd.pred = exp(gpd.pred[,1]);shape.gpd.pred = gpd.pred[1,2]
#     return.level = (1-1/(100*365)) 
#     prob.gpd <- (return.level - (1-prob.exceed.pred))/prob.exceed.pred
    
#     return.value <- y.thres + quantile.pred + qgpd(prob.gpd,loc=0,scale=scale.gpd.pred,shape=shape.gpd.pred)
#     by.list = list(season=getSeason(days)[ind],year=sapply(days,getYear)[ind],type=data$type[ind])
#     data <-aggregate(return.value,by=by.list,FUN=mean)
# 	data$season <- as.factor(data$season)
# 	data$type <- as.factor(data$type)
# 	data$x[data$type == "Obs" & data$season == "Winter" & data$year == 2016] = NA
#     p <- ggplot(data,aes(x=year,y=x,color=season,linetype=type)) + geom_line(alpha=1,size=0.6) 
#     p <- p + scale_linetype_manual(values=c("solid","dashed","dotted"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
#     p <- p +  xlab("Year") + ylab ("Return level (mm)") 
#     p <- p + labs(color="Season",linetype="Group") 
#     #p <- p + scale_color_manual(values=hcl.colors(4, "Berlin")) 
#     p <- p + ggtitle(main)
#     p <- p + theme(axis.text = element_text(size=10), 
#                    axis.title.x = element_text(size=14), 
#                   axis.title.y = element_text(size=14),
#                    plot.title = element_text(hjust = 0.5),
#                    panel.border = element_rect(fill = "transparent", # Needed to add the border
#                                                color = "black",            # Color of the border
#                                                size = 1),
#                    panel.background = element_rect(fill = "transparent")) 
    
#     #if(count %% 4 != 0){p <- p + theme(legend.position="none")}
#     p.list[[count]] <- p;count = count + 1
#     print(count)
#     }
#     }

#   pdf("figures/return_level_margins.pdf",width = 6,height = 3,onefile = TRUE)
#   for(i in 1:8){
#     show(p.list[i])
#   }
#   dev.off()