rm(list=ls())
library(ggplot2)
library(lubridate)
library(mgcv)
library(evgam)
library(gridExtra)
library(ggpubr)
library(evd)
load("data/dep.fit.boot.results.RData")
load("data/temperature.RData")
load("data/temperature_pred.RData")
load("data/precip.RData")
load("data/era5_geoinfo.RData")
load("data/transformed_coordinates.RData")
source("code/utility.R")

## Present the results ##
# 4 seasons * 8 regions * 2 risk functions
# create a data frame with the results 
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
idx.grid <- as.matrix(expand.grid(1:8, 1:4, 1:2))
boot.result.df <- data.frame(
    region = idx.grid[,1],
    season = idx.grid[,2],
    risk = idx.grid[,3],
    shape = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[1]}),
    lambda0 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[2]}),
    lambda1 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[3]}),
    sd.shape = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[1]}),
    sd.lambda0 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[2]}),
    sd.lambda1 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$sd[3]}))

# plot the results
# Assuming boot.result.df has columns for x, y, ymin, and ymax
p1 <- ggplot(boot.result.df, aes(x=factor(region), y=shape, color=factor(risk), group=season)) +
geom_point(size=1) +
geom_errorbar(aes(ymin=shape - 1.96*sd.shape, ymax=shape + 1.96*sd.shape), width=0.5) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(nu))  + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p2 <- ggplot(boot.result.df, aes(x=factor(region), y=lambda0, color=factor(risk), group=season)) +
geom_point(size=1) +
geom_errorbar(aes(ymin=lambda0 - 1.96*sd.lambda0, ymax=lambda0 + 1.96*sd.lambda0), width=0.5) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[0])) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p3 <- ggplot(boot.result.df, aes(x=factor(region), y=lambda1, color=factor(risk), group=season)) +
geom_point(size=1) +
geom_errorbar(aes(ymin=lambda1 - 1.96*sd.lambda1, ymax=lambda1 + 1.96*sd.lambda1), width=0.5,position=position_dodge(width=0)) +
ggh4x::facet_grid2(~season,scales="free_y",independent="y", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[1]))  + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

p3

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

pdf(file="figures/boot_plot.pdf",width=4*2+1,height=3)
grid.arrange(p1, p2, legend, nrow = 1, widths = c(4, 4, 1))
p3
dev.off()

save(boot.result.df,p1,p2,p3,file="data/boot_plot.RData")

## the prediciton into the future: marginal return level and the dependence range ## 

## plot the temperature covariate ##
model.selected = c(1,3,4)
align <- function(m,r,data){
    intercep.1 = as.Date("2015-01-01") <= date.df[,1] & date.df[,1] <= as.Date("2020-12-31") 
    intercep.2 = as.Date("2015-01-01") <= date.245 & date.245 <= as.Date("2020-12-31") 
    data.1 = data.frame(season=date.df[intercep.1,3],temp=temperature.covariate[[r]][intercep.1])
    data.2 = data.frame(season=getSeason(date.245[intercep.2]),temp=data[[m]][[r]][intercep.2])
    data_avg.1 = aggregate(temp ~ season, data.1, mean)
    data_avg.2 = aggregate(temp ~ season, data.2, mean)
    data_avg = data_avg.1$temp - data_avg.2$temp
    data.out = data.frame(season=getSeason(date.245),temp=data[[m]][[r]])
    data.out$temp = data.out$temp + data_avg[sapply(data.out$season,function(x){which(data_avg.1$season==x)})]
    return(data.out$temp)
}

for(m in model.selected){
    for(r in 1:8){
        temperature.245.avg[[m]][[r]] = align(m,r,temperature.245.avg)
        temperature.585.avg[[m]][[r]] = align(m,r,temperature.585.avg)
    }
}

count = 1;p.list <- list()
for(r in 1:8){
        data.1 = temperature.covariate[[r]]
        data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        
        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[,1],date.245,date.585)
        data.type = c(rep("Obs",length(data.1)),
                    rep("SSP 2-4.5",length(data.2)),
                    rep("SSP 5-8.5",length(data.3)))
        data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))

        data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
        p <- ggplot(data.df.avg, aes(x=year, y=tep, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5)
        p <- p  + xlab(NULL) + ylab(NULL) #xlab("Year") + ylab ("Temperature (°C)")
        p <- p + labs(color='Season',linetype='Group') 
        p <- p + scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) 
        p <- p + ggtitle(paste0(region.name[r]))
        p <- p + theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=16,face="bold",hjust=0.5),
                        axis.ticks =  element_line(linewidth = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) 
        p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),,keywidth = unit(1.5,"cm")))
        p.list[[count]] <- p;count = count + 1
}
  
pdf("figures/temperature_covariate.pdf",width = 20,height = 9)
ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()

save(p.list,file = "data/plot_temperature_covariate.RData")

## plot the marginal return level ##
model.selected = c(1,3,4)
y.thres=10
count = 1
p.list <- list()
for(r in 1:8){
    load(paste0("data/marginal_fit_301_",r,".RData"),e<-new.env())
    ### prepare the data frame to predict ###
    # print(r)
    # print(summary(e$results.gam))
    data.1 = temperature.covariate[[r]]
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[,1],date.245,date.585)
    data.type = c(rep("Obs",length(data.1)),
                rep("SSP 2-4.5",length(data.2)),
                rep("SSP 5-8.5",length(data.3)))
    data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))
    data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
    
    days <- c(315,135,225,45) #(Fall,Spring,Summer,Winter)
    data.df.avg$day = days[as.factor(data.df.avg$season)]
    idx.loc1 = which.max(unlist(lapply(precip[[r]],function(x){sum(na.omit(x)>0)})))
    idx.loc2 = station$group.id == region.id[r]
    alt <- station$elev[idx.loc2][idx.loc1]/1000
    lon <- station$Y[idx.loc2][idx.loc1]
    lat <- station$X[idx.loc2][idx.loc1]

    #main = paste0(region.name[r],"; ","Station ",idx.loc1," (", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
    main = paste0(region.name[r],"; ", "(", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
    D = length(lat);Dt = length(data.df.avg$day)
    data.pred <- data.frame(temp = rep(data.df.avg$tep,times=D),
                            day = rep(data.df.avg$day,times=D),
                            year = rep(data.df.avg$year,times=D),
                            alt = rep(alt,each=Dt),
                            lon = rep(lon,each=Dt),
                            lat = rep(lat,each=Dt))

    mean.pred  <- exp(predict.gam(e$results.gam,newdata = data.pred))
    sig2.pred <- e$results.gam$sig2
    shape.pred = 1/sig2.pred;scale.pred <- mean.pred/shape.pred;
    quantile.pred <- qgamma(0.9,shape = shape.pred,scale=scale.pred)
    data.pred$est.quantile = quantile.pred
    prob.exceed.pred <- predict.gam(e$results.bin,newdata=data.pred,type="response")
    gpd.pred <- predict(e$results.gpd,newdata=data.pred)
    scale.gpd.pred = exp(gpd.pred[,1]);shape.gpd.pred = gpd.pred[1,2]
    return.level = (1-1/(100*365)) 
    prob.gpd <- (return.level - (1-prob.exceed.pred))/prob.exceed.pred

    return.value <- y.thres + quantile.pred + qgpd(prob.gpd,loc=0,scale=scale.gpd.pred,shape=shape.gpd.pred)
    data.df.avg$return.value = return.value

    p <- ggplot(data.df.avg, aes(x=year, y=return.value, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.8,linewidth=1.5)
    p <- p + scale_linetype_manual(values=c("dotted","dashed","solid"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
    p <- p + xlab(NULL) + ylab(NULL) # xlab("Year") + ylab ("Return level (mm)") 
    p <- p + labs(color="Season",linetype="Group") 
    #p <- p + scale_color_manual(values=hcl.colors(4, "Berlin")) 
    p <- p + ggtitle(main)
    p <- p + theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=16,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) 
    p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),,keywidth = unit(1.5,"cm")))#if(count %% 4 != 0){p <- p + theme(legend.position="none")}
    p.list[[count]] <- p;count = count + 1
    print(count)
}

pdf("figures/return_level_margins.pdf",width = 20,height = 9,onefile = TRUE)
ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()
save(p.list,file = "data/plot_return_level_margins.RData")


## plot qqplot for random locations ##
for(idx in 1:8){
    load(paste0("data/marginal_fit_301_",idx,".RData"),e<-new.env())
    # sig2.pred <- e$results.gpd$sig2
    # shape.pred = 1/sig2.pred
    shape.pred = fitted(e$results.gpd)[1,2]
    print(shape.pred)
#    }
    set.seed(1000)
    idx.list = sample(1:sum(station$group.id==region.id[idx]),2,replace = F,prob=apply(e$U,2,function(x){sum(!is.na(x))}))
    png(file = paste0("figures/qqplot_marginal_",idx,".png"),height=6,width=6*3,units="cm",res=300, pointsize=6)
    par(mfrow=c(1,3),mar=c(3,4,3,1),mgp=c(2.5,2,0),cex.lab=3,cex.axis=3,cex.main=3)
    theoretical.quantiles <- qgpd(e$U[,idx.list],loc=1,scale=1,shape=0)
    empirical.quantiles <- qgpd(1:1000/(1+1000),loc=1,scale=1,shape=0)
    theoretical.quantiles <- split(theoretical.quantiles,col(theoretical.quantiles))
    theoretical.quantiles <- sapply(theoretical.quantiles,function(x) quantile(x,prob=1:1000/(1+1000),na.rm=T),simplify = F)
    qqplot(unlist(theoretical.quantiles),empirical.quantiles,cex=1.5,pch=20,
        xlab="",ylab="",main="All stations")
        #xlab="Theoretical quantiles",ylab="Empirical quantiles",main="All stations")
    abline(0,1,col=2,lwd=2)
    for(ind in 1:length(idx.list)){
    qqplot(theoretical.quantiles[[ind]],empirical.quantiles,
        cex=1.5,pch=20,
        xlab="",ylab="",main=paste("Station",idx.list[ind]))
    abline(0,1,col=2,lwd=2)
    }        
    dev.off() 
}

## plot the tail-correlation range ##
model.selected <- c(1,3,4)
count = 1;p.list1 <- list()
p.list2 <- list()
for(r in 1:8){
    data.1 = temperature.covariate[[r]]
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[,1],date.245,date.585)
    data.type = c(rep("Obs",length(data.1)),
                rep("SSP 2-4.5",length(data.2)),
                rep("SSP 5-8.5",length(data.3)))
    data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))
    data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
    risk = rep(1:2,nrow(data.df.avg));data.df.avg <- rbind(data.df.avg,data.df.avg)
    data.df.avg$risk = risk
    data.df.avg$season = as.numeric(factor(data.df.avg$season,season))
    data.df.avg = merge(data.df.avg,subset(boot.result.df,region==r),by=c("season","risk"))
    data.df.avg$season = season[data.df.avg$season]
    data.df.avg$range = sapply(1:nrow(data.df.avg),function(i){solve.h.BR(c(data.df.avg$shape[i],data.df.avg$lambda0[i],data.df.avg$lambda1[i]),temp=data.df.avg$tep[i],logval=TRUE)})

    p.list1[[r]] <- ggplot(subset(data.df.avg,risk==1), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5) + 
        xlab(NULL) + ylab (NULL) + labs(color='Season',linetype='Group') + 
        scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
        ggtitle(paste0(region.name[r]," with risk functional 1")) +
        theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=16,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) + 
        guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),,keywidth = unit(1.5,"cm")))

    
    p.list2[[r]] <- ggplot(subset(data.df.avg,risk==2), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5)  + 
        xlab(NULL) + ylab (NULL) +
        labs(color='Season',linetype='Group')  + 
        scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
        ggtitle(paste0(region.name[r]," with risk functional 2")) +
        theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=16,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) + 
        guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),,keywidth = unit(1.5,"cm")))    
}

pdf("figures/tail_correlation_range.pdf",width = 20,height = 9,onefile = TRUE)
ggarrange(plotlist=p.list1,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
ggarrange(plotlist=p.list2,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()

save(p.list1,p.list2,file = "data/tail_correlation_range.RData")

## Simulations ##
library(mvPotST)
loc = as.matrix(expand.grid((1:50)/5,(1:50)/5),ncol=2)
param_mat <- as.matrix(expand.grid(1,log(c(2,5,10)),0))
dep_range <- apply(param_mat,1,solve.h.BR,temp=0,logval=FALSE)
simu <- list()
for(i in 1:nrow(param_mat)){
	param = param_mat[i,-3]
	vario <- function(h,par=param){ 
  		## reparametrization
  		alpha = par[1]
  		lambda <- exp(par[2])
  		val=(sqrt(sum(h^2))/lambda)^alpha
  	    return(val)
	}
	set.seed(245254325)
    simu[[i]] <- unlist(simulPareto(1,loc=loc,vario=vario,nCores=4,u=1,shape=1/5))
}

pic <- list()
#MaxLimits = evd::qgpd(0.99,loc=1,scale=1,shape=1/5)
for(i in 1:length(simu)){
 MaxLimits = max(simu[[i]])
 data <- data.frame(x=loc[,1],y=loc[,2],z=pmin(simu[[i]],MaxLimits))
 pic1 <- ggplot(aes(x=x,y=y,fill=z),data=data) + geom_tile()
 pic1 <- pic1 + scale_fill_gradientn(name="Value",colors=hcl.colors(12,"BluYl",rev=TRUE,alpha=1),limits=c(0,MaxLimits))
 pic1 <- pic1 + theme(axis.text = element_text(size=10), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          plot.title = element_text(hjust = 0.5),aspect.ratio=1) + coord_fixed() 
 pic[[i]] <- pic1 + ggtitle(
  label=bquote(paste("Tail-correlation range:", ~.(round(dep_range[i],1)),";",~lambda==.(exp(param_mat[i,2]))))
  )
}


pdf("figures/simulation.pdf",width=4.5*3,height=4,onefile=TRUE)
combined_pic <- plot_grid(pic[[1]],pic[[2]],pic[[3]],
nrow=1,rel_widths=c(1,1,1),label_y = expression(alpha==0.5),hjust=-1,align="vh")
show(combined_pic)
dev.off()

save(pic,param_mat,vario,loc,simu,file="data/plot_simulation.RData")

## plot the location of 
library(sf)
load("data/temperature_pred.RData")
load("data/temperatures-mississippi.RData")
load("data/era5_geoinfo.RData")
g<-list()
colors = c(brewer.pal(6,"Set2"),brewer.pal(12,"Paired"))
for(idx in 1:5){
    i = which(models==idx.models[idx] & periods=="historical") ## 
    idx.group.id = locate.idx.list[[idx]][locate.idx.list[[idx]][,2]!=19,1]
    data = data.frame(x = xyt[[i]]$lon[idx.grid.list[[idx]]$miss[,1]],y=xyt[[i]]$lat[idx.grid.list[[idx]]$miss[,2]])[idx.group.id,]
    data$group.id = locate.idx.list[[idx]][locate.idx.list[[idx]][,2]!=19,2]
    g[[idx]] <- ggplot() + geom_sf(data=shape1,aes(fill=names),alpha=0.2) + geom_point(data=data,aes(x=x,y=y,col=as.factor(group.id)),size=0.5) + coord_sf(xlim = c(-79, -114), ylim = c(29, 49)) + scale_fill_manual(values=colors) + scale_color_brewer(palette="Dark2") + ggtitle(idx.models[idx]) + theme(plot.title = element_text(hjust = 0.5)) + labs(col="Regions",fill="Region names",x="Longitude",y="Latitude")
    #g[[idx]] <- g[[idx]] + guides(fill=FALSE,colour=FALSE)
}

png("figures/temperature-mississippi_location.png",width=10*5,height=10,res=300,units="cm")
ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],ncol=5,common.legend=TRUE,legend="bottom")
dev.off()

load("data/temperatures-danube.RData")
g1<-list()
for(idx in 1:5){
    i = which(models==idx.models[idx] & periods=="historical") ## 
    idx.group.id = locate.idx.list[[idx]][locate.idx.list[[idx]][,2]==19,1]
    data = data.frame(x = xyt[[i]]$lon[idx.grid.list[[idx]]$danube[,1]],y=xyt[[i]]$lat[idx.grid.list[[idx]]$danube[,2]],group.id=19)[idx.group.id,]
    g1[[idx]] <- ggplot() + geom_sf(data=shape3,aes(fill="blue"),alpha=0.2) + geom_point(data=data,aes(x=x,y=y),color="red",size=0.5) + ggtitle(idx.models[idx]) + theme(plot.title = element_text(hjust = 0.5))
    g1[[idx]] <- g1[[idx]] + guides(fill=FALSE,colour=FALSE) + labs(x="Longitude",y="Latitude")
}

png("figures/temperature-danube_location.png",width=10*5,height=10,res=300,units="cm")
ggarrange(g1[[1]],g1[[2]],g1[[3]],g1[[4]],g1[[5]],ncol=5,common.legend=TRUE,legend="bottom")
dev.off()

save(g,g1,file="data/plot_location_climate_models.RData")

##plot the location of the stations ##
load("data/era5_geoinfo.RData")
colors =c(brewer.pal(7,"Paired"),"grey50")
shape1$names2 = shape1$names;shape1$names2[-region.id[-1]] = "Others";shape1$names2 = factor(shape1$names2,levels=c(region.name[-1],"Others")) 
data = data.frame(x = station$Y,y=station$X,region=station$group.id)[station$group.id!=19,]
p1 <- ggplot() + geom_sf(data=shape1,aes(fill=names2),alpha=0.5) + geom_point(data=data,aes(x=x,y=y,col=as.factor(region)),size=0.5) + coord_sf(xlim = c(-79, -114), ylim = c(29, 49)) + scale_fill_manual(values=colors) + scale_color_brewer(palette="Dark2") + labs(fill="Subregions",x="Longitude",y="Latitude") + guides(colour=FALSE) + theme(legend.position = "right",plot.title = element_text(hjust = 0.5)) + ggtitle("Mississippi")


data = data.frame(x = station$Y,y=station$X,region=station$group.id)[station$group.id==19,]
p2 <- ggplot() + geom_sf(data=shape3,aes(fill=colors[1]),alpha=0.5) + geom_point(data=data,aes(x=x,y=y),size=0.5) + scale_fill_manual(values=colors) + scale_color_brewer(palette="Dark2") + labs(x="Longitude",y="Latitude") + guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Danube")

pdf("figures/precip_loc.pdf",width=11,height=3.5,onefile=TRUE)
grid.arrange(p1,p2,nrow=1,widths=c(6,5))
dev.off()

