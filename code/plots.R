rm(list=ls())
library(ggplot2)
library(gridExtra)
library(cowplot)
library(lubridate)
library(mgcv)
library(evgam)
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
    lambda0 = apply(idx.grid,1,function(x){boot.result.list[[x[3]]][[x[2]]][[x[1]]]$true[2]}+log(2)),
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
    data.1 = data.frame(season=date.df[intercep.1,3],temp=temperature.covariate[[r]][intercep.1]/10)
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
        data.1 = temperature.covariate[[r]]/10
        data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        
        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[,1],date.245,date.585)
        data.type = c(rep("Obs",length(data.1)),
                    rep("SSP 2-4.5",length(data.2)),
                    rep("SSP 5-8.5",length(data.3)))
        data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))

        data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
        p <- ggplot(data.df.avg, aes(x=year, y=tep, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9)
        p <- p  + xlab("Year") + ylab ("Temperature (°C)")
        p <- p + labs(color='Season',linetype='Group') 
        p <- p + scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) 
        p <- p + ggtitle(paste0(region.name[r],"; ","Average"))
        p <- p + theme(axis.text = element_text(size=10), 
                        axis.title.x = element_text(size=14), 
                        axis.title.y = element_text(size=14),
                        plot.title = element_text(hjust = 0.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent"))
        p.list[[count]] <- p;count = count + 1
}
  
pdf("figures/temperature_covariate.pdf",width = 6,height = 3,onefile = TRUE)
for(i in 1:8){
    show(p.list[i])
}
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
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)*10
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)*10
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[,1],date.245,date.585)
    data.type = c(rep("Obs",length(data.1)),
                rep("SSP 2-4.5",length(data.2)),
                rep("SSP 5-8.5",length(data.3)))
    data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))
    data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
    
    days <- c(315,135,225,45) #(Fall,Spring,Summer,Winter)
    data.df.avg$day = days[as.factor(data.df.avg$season)]
    idx.loc1 = which.max(unlist(lapply(precip[[r]],mean,na.rm=TRUE)))
    idx.loc2 = station$group.id == region.id[r]
    alt <- station$elev[idx.loc2][idx.loc1]/1000
    lon <- station$Y[idx.loc2][idx.loc1]
    lat <- station$X[idx.loc2][idx.loc1]

    main = paste0(region.name[r],"; ","Station ",idx.loc1," (", round(lat,2) ,"N°,",round(lon,2),"E°" ,") ")
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
    quantile.pred <- qgamma(0.8,shape = shape.pred,scale=scale.pred)
    data.pred$est.quantile = quantile.pred
    prob.exceed.pred <- predict.gam(e$results.bin,newdata=data.pred,type="response")
    gpd.pred <- predict(e$results.gpd,newdata=data.pred)
    scale.gpd.pred = exp(gpd.pred[,1]);shape.gpd.pred = gpd.pred[1,2]
    return.level = (1-1/(100*365)) 
    prob.gpd <- (return.level - (1-prob.exceed.pred))/prob.exceed.pred

    return.value <- y.thres + quantile.pred + qgpd(prob.gpd,loc=0,scale=scale.gpd.pred,shape=shape.gpd.pred)
    data.df.avg$return.value = return.value

    p <- ggplot(data.df.avg, aes(x=year, y=return.value, group=interaction(season,type), color=season, linetype=type)) + geom_line()
    p <- p + scale_linetype_manual(values=c("dotted","dashed","solid"),labels=c("Obs","SSP 2-4.5","SSP 5-8.5"))
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

pdf("figures/return_level_margins.pdf",width = 6,height = 3,onefile = TRUE)
for(i in 1:8){
    show(p.list[i])
}
dev.off()
save(p.list,file = "data/plot_return_level_margins.RData")

## plot the tail-correlation range ##
model.selected <- c(1,3,4)
count = 1;p.list <- list()
for(r in 1:8){
    data.1 = temperature.covariate[[r]]/10
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
    data.df.avg$range = sapply(1:nrow(data.df.avg),function(i){solve.h.BR(c(data.df.avg$shape[i],data.df.avg$lambda0[i],data.df.avg$lambda1[i]),temp=data.df.avg$tep[i]*10,logval=TRUE)})

    p[[1]] <- ggplot(subset(data.df.avg,risk==1), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9) + 
        xlab("Year") + ylab ("Logarithmic range") + labs(color='Season',linetype='Group') + 
        scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
        ggtitle(paste0(region.name[r]," with risk functional 1")) +
        theme(axis.text = element_text(size=10), 
                        axis.title.x = element_text(size=14), 
                        axis.title.y = element_text(size=14),
                        plot.title = element_text(hjust = 0.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent"))
    
    p[[2]] <- ggplot(subset(data.df.avg,risk==2), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9)  + 
        xlab("Year") + ylab ("Logarithmic range") +
        labs(color='Season',linetype='Group')  + 
        scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
        ggtitle(paste0(region.name[r]," with risk functional 2")) +
        theme(axis.text = element_text(size=10), 
                        axis.title.x = element_text(size=14), 
                        axis.title.y = element_text(size=14),
                        plot.title = element_text(hjust = 0.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent"))
    
    legend <- get_legend(p[[1]])
    p[[1]] <- p[[1]] + theme(legend.position = "none")
    p[[2]] <- p[[2]] + theme(legend.position = "none")
    p[[3]] <- legend
    p.list[[r]] <- p
}

pdf("figures/tail_correlation_range.pdf",width = 4*2+1,height = 3,onefile = TRUE)
for(i in 1:length(p.list)){
    grid.arrange(p.list[[i]][[1]], p.list[[i]][[2]], p.list[[i]][[3]], nrow = 1, widths = c(4, 4, 1))
}
dev.off()

save(p.list,file = "data/tail_correlation_range.RData")

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
  label=bquote(paste("Tail-correlation range:", ~.(round(dep_range[i],1)),";",~lambda==.(param_mat[i,1])))
  )
}


pdf("figures/simulation.pdf",width=4.5*3,height=4,onefile=TRUE)
combined_pic <- plot_grid(pic[[1]],pic[[2]],pic[[3]],
nrow=1,rel_widths=c(1,1,1),label_y = expression(alpha==0.5),hjust=-1,align="vh")
show(combined_pic)
dev.off()

save(pic,param_mat,vario,loc,simu,file="data/plot_simulation.RData")
