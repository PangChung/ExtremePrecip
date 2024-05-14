rm(list=ls())
library(ggplot2)
library(lubridate)
library(mgcv)
library(evgam)
library(gridExtra)
library(ggpubr)
library(evd)
DataPath <- "data/boot4/"
load("data/dep.fit.boot.results3.RData")
load("data/temperature.RData")
load("data/temperature_pred.RData")
load("data/precip.RData")
load("data/era5_geoinfo.RData")
load("data/transformed_coordinates.RData")
source("code/utility.R")
season = c("Winter" ,"Spring" ,"Summer" ,"Fall")
## Present the results ##
# 4 seasons * 8 regions * 2 risk functions
# create a data frame with the results 
# plot the results
# Assuming boot.result.df has columns for x, y, ymin, and ymax
p1 <- p2 <- list()
p1[[1]] <- ggplot(boot.result.df, aes(x=factor(region), y=shape, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=shape - 1.96*sd.shape, ymax=shape + 1.96*sd.shape), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(nu))

p1[[2]] <- ggplot(boot.result.df, aes(x=factor(region), y=lambda0, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=lambda0 - 1.96*sd.lambda0, ymax=lambda0 + 1.96*sd.lambda0), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[0]))

p1[[3]] <- ggplot(boot.result.df, aes(x=factor(region), y=lambda1, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=lambda1 - 1.96*sd.lambda1, ymax=lambda1 + 1.96*sd.lambda1), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season,scales="free_y",independent="y", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[1])) + geom_hline(yintercept = 0, linetype="dashed", color = "black")

p1[[1]] <- p1[[1]] + theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")
p1[[2]] <- p1[[2]]  + theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")
p1[[3]] <- p1[[3]]  +  theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")

pdf(file="figures/boot_plot.pdf",width=16,height=4)
p1[[1]]
p1[[2]]
p1[[3]]
dev.off()

p2[[1]] <- ggplot(boot.result.df, aes(x=factor(region), y=shape, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=low.shape, ymax=high.shape), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(nu)) 
p2[[2]] <- ggplot(boot.result.df, aes(x=factor(region), y=lambda0, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=low.lambda0, ymax=high.lambda0), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season, space = "free_x", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[0])) 

p2[[3]] <- ggplot(boot.result.df, aes(x=factor(region), y=lambda1, color=factor(risk), group=season)) +
geom_point(size=1,position=position_dodge2(width=0.5)) +
geom_errorbar(aes(ymin=low.lambda1, ymax=high.lambda1), width=0.5,linewidth=1.5,position=position_dodge2(width=1)) +
ggh4x::facet_grid2(~season,scales="free_y",independent="y", labeller = labeller(season = as_labeller(c("1" = "Winter", "2" = "Spring", "3" = "Summer", "4" = "Fall")))) + labs(color="risk functional",x="Region",y=expression(lambda[1])) + geom_hline(yintercept = 0, linetype="dashed", color = "black")

p2[[1]] <- p2[[1]] +   theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")
p2[[2]] <- p2[[2]] +   theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")
p2[[3]] <- p2[[3]] +   theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.position = "none")

pdf(file="figures/boot_plot2.pdf",width=16,height=4)
p2[[1]]
p2[[2]]
p2[[3]]
dev.off()

est.mat <- NULL
col.names <- c("region","season","risk","nu","lambda[0]","lambda[1]")
for(i in 1:nrow(boot.result.df)){
    data.temp = boot.result.list[[boot.result.df$risk[i]]][[boot.result.df$season[i]]][[boot.result.df$region[i]]]$jack
    data.id = matrix(as.numeric(boot.result.df[i,1:3]),ncol=3,nrow=nrow(data.temp),byrow=TRUE)
    est.mat = rbind(est.mat,cbind(data.id,data.temp))
}
est.df <- as.data.frame(matrix(est.mat,ncol=6))
names(est.df) <- col.names
est.df$season = factor(est.df$season,labels=c("Winter" ,"Spring" ,"Summer" ,"Fall"))
est.df$risk = factor(est.df$risk,labels=c("Risk functional 1","Risk functional 2"))
colnames(boot.result.df)[1:6] <- col.names
boot.result.df$season = factor(boot.result.df$season,labels=c("Winter" ,"Spring" ,"Summer" ,"Fall"))
boot.result.df$risk = factor(boot.result.df$risk,labels=c("Risk functional 1","Risk functional 2"))
p.list <- list()

p.list[[1]] <- ggplot(est.df, aes(x = factor(region), y = nu, fill = risk)) +
  geom_violin(position = position_dodge(width=0.75),draw_quantiles = c(0.025,0.975)) + 
  facet_wrap(~ season,scales = "free",nrow=2,ncol=4,labeller = label_parsed) +
  geom_point(data=boot.result.df[,1:6],color="black",size=1,position=position_dodge(width = 0.75)) +
  labs(title = "",
       x = "Region",
       y = bquote(nu),
       fill = "Risk") + 
  theme(axis.text = element_text(size = 16,face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        strip.text = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"),
        legend.title = element_text(size = 16),
        legend.position = "none")

p.list[[2]] <- ggplot(est.df, aes(x = factor(region), y = `lambda[0]`, fill = risk)) +
  geom_violin(position = position_dodge(width=0.75),draw_quantiles = c(0.025,0.975)) + 
  facet_wrap(~ season,scales = "free",nrow=2,ncol=4,labeller = label_parsed) +
  geom_point(data=boot.result.df[,1:6],color="black",size=1,position=position_dodge(width = 0.75)) +
  labs(title = "",
       x = "Region",
       y = bquote(lambda[0]),
       fill = "Risk") + 
  theme(axis.text = element_text(size = 16,face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        strip.text = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"),
        legend.title = element_text(size = 16),
        legend.position = "none")

p.list[[3]] <- ggplot(est.df, aes(x = factor(region), y = `lambda[1]`, fill = risk)) +
  geom_violin(position = position_dodge(width=0.75),draw_quantiles = c(0.05,0.95)) + 
  facet_wrap(~ season,scales = "free",nrow=2,ncol=4,labeller = label_parsed) +
  geom_point(data=boot.result.df[,1:6],color="black",size=1,position=position_dodge(width = 0.75)) + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(title = "",
       x = "Region",
       y = bquote(lambda[1]),
       fill = "Risk") +
  theme(axis.text = element_text(size = 16,face="bold"),
        axis.title.x = element_text(size = 16,face="bold"),
        strip.text = element_text(size = 16,face="bold"),
        axis.title.y = element_text(size = 16,face="bold"),
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"),
        legend.title = element_text(size = 16),
        legend.position = "none")

pdf(file="figures/boot_plot3.pdf",width=4*4,height=4,onefile = TRUE)
p.list[[1]]
p.list[[2]]
p.list[[3]]
dev.off()

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
        idx.obs = date.df$date>as.Date(START.date)
        data.1 = temperature.covariate[[r]][idx.obs]
        data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        
        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[idx.obs,1],date.245,date.585)
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
                        legend.title = element_text(size = 16, face = "bold"),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) 
        p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16, face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))
        p.list[[count]] <- p;count = count + 1
}
  
pdf("figures/temperature_covariate.pdf",width = 16,height = 8)
ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()
save(p.list,file = paste0(DataPath,"plot_temperature_covariate_.RData"))

## for individual climate outputs ##
pdf(paste0("figures/temperature_covariate_i.pdf"),width = 16,height = 8,onefile=TRUE)
for(i in model.selected){
    count = 1;p.list <- list()
    for(r in 1:8){
        idx.obs = date.df$date>as.Date(START.date)
        data.1 = temperature.covariate[[r]][idx.obs]
        data.2 = temperature.245.avg[[i]][[r]] #apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = temperature.585.avg[[i]][[r]] #apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)

        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[idx.obs,1],date.245,date.585)
        data.type = c(rep("Obs",length(data.1)),
            rep("SSP 2-4.5",length(data.2)),
            rep("SSP 5-8.5",length(data.3)))
        data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))

        data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
        p <- ggplot(data.df.avg, aes(x=year, y=tep, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5)
        p <- p  + xlab(NULL) + ylab(NULL) #xlab("Year") + ylab ("Temperature (°C)")
        p <- p + labs(color='Season',linetype='Group') 
        p <- p + scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) 
        p <- p + ggtitle(paste0(region.name[r],"; ",idx.models[i]))
        p <- p + theme(axis.text = element_text(size=16,face="bold"),
                plot.title = element_text(size=14,face="bold",hjust=0.5),
                axis.ticks =  element_line(linewidth = 1.5),
                panel.border = element_rect(fill = "transparent", # Needed to add the border
                                            color = "black",            # Color of the border
                                            linewidth = 1),
                panel.background = element_rect(fill = "transparent")) 
        p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),,keywidth = unit(1.5,"cm")))
        p.list[[count]] <- p;count = count + 1
    }
    show(ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom"))
}
dev.off()  

## plot the marginal return level ##
model.selected = c(1,3,4)
y.thres=0
count = 1
p.list <- list()
for(r in 1:8){
    load(paste0(DataPath,"marginal_fit_0_",r,".RData"),e<-new.env())
    ### prepare the data frame to predict ###
    # print(r)
    # print(summary(e$results.gam))
    idx.obs = date.df$date>as.Date(START.date)
    data.1 = temperature.covariate[[r]][idx.obs]
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[idx.obs,1],date.245,date.585)
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
                       plot.title = element_text(size=14,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) 
    p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))#if(count %% 4 != 0){p <- p + theme(legend.position="none")}
    p.list[[count]] <- p;count = count + 1
    print(count)
}

pdf("figures/return_level_margins.pdf",width = 16,height = 8,onefile = TRUE)
ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()

save(p.list,file = paste0(DataPath,"plot_return_level_margins.RData"))

## for individual climate outputs ##
model.selected = c(1,3,4)
pdf("figures/return_level_margins_i.pdf",width = 18,height = 8,onefile = TRUE)
for(i in model.selected){
y.thres=0
count = 1
p.list <- list()
    for(r in 1:8){
        load(paste0(DataPath,"marginal_fit_0_",r,".RData"),e<-new.env())
        ### prepare the data frame to predict ###
        # print(r)
        # print(summary(e$results.gam))
        idx.obs = date.df$date>as.Date(START.date)
        data.1 = temperature.covariate[[r]][idx.obs]
        data.2 = temperature.245.avg[[i]][[r]] #apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = temperature.585.avg[[i]][[r]] #apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[idx.obs,1],date.245,date.585)
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
        main = paste0(region.name[r],"; ", "(", round(lat,2) ,"N°,",round(lon,2),"E°" ,")","; ",idx.models[i])
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
                        plot.title = element_text(size=13,face="bold",hjust=0.5),
                            axis.ticks =  element_line(size = 1.5),
                            panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                        color = "black",            # Color of the border
                                                        linewidth = 1),
                            panel.background = element_rect(fill = "transparent")) 
        p <- p + guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))#if(count %% 4 != 0){p <- p + theme(legend.position="none")}
        p.list[[count]] <- p;count = count + 1
        print(count)
    }
    show(ggarrange(plotlist=p.list,nrow=2,ncol=4,common.legend=TRUE,legend="bottom"))
}
dev.off()

## plot qqplot for random locations ##
for(idx in 1:8){
    load(paste0(DataPath,"marginal_fit_0_",idx,".RData"),e<-new.env())
    # sig2.pred <- e$results.gpd$sig2
    # shape.pred = 1/sig2.pred
    shape.pred = fitted(e$results.gpd)[1,2]
    set.seed(1234435)
    idx.list = sample(1:sum(station$group.id==region.id[idx]),2,replace = F,prob=apply(e$U,2,function(x){sum(!is.na(x))}))
    print(idx.list)
    png(file = paste0("figures/qqplot_marginal_",idx,".png"),height=6,width=6*3,units="cm",res=300, pointsize=6)
    par(mfrow=c(1,3),mar=c(3,4,3,1),mgp=c(2.5,2,0),cex.lab=3,cex.axis=3,cex.main=3,pty="s")
    theoretical.quantiles <- qgpd(e$U,loc=0,scale=1,shape=shape.pred)
    empirical.quantiles <- qgpd(1:1000/(1+1000),loc=0,scale=1,shape=shape.pred)
    theoretical.quantiles <- split(theoretical.quantiles,col(theoretical.quantiles))
    theoretical.all <- quantile(unlist(theoretical.quantiles),prob=1:1000/(1+1000),na.rm=T)
    theoretical.quantiles <- sapply(theoretical.quantiles,function(x) quantile(x,prob=1:1000/(1+1000),na.rm=T),simplify = F)
    xlim = range(c(theoretical.all,empirical.quantiles)) + c(-0.1,0.1)
    plot(theoretical.all,empirical.quantiles,cex=1.5,pch=20,
        xlab="",ylab="",main="All stations",asp=1,xlim=xlim,ylim=xlim)
        #xlab="Theoretical quantiles",ylab="Empirical quantiles",main="All stations")
    abline(0,1,col=2,lwd=2)
    for(ind in 1:length(idx.list)){
    xlim = range(c(theoretical.quantiles[[idx.list[ind]]],empirical.quantiles)) + c(-0.1,0.1)
    plot(theoretical.quantiles[[idx.list[ind]]],empirical.quantiles,
        cex=1.5,pch=20,xlab="",ylab="",main=paste("Station",idx.list[ind]),asp=1,xlim=xlim,ylim=xlim)
    abline(0,1,col=2,lwd=2)
    }        
    dev.off() 
}

## plot the tail-correlation range ##
load("data/dep.fit.boot.results3.RData")
model.selected <- c(1,3,4)
count = 1;p.list1 <- list()
p.list2 <- list()
for(r in 1:8){
    idx.obs = date.df$date>as.Date(START.date)
    data.1 = temperature.covariate[[r]][idx.obs]
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[idx.obs,1],date.245,date.585)
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
        ggtitle(paste0(region.name[r],"; risk functional 1")) +
        theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=14,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) + 
        guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))

    
    p.list2[[r]] <- ggplot(subset(data.df.avg,risk==2), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5)  + 
        xlab(NULL) + ylab (NULL) +
        labs(color='Season',linetype='Group')  + 
        scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
        ggtitle(paste0(region.name[r],"; risk functional 2")) +
        theme(axis.text = element_text(size=16,face="bold"),
                       plot.title = element_text(size=14,face="bold",hjust=0.5),
                        axis.ticks =  element_line(size = 1.5),
                        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                    color = "black",            # Color of the border
                                                    linewidth = 1),
                        panel.background = element_rect(fill = "transparent")) + 
        guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))    
}

pdf("figures/tail_correlation_range.pdf",width = 18,height = 8,onefile = TRUE)
ggarrange(plotlist=p.list1,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
ggarrange(plotlist=p.list2,nrow=2,ncol=4,common.legend=TRUE,legend="bottom")
dev.off()

save(p.list1,p.list2,file = paste0(DataPath,"tail_correlation_range.RData"))

## for individual climate outputs ##
load("data/dep.fit.boot.results3.RData")
pdf("figures/tail_correlation_range_i.pdf",width = 18,height = 8,onefile = TRUE)
model.selected <- c(1,3,4)
for(i in model.selected){
    count = 1;p.list1 <- list()
    p.list2 <- list()
    for(r in 1:8){
        idx.obs = date.df$date>as.Date(START.date)
        data.1 = temperature.covariate[[r]][idx.obs]
        data.2 = temperature.245.avg[[i]][[r]]#apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.3 = temperature.585.avg[[i]][[r]]#apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
        data.temp = c(data.1, data.2, data.3)
        date.temp = c(date.df[idx.obs,1],date.245,date.585)
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
            ggtitle(paste0(region.name[r],"; risk functional 1","; ",idx.models[i])) +
            theme(axis.text = element_text(size=16,face="bold"),
                        plot.title = element_text(size=13,face="bold",hjust=0.5),
                            axis.ticks =  element_line(size = 1.5),
                            panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                        color = "black",            # Color of the border
                                                        linewidth = 1),
                            panel.background = element_rect(fill = "transparent")) + 
            guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))

        
        p.list2[[r]] <- ggplot(subset(data.df.avg,risk==2), aes(x=year, y=range, group=interaction(season,type), color=season, linetype=type)) + geom_line(alpha=0.9,linewidth=1.5)  + xlab(NULL) + ylab (NULL) +
            labs(color='Season',linetype='Group')  + 
            scale_color_manual(values=hcl.colors(4,"Dynamic")) + scale_linetype_manual(values=c("dotted","dashed","solid")) +
            ggtitle(paste0(region.name[r],"; risk functional 2","; " ,idx.models[i])) +
            theme(axis.text = element_text(size=16,face="bold"),
                        plot.title = element_text(size=13,face="bold",hjust=0.5),
                            axis.ticks =  element_line(size = 1.5),
                            panel.border = element_rect(fill = "transparent", # Needed to add the border
                                                        color = "black",            # Color of the border
                                                        linewidth = 1),
                            panel.background = element_rect(fill = "transparent")) + 
            guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))    
    }
    show(ggarrange(plotlist=p.list1,nrow=2,ncol=4,common.legend=TRUE,legend="bottom"))
    show(ggarrange(plotlist=p.list2,nrow=2,ncol=4,common.legend=TRUE,legend="bottom"))
}
dev.off()

## plot the tail correlation range for the bootstrap estimates to ##
## check the uncertainty of the estimates ##
solve.h.boot <- function(i,j,k,temp){
    est = boot.result.list[[i]][[j]][[k]]$jack
    return(quantile(apply(est,1,solve.h.BR,temp=temp,logval=TRUE),c(0.025,0.975)))
}

load("data/dep.fit.boot.results3.RData")
model.selected <- c(1,3,4)
count = 1;
p.list1 <- p.list2 <- list()
for(r in 1:8){
    idx.obs = date.df$date>as.Date(START.date)
    data.1 = temperature.covariate[[r]][idx.obs]
    data.2 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.245.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.3 = apply(matrix(unlist(lapply(model.selected,function(i){temperature.585.avg[[i]][[r]]})),ncol=length(model.selected),byrow=FALSE),1,mean)
    data.temp = c(data.1, data.2, data.3)
    date.temp = c(date.df[idx.obs,1],date.245,date.585)
    data.type = c(rep("Obs",length(data.1)),
                rep("SSP 2-4.5",length(data.2)),
                rep("SSP 5-8.5",length(data.3)))
    data.df = data.frame(season=getSeason(date.temp),tep=data.temp,type=data.type,year=getYear(date.temp))
    data.df.avg <- aggregate(tep ~ season + year + type, data.df, mean)
    risk = rep(1:2,nrow(data.df.avg));data.df.avg <- rbind(data.df.avg,data.df.avg)
    data.df.avg$risk = risk
    data.df.avg$season = as.numeric(factor(data.df.avg$season,season))
    data.df.avg = merge(data.df.avg,subset(boot.result.df,region==r),by=c("season","risk"))
    data.df.avg$range = sapply(1:nrow(data.df.avg),function(i){solve.h.BR(c(data.df.avg$shape[i],data.df.avg$lambda0[i],data.df.avg$lambda1[i]),temp=data.df.avg$tep[i],logval=TRUE)})
    data.df.avg[,c("low","up")] = t(sapply(1:nrow(data.df.avg),function(i){solve.h.boot(data.df.avg$risk[i],data.df.avg$season[i],r,temp=data.df.avg$tep[i])}))
    
    data.df.avg$season = season[data.df.avg$season]
    p.list1[[r]] <- ggplot(subset(data.df.avg,risk==1 & type=="SSP 5-8.5"), aes(x=year, y=range)) + geom_line(alpha=0.9,linewidth=1.5) + 
        xlab(NULL) + ylab (NULL) + geom_ribbon(aes(ymin=low,ymax=up),alpha=0.3,linewidth=1.5) + facet_wrap(~ season,scales = "free",nrow=1,ncol=4,labeller = label_parsed) +
        ggtitle(paste0(region.name[r],"; risk functional 1")) +
        theme(axis.text = element_text(size=16,face="bold"),
                strip.text = element_text(size = 16),
                plot.title = element_text(size=16,face="bold",hjust=0.5),
                axis.ticks =  element_line(size = 1.5),
                panel.border = element_rect(fill = "transparent", # Needed to add the border
                                        color = "black",            # Color of the border
                                        linewidth = 1),
                panel.background = element_rect(fill = "transparent")) + 
                guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))

    
    p.list2[[r]] <- ggplot(subset(data.df.avg,risk==2 & type=="SSP 5-8.5"), aes(x=year, y=range)) + geom_line(alpha=0.9,linewidth=1.5) + 
        xlab(NULL) + ylab (NULL) + geom_ribbon(aes(ymin=low,ymax=up),alpha=0.3,linewidth=1.5) + facet_wrap(~ season,scales = "free",nrow=1,ncol=4,labeller = label_parsed) +
        ggtitle(paste0(region.name[r],"; risk functional 2")) +
        theme(axis.text = element_text(size=16,face="bold"),
                strip.text = element_text(size = 16),
                plot.title = element_text(size=16,face="bold",hjust=0.5),
                axis.ticks =  element_line(size = 1.5),
                panel.border = element_rect(fill = "transparent", # Needed to add the border
                                            color = "black",            # Color of the border
                                            linewidth = 1),
                panel.background = element_rect(fill = "transparent")) + 
        guides(colour = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")),linetype = guide_legend(title.theme = element_text(size = 16, face = "bold"), label.theme = element_text(size = 16,face="bold"),override.aes = list(size = 2),keywidth = unit(1.5,"cm")))

}

pdf("figures/tail_correlation_range_uncertainty.pdf",width = 5*4,height = 3*8,onefile = TRUE)
ggarrange(plotlist=p.list1,nrow=8,ncol=1,common.legend=TRUE,legend="bottom")
ggarrange(plotlist=p.list2,nrow=8,ncol=1,common.legend=TRUE,legend="bottom")
dev.off()

##############################################################
######### Simulations ########################################
##############################################################
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

