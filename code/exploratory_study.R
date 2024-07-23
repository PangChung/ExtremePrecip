## code to generate the heatmaps for the precipitation data ##
rm(list=ls())
load("data/precip.RData")
load("data/temperature.RData")
date = seq(START.date,END.date, 1)

### count the number of missing values ###
# na.count <-sapply(precip, function(data) apply(sapply(data,function(x){is.na(x)}),
# 			1,sum,na.rm=TRUE)/length(data)*100)
# plot(date,na.count[,1],type='l',ylim=c(0,100),ylab="Percentage of Missing Values")
# legend("left",region.name,col=1:8,lwd=1,bty="n")
# for(i in 2:8){
# lines(date,na.count[,i],col=i)
# }

### plot the summaries ### 
library(ggplot2)
library(cowplot)
library(dplyr)
library(lubridate)
library(viridis)
library(RColorBrewer)

p1.list <- p2.list <- list()
mean.func <- function(x,y=TRUE){
    n <- sum(!is.na(x))
    if(y & n>=36){
        return(mean(x,na.rm=TRUE))
    }
    if(!y & n>=5){
        return(mean(x,na.rm=TRUE))
    }
    return(NA)
}

## prepare for the data.frame ##
for(idx in c(1:length(region.id))){
    idx.region.id = region.id[idx]
    idx.region.name = region.name[idx]
    idx.elev <- station$elev[station$group.id==idx.region.id]
    nD = length(idx.elev)
    idx.data <- data.frame(precip = unlist(precip[[idx]][order(idx.elev)]),    
            date=rep(date,nD),
            year=rep(year(date),nD),
            yday=rep(yday(date),nD),
            id = rep(1:nD,each=length(date)))
    idx.valid = !is.na(idx.data$precip) & idx.data$precip < 1000
    idx.data = idx.data[idx.valid, ]
    mean.year <- aggregate(idx.data$precip,by=list(idx.data$year,idx.data$id),FUN=mean.func,y=TRUE)
    mean.day <- aggregate(idx.data$precip,by=list(idx.data$yday,idx.data$id),FUN=mean.func,y=FALSE)  
    mean.day$elev <- idx.elev[mean.day$Group.2]
    names(mean.day) <- c("day","ID","precip","elev")
    mean.year$elev <- idx.elev[mean.year$Group.2]
    names(mean.year) <- c("year","ID","precip","elev")

    color_breaks = quantile(mean.year$precip,seq(1,10,length.out = 9)/10,na.rm=TRUE)+1
    main = paste0("Yearly average: ",idx.region.name)

    p1 <- ggplot(mean.year) + geom_tile(aes(x=year,y=ID,fill=precip)) 
    p1 <- p1 + scale_y_continuous(name=NULL,breaks=seq(1,nD,length.out=5),
    labels=as.character(idx.elev[order(idx.elev)][seq(1,nD,length.out=5)]),limits=c(0,nD))
    p1 <- p1 + scale_fill_gradientn(name=NULL,colours=brewer.pal(length(color_breaks),"BuPu"),breaks=ceiling(exp(seq(0,log1p(max(color_breaks)),length.out=5))-1),
    limits=c(0,max(color_breaks)),trans="log1p") 
    p1 <- p1 + xlab(NULL) + ylab(NULL)  + ggtitle(main) +
    theme(axis.text = element_text(size=16,face="bold"), 
        plot.title = element_text(size=16,face="bold",hjust=0.5),
        legend.position = "bottom",
        axis.ticks =  element_line(size = 2),
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "transparent",            # Color of the border
                                    linewidth = 0.5)) + guides(fill = guide_colourbar(barheight = unit(0.5, "cm"), barwidth = unit(5, "cm"), title.theme = element_text(size = 14, face = "bold"), label.theme = element_text(size = 14)))

    color_breaks = quantile(mean.day$precip,seq(1,10,length.out = 9)/10,na.rm=TRUE)+1    
    main = paste0("Daily average: ",idx.region.name)
    p2 <- ggplot(mean.day) + geom_tile(aes(x=day,y=ID,fill=precip)) 
    p2 <- p2 + scale_y_continuous(name=NULL,breaks=seq(1,nD,length.out=5),
    labels=as.character(idx.elev[order(idx.elev)][seq(1,nD,length.out=5)]),
    limits=c(0,nD)) + scale_x_continuous(name=NULL,breaks=seq(15,350,length.out=12),labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))
    p2 <- p2 + scale_fill_gradientn(name=NULL,colours=brewer.pal(length(color_breaks),"BuPu"),breaks=ceiling(exp(seq(0,log1p(max(color_breaks)),length.out=5))-1),
    limits=c(0,max(color_breaks)), trans="log1p") 
    p2 <- p2 + xlab(NULL) + ylab(NULL) + ggtitle(main) +
    theme(axis.text = element_text(size=16,face="bold"),
            plot.title = element_text(size=16,face="bold",hjust=0.5),
            axis.ticks =  element_line(size = 2),
            legend.position = "bottom",
            panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "transparent",            # Color of the border
                                    linewidth = 0.5)) + guides(fill = guide_colourbar(barheight = unit(0.5, "cm"), barwidth = unit(5, "cm"), title.theme = element_text(size = 14, face = "bold"), label.theme = element_text(size = 14)))
    p1.list[[idx]] <- p1; p2.list[[idx]] <- p2
}

pdf("figures/heatmaps.pdf",width = 16 ,height = 8,onefile = TRUE)
ggarrange(plotlist=p1.list,nrow=2,ncol=4,common.legend=FALSE)
ggarrange(plotlist=p2.list,nrow=2,ncol=4,common.legend=FALSE)
dev.off()  

save(p1.list,p2.list,file="data/plot_heatmaps.RData")
 







