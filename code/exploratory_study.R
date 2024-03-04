rm(list=ls())
load("data/precip.RData")
load("data/temperature.RData")
date = seq(START.date,END.date, 1)

### count the number of missing values ###
na.count <-sapply(precip, function(data) apply(sapply(data,function(x){is.na(x)}),
			1,sum,na.rm=TRUE)/length(data)*100)
plot(date,na.count[,1],type='l',ylim=c(0,100),ylab="Percentage of Missing Values")
legend("left",region.name,col=1:8,lwd=1,bty="n")
for(i in 2:8){
lines(date,na.count[,i],col=i)
}

### plot the summaries ### 
library(ggplot2)
library(cowplot)
library(dplyr)
library(lubridate)
library(viridis)

p1.list <- p2.list <- list()
mean.func <- function(x,y=TRUE){
    n <- sum(!is.na(x))
    if(y & n >= 52){
        return(mean(x,na.rm=TRUE))
    }
    if(!y & n >= 10){
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
    mean.year <- aggregate(idx.data$precip,by=list(idx.data$year,idx.data$id),FUN=mean.func,y=TRUE)
    mean.day <- aggregate(idx.data$precip,by=list(idx.data$yday,idx.data$id),FUN=mean.func,y=FALSE)  
    mean.day$elev <- rep(idx.elev,each=366)
    names(mean.day) <- c("day","ID","precip","elev")
    mean.year$elev <- rep(idx.elev,each=length(unique(mean.year$Group.1)))
    names(mean.year) <- c("year","ID","precip","elev")
    color_breaks = quantile(mean.year$precip,c(1:10)/10,na.rm=TRUE)
    main = paste0("Yearly average: ",idx.region.name)

    p1 <- ggplot(mean.year) + geom_tile(aes(x=year,y=ID,fill=precip)) 

    p1 <- p1 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,nD,length.out=5),
    labels=as.character(idx.elev[order(idx.elev)][seq(1,nD,length.out=5)]),limits=c(0,nD))

    p1 <- p1 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev = TRUE),
    limits=c(0,max(color_breaks)),
    values=scales::rescale(color_breaks, to = c(0,1), 
    from = c(0,max(color_breaks)))) 
    p1 <- p1 + ggtitle(main) + xlab("Year") + 
    theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "transparent",            # Color of the border
                                    linewidth = 0.5))

    color_breaks = quantile(mean.day$precip,c(1:10)/10,na.rm=TRUE)
    main = paste0("Daily average: ",idx.region.name)
    p2 <- ggplot(mean.day) + geom_tile(aes(x=day,y=ID,fill=precip)) 
    p2 <- p2 + scale_y_continuous(name="Elevation (m)",breaks=seq(1,nD,length.out=5),
    labels=as.character(idx.elev[order(idx.elev)][seq(1,nD,length.out=5)]),
    limits=c(0,nD))
    p2 <- p2 + scale_fill_gradientn(name="Precipitation",colours=topo.colors(length(color_breaks),rev=TRUE),
    limits=c(0,max(color_breaks)),
    values=scales::rescale(color_breaks, to = c(0,1), 
    from = c(0, max(color_breaks)))) 
    p2 <- p2 + ggtitle(main) + xlab("Day") + 
    theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "transparent",            # Color of the border
                                    linewidth = 0.5))
    p1.list[[idx]] <- p1; p2.list[[idx]] <- p2
}

pdf("figures/heatmaps.pdf",width = 6,height = 5,onefile = TRUE)
for(i in 1:length(region.id)){
    show(p1.list[[i]])
    show(p2.list[[i]])
}
dev.off()  

save(p1.list,p2.list,file="data/plot_heatmaps.RData")
 






