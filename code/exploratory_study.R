rm(list=ls())
load("data/precip.RData")
load("data/temperature.RData")
date = seq(START.date,END.date, 1)

### count the number of missing values ###
na.count <-sapply(precip, function(data) apply(sapply(data,function(x){is.na(x)}),
			1,sum,na.rm=TRUE)/length(data)*100)
plot(date,na.count[,1],type='l',ylim=c(0,100),ylab="Percentage of Missing Values")
legend("left",region.name,col=1:8,lwd=1)
for(i in 2:8){
lines(date,na.count[,i],col=i)
}

### plot the summaries ### 
library(ggplot2)
library(cowplot)








