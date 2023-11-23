rm(list=ls())
load("data/precip.RData")
library(sp)
zones = c(32,16,13,14,14)
loc.trans.list <- list()
for(idx in 1:5){
idx.loc = station$group.id == region.id[idx]
loc = cbind(station$Y[idx.loc],station$X[idx.loc])
points <- SpatialPoints(coords = loc, proj4string = CRS("+proj=longlat +datum=WGS84"))
# Transform the coordinates to Euclidean coordinates
euclidean_coords <- spTransform(points, CRS(paste0("+proj=utm ","+zone=",zones[idx]," +datum=WGS84")))
loc.trans.list[[idx]] = cbind(euclidean_coords@coords[,1],euclidean_coords@coords[,2])
}

save(loc.trans.list,zones,file="data/transformed_coordinates.RData")
