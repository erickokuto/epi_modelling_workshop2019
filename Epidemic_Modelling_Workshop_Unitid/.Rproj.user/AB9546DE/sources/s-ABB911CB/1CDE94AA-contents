
require(rgdal)
require(maptools)
rusinga <- readShapePoly("D:/workshops/UNITID_Modelling_workshop/day2/data/rusinga/rusinga.shp")
plot(rusinga)

points <- read.csv("D:/workshops/UNITID_Modelling_workshop/day2/data/points/hh.csv")
points = points[, c("latmean",  "lngmean")] 


# libraries ---------------------------------------------------------------



library(ggplot2)
library(dplyr)
library(sp)
library(rgdal)

plot(rusinga)
points <- SpatialPoints(points)
plot(points, add= TRUE , col = 'red')
