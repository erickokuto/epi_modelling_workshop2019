###

# libraries ---------------------------------------------------------------

library(rgdal)
library(plyr)
library(dplyr)

#import the shapefile
rusinga<- readOGR("D:/workshops/UNITID_Modelling_workshop/day2/data/rusinga", "rusinga")
plot(rusinga)

points.data<-read.csv("D:/R_projects/Epidemic_Modelling_Workshop_Unitid/data/hh.csv")
#write the points to a csv
#points.data<-write.csv(hh, "hh.csv")

#import the csv containing the points
#points<-read.csv("hh.csv")
points <-read.csv("D:/R_projects/Epidemic_Modelling_Workshop_Unitid/data/hh.csv")

require(ggplot2)
#fortify the shapefile
rusingafort<-fortify(rusinga, region="ADM4")
?fortify

#rename the fortified shapefile data
rusingafort$id<-ifelse(rusingafort$id%in%"RISINGA", "RUSINGA", rusingafort$id)

#plot the shapefile and the data
ggplot()+theme_bw()+
  geom_polygon(data=rusingafort, aes(y=lat, x=long), fill=NA, colour="black")+
  geom_point(data=points, aes(y=latmean, x=lngmean), color="red")+
  expand_limits(y=rusingafort$lat, x=rusingafort$long)

 
