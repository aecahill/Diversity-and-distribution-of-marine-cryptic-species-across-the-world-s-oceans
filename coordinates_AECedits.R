rm(list=ls())


library(dplyr)
library(robis)
library(readr)
library(ggplot2)
library(vegan)
library(mregions)
library(rgeos)



eez <- mr_names("MarineRegions:eez")
iho <- mr_names('MarineRegions:iho')
fao <- mr_names('MarineRegions:fao')
lme <- mr_names('MarineRegions:lme')


#Get the polygon for the Sea of Japan
SJ<-mr_shp(key = "MarineRegions:lme", filter= "Sea of Japan", read= TRUE)
SeaofJapan<-mr_as_wkt(SJ,fmt=3)
p = readWKT(paste(SeaofJapan))  #I don't know what this is doing - I think moving it into a different format? Look up the function.
SeaJapan2<-gSimplify(p,tol=5) #Simplifying the polygon. Again, look up the function and play around with this to see what TOL values work.
SeaJapan3<-writeWKT(SeaJapan2) #Putting back into WKT so that it works with occurrence

#Get the occurrence in that polygon but failed -- AEC notes: it worked!!!
occurrence(geometry= SeaJapan3)

#try mapping the polygon
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = SeaJapan2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)

#note that I had to use SeaJapan2, not SeaJapan3 (not sure why - wrong data format I)
