rm(list=ls())


library(dplyr)
library(robis)
library(readr)
library(ggplot2)
library(vegan)
library(mregions)
library(rgeos)
library(rgdal) # I installed this package and now the iho (line 15) works for me? who knows why


eez <- mr_names("MarineRegions:eez")
iho <- mr_names('MarineRegions:iho')
fao <- mr_names('MarineRegions:fao')
lme <- mr_names('MarineRegions:lme')


#Get the polygon for the Bass Strait
BS<-mr_shp(key = "MarineRegions:iho", filter= "Bass Strait", read= TRUE)
BassStrait<-mr_as_wkt(SJ,fmt=3)
p = readWKT(paste(BassStrait))  #I don't know what this is doing - I think moving it into a different format? Look up the function.
BS2<-gSimplify(p,tol=3) #Simplifying the polygon. Again, look up the function and play around with this to see what TOL values work.
#note that in the previous version I had tol = 5 and that didn't work for me; tol=3 did.
BS3<-writeWKT(BS2) #Putting back into WKT so that it works with occurrence

#Get the occurrence in that polygon but failed -- AEC notes: it worked!!!
occurrence(geometry= BS3)

#try mapping the polygon
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = BS2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)

