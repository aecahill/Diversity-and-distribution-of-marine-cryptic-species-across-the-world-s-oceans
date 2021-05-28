rm(list=ls())


library(dplyr)
library(robis)
library(readr)
library(ggplot2)
library(vegan)
library(mregions)



eez <- mr_names("MarineRegions:eez")
iho <- mr_names('MarineRegions:iho')
fao <- mr_names('MarineRegions:fao')
lme <- mr_names('MarineRegions:lme')


#Get the polygon for the North sea
NS<-mr_shp(key = "MarineRegions:lme", filter= "North Sea", read= TRUE)
Northsea<-mr_as_wkt(NS)

#Get the occurrence in that polygon but failed
occurrence(geometry= Northsea )

#try mapping the polygon
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = NS, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)
