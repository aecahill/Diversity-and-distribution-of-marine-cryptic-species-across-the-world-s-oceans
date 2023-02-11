rm(list=ls())


library(dplyr)
library(robis)
library(readr)
library(ggplot2)
library(vegan)
library(mregions)
library(rgeos)
library(sf)

setwd("~/Desktop/Summer research/2021")

#Notes: There are some parts of the code where I extract the coordinate from IHO and FAO
#then I changed the coordinates based on mregions map and google earth to cover the desired areas
#so the results you get from the code will not be the same as the one I got.
#I'll leave a note says "changed" at those places.



#Import the shp file of oceans
Oceans_boundary <- st_read("World_Seas_IHO_v3/World_Seas_IHO_v3.shp")


#Extract the coordinates for the Southern ocean

SouthOcean<-writeWKT(gSimplify(as_Spatial(Oceans_boundary[63,]), tol=1))

SouthOcean_coor<-cbind("South Ocean",SouthOcean)
SouthOcean_coor<- as.data.frame(SouthOcean_coor)
colnames(SouthOcean_coor)<- c("Regions", "coor")

write.csv(SouthOcean_coor,"~/Desktop/Summer research/2021/Github/SouthOcean_coor.csv")

SouthOcean_coor <- read_csv("Github/SouthOcean_coor.csv")

#Extract the coordinates for the Atlantic ocean
#South Atlantic changed

AtlanticRegionNames<-c("South Atlantic","Caribbean Sea","Labrador Sea","Irish Sea and St. Georges Channel",
                       "Inner Seas off the West Coast of Scotland","Adriatic Sea","Gulf of Riga","Baltic Sea",
                       "Gulf of Finland","Gulf of Bothnia","Gulf of Mexico","North Atlantic","Gulf of St. Lawrence","Bay of Biscay",
                       "Celtic Sea","English Channel","North Sea","Kattegat","Skagerrak","Gulf of Guinea ","Mediterranean and Black Sea" )
AtlanticRegions<-c(62,22,27,33,34,42,56,57,58,59,75,76,77,79,80,85,88,96,98,101)
Atlantic_coor=NULL
for(x in AtlanticRegions){
  A<-gSimplify(as_Spatial(Oceans_boundary[x,]), tol=1)
  Atlantic_coor<-rbind(Atlantic_coor, writeWKT(A))
}

NAtlantic_Mediterranean_BlackSea<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Mediterranean and Black Sea", read= TRUE), tol=1))
Atlantic_coor<-rbind(Atlantic_coor, NAtlantic_Mediterranean_BlackSea)

Atlantic_coor<-cbind(AtlanticRegionNames,Atlantic_coor)
Atlantic_coor<- as.data.frame(Atlantic_coor)
colnames(Atlantic_coor)<- c("Regions", "coor")
write.csv(Atlantic_coor,"~/Desktop/Summer research/2021/Github/Atlantic_coor.csv")

Atlantic_coor <- read_csv("Github/Atlantic_coor.csv")

#Extract the coordinates for Pacific Oceans

PacificRegions<- c("Pacific, Northeast","Pacific, Eastern Central","Pacific, Southeast","Pacific, Southwest",
                   "Pacific, Western Central","Pacific, Nothwest")
Pacific_coor=NULL
for (x in PacificRegions) {
  A<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= x, read= TRUE), tol=1)
  Pacific_coor<-rbind(Pacific_coor, writeWKT(A))
}
Pacific_coor<-cbind(PacificRegions,Pacific_coor)
Pacific_coor<-as.data.frame(Pacific_coor)
colnames(Pacific_coor)<- c("Regions", "coor")

write.csv(Pacific_coor,"~/Desktop/Summer research/2021/Github/Pacific_coor.csv")

Pacific_coor <- read_csv("Github/Pacific_coor.csv")

#Extract the coordinates for Indian Ocean
#Changed

IndianRegions<-c("Indian Ocean, Western","Indian Ocean, Eastern")
ID1<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Indian Ocean, Western", read= TRUE), tol=1))
ID2<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Indian Ocean, Eastern", read= TRUE), tol=1))
Indian_coor<-rbind(ID1, ID2)
Indian_coor<-cbind(IndianRegions,Indian_coor)
Indian_coor<-as.data.frame(Indian_coor)
colnames(Indian_coor)<- c("Regions", "coor")

Indian_coor <- read_csv("Github/Indian_coor.csv")

#Extract the coordinates for Arctic Ocean
#changed

ArcticRegions<-c("Arctic, Central West","Arctic, Central East", "Arctic Sea")

ArcticWest<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Atlantic, Northwest", read= TRUE), tol=1))
ArcticEast<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Atlantic, Northeast", read= TRUE), tol=1))
ArcticSea<-writeWKT(gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Arctic Sea", read= TRUE), tol=1))

Arctic_coor<-rbind(ArcticWest, ArcticEast, ArcticSea)
Arctic_coor<-cbind(ArcticRegions,Arctic_coor)
Arctic_coor<-as.data.frame(Arctic_coor)
colnames(Arctic_coor)<- c("Regions", "coor")

Arctic_coor <- read_csv("Github/Arctic_coor.csv")


#Put the oceans together
World_Oceans<-rbind(Atlantic_coor,Arctic_coor,Indian_coor,Pacific_coor,SouthOcean_coor)
write.csv(World_Oceans,"~/Desktop/Summer research/2021/Github/World_Oceans.csv")

#Try map the oceans to double check
SO = readWKT(paste(SouthOcean_coor$coor[1]))
IN1 = readWKT(paste(Indian_coor$coor[1]))
IN2 = readWKT(paste(Indian_coor$coor[2]))
ARC1= readWKT(paste(Arctic_coor$coor[1]))
ARC2= readWKT(paste(Arctic_coor$coor[2]))
ARC3= readWKT(paste(Arctic_coor$coor[3]))
PAC1= readWKT(paste(Pacific_coor$coor[1]))
PAC2= readWKT(paste(Pacific_coor$coor[2]))
PAC3= readWKT(paste(Pacific_coor$coor[3]))
PAC4= readWKT(paste(Pacific_coor$coor[4]))
PAC5= readWKT(paste(Pacific_coor$coor[5]))
PAC6= readWKT(paste(Pacific_coor$coor[6]))
ATL1= readWKT(paste(Atlantic_coor$coor[1]))
ATL2= readWKT(paste(Atlantic_coor$coor[2]))
ATL3= readWKT(paste(Atlantic_coor$coor[3]))
ATL4= readWKT(paste(Atlantic_coor$coor[4]))
ATL5= readWKT(paste(Atlantic_coor$coor[5]))
ATL6= readWKT(paste(Atlantic_coor$coor[6]))
ATL7= readWKT(paste(Atlantic_coor$coor[7]))
ATL8= readWKT(paste(Atlantic_coor$coor[8]))
ATL9= readWKT(paste(Atlantic_coor$coor[9]))
ATL10= readWKT(paste(Atlantic_coor$coor[10]))
ATL11= readWKT(paste(Atlantic_coor$coor[11]))
ATL12= readWKT(paste(Atlantic_coor$coor[12]))
ATL13= readWKT(paste(Atlantic_coor$coor[13]))
ATL14= readWKT(paste(Atlantic_coor$coor[14]))
ATL15= readWKT(paste(Atlantic_coor$coor[15]))
ATL16= readWKT(paste(Atlantic_coor$coor[16]))
ATL17= readWKT(paste(Atlantic_coor$coor[17]))
ATL18= readWKT(paste(Atlantic_coor$coor[18]))
ATL19= readWKT(paste(Atlantic_coor$coor[19]))
ATL20= readWKT(paste(Atlantic_coor$coor[20]))
ATL21= readWKT(paste(Atlantic_coor$coor[21]))

world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = SO, aes(x = long, y = lat, group = group), fill = "#b400cc") +
  geom_polygon(data = IN1, aes(x = long, y = lat, group = group), fill = "#cccc00") +
  geom_polygon(data = IN2, aes(x = long, y = lat, group = group), fill = "#cccc00") +
  geom_polygon(data = ARC1, aes(x = long, y = lat, group = group), fill = "#58cc00") +
  geom_polygon(data = ARC2, aes(x = long, y = lat, group = group), fill = "#58cc00") +
  geom_polygon(data = ARC3, aes(x = long, y = lat, group = group), fill = "#58cc00") +
  geom_polygon(data = PAC1, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = PAC2, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = PAC3, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = PAC4, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = PAC5, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = PAC6, aes(x = long, y = lat, group = group), fill = "#00c5cc") +
  geom_polygon(data = ATL1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL3, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL4, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL5, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL6, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL7, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL8, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL9, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL10, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL11, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL12, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL13, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL14, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL15, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL16, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL17, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL18, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL19, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL20, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = ATL21, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)



