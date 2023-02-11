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

eez <- mr_names("MarineRegions:eez")
iho <- mr_names('MarineRegions:iho')
fao <- mr_names('MarineRegions:fao')
lme <- mr_names('MarineRegions:lme')
?mr_names()

#Import the shp file of oceans
Oceans_boundary <- st_read("World_Seas_IHO_v3/World_Seas_IHO_v3.shp")


#Get the polygon for the Pacific


SW1<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Northeast", read= TRUE), tol=1)

SW2<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Eastern Central", read= TRUE), tol=1)

SW3<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Southeast", read= TRUE), tol=1)

SW4<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Western Central", read= TRUE), tol=1)

SW5<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Nothwest", read= TRUE), tol=1)

SW6<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Pacific, Southwest", read= TRUE), tol=1)


Pacific_coor_original <- read_csv("Github/Pacific_coor.csv")
PA1 = readWKT(paste(Oceans_coor_fixed$coor[1]))
PA2 = readWKT(paste(Oceans_coor_fixed$coor[2]))
PA3 = readWKT(paste(Oceans_coor_fixed$coor[3]))
PA4 = readWKT(paste(Oceans_coor_fixed$coor[4]))
PA5 = readWKT(paste(Oceans_coor_fixed$coor[5]))
PA6 = readWKT(paste(Oceans_coor_fixed$coor[6]))

SE<-mr_shp(key = "MarineRegions:fao", filter= "Pacific, Southwest", read= TRUE)
SE1<-gSimplify(SE, tol=1)
SouthWestAtlantic<-writeWKT(SE1)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = id1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = id2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)

#Area of north atlantic
NA1<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Mediterranean and Black Sea", read= TRUE), tol=1)

#try mapping the polygon
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_polygon(data = SW1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SW2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SW3, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SW4, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SW5, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SW6, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = id1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = id2, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = SO, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = A1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = G1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  geom_polygon(data = NA1, aes(x = long, y = lat, group = group), fill = "#cc3300") +
  coord_fixed(1)


Indiann <- read_csv("Github/Indian_coor.csv")
id1 = readWKT(paste(Indiann$coor[1]))
id2 = readWKT(paste(Indiann$coor[2]))



ID1<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Indian Ocean, Western", read= TRUE), tol=1)


ID2<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Indian Ocean, Eastern", read= TRUE), tol=1)


#Get the polygons for the Atlantic ocean
AtlanticRegions<- c("Atlantic, Antarctic","Atlantic, Southwest","Atlantic, Eastern Central","Atlantic, Southeast",
                  "Atlantic, Western Central","Atlantic, Northwest","Atlantic, Northeast")
coor_Atlantic=NULL
for (x in AtlanticRegions) {
  A<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= x, read= TRUE), tol=1)
  coor_Atlantic<-rbind(coor_Atlantic, writeWKT(A))
}
coor_Atlantic<-cbind(AtlanticRegions,coor_Atlantic)
coor_Atlantic<-as.data.frame(coor_Atlantic)
colnames(coor_Atlantic)<- c("Regions", "coor")

#Get the polygons for the Pacific ocean
PacificRegions<- c("Pacific, Northeast","Pacific, Eastern Central","Pacific, Southeast","Pacific, Southwest",
                    "Pacific, Antarctic","Pacific, Western Central","Pacific, Nothwest")
coor_Pacific=NULL
for (x in AtlanticRegions) {
  A<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= x, read= TRUE), tol=1)
  coor_Pacific<-rbind(coor_Pacific, writeWKT(A))
}
coor_Pacific<-cbind(PacificRegions,coor_Pacific)
coor_Pacific<-as.data.frame(coor_Pacific)
colnames(coor_Pacific)<- c("Regions", "coor")

#Get the polygons for the Indian ocean
IndianRegions<- c("Indian Ocean, Antarctic and Southern","Indian Ocean, Western","Indian Ocean, Eastern")
coor_Indian=NULL
for (x in IndianRegions) {
  A<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= x, read= TRUE), tol=1)
  coor_Indian<-rbind(coor_Indian, writeWKT(A))
}
coor_Indian<-cbind(IndianRegions,coor_Indian)
coor_Indian<-as.data.frame(coor_Indian)
colnames(coor_Indian)<- c("Regions", "coor")

#Get the polygon for the Arctic sea and Mediterranean
AS<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Arctic Sea", read= TRUE), tol=1)
ArcticSea<-writeWKT(AS)

M<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= "Mediterranean and Black Sea", read= TRUE), tol=1)
Mediterranean<-writeWKT(M)

SeaRegions<- c(ArcticSea,Mediterranean)
SeaRegionsname<- c("Arctic Sea","Mediterranean and Black Sea")
SeaRegions<-cbind(SeaRegionsname,SeaRegions)
SeaRegions<- as.data.frame(SeaRegions)
colnames(SeaRegions)<- c("Regions", "coor")

#Coor of all regions
Oceans_coor<-rbind(coor_Atlantic,coor_Pacific,coor_Indian, SeaRegions)
write.csv(Oceans_coor,"~/Desktop/Summer research/2021/Github/Oceans_coor.csv") 

coor_fao=NULL
for(x in fao$name){
  A<-gSimplify(mr_shp(key = "MarineRegions:fao", filter= x, read= TRUE), tol=1)
  coor_fao<-rbind(coor_fao, writeWKT(A))
}
coor_fao<-cbind(fao$name,coor_fao)
coor_fao<- as.data.frame(coor_fao)
colnames(coor_fao)<- c("Regions", "coor")
coor_fao<-arrange(coor_fao, Regions)
coor_fao <- coor_fao[-13,]

AA<-coor_fao[4,2]

