library(vegan)
library(mregions)
library(rgdal)
library(rgeos)
library(sf)
library(robis)
library(dplyr)

World_Oceans<-read_csv("World_Oceans.csv")

#Atlantic
Geometry_AtlanticOcean<-World_Oceans$coor[1:21]
Checklist_AtlanticOcean=NULL
for(x in Geometry_AtlanticOcean){
  ATL<-checklist(geometry = x)
  ATL<- as.data.frame(ATL)
  ATL<-subset(ATL, select = c(scientificName,records))
  Checklist_AtlanticOcean<-rbind(Checklist_AtlanticOcean,ATL)
}

Checklist_AtlanticOcean<-group_by(Checklist_AtlanticOcean, scientificName)
Checklist_AtlanticOcean<-Checklist_AtlanticOcean %>% summarise(
  records = sum(records)
  )
write_csv(Checklist_AtlanticOcean, "Checklist_AtlanticOcean.csv")

#Pacific
Geometry_PacificOcean<-World_Oceans$coor[27:32]
Checklist_PacificOcean=NULL
for(x in Geometry_PacificOcean){
  PA<-checklist(geometry = x)
  PA<- as.data.frame(PA)
  PA<-subset(PA, select = c(scientificName,records))
  Checklist_PacificOcean<-rbind(Checklist_PacificOcean,PA)
}

Checklist_PacificOcean<-group_by(Checklist_PacificOcean, scientificName)
Checklist_PacificOcean<-Checklist_PacificOcean %>% summarise(
  records = sum(records)
)
write_csv(Checklist_PacificOcean, "Checklist_PacificOcean.csv")


#Indian
Geometry_IndianOcean<-World_Oceans$coor[25:26]
Checklist_IndianOcean=NULL
for(x in Geometry_IndianOcean){
  ID<-checklist(geometry = x)
  ID<- as.data.frame(ID)
  ID<-subset(ID, select = c(scientificName,records))
  Checklist_IndianOcean<-rbind(Checklist_IndianOcean,ID)
}

Checklist_IndianOcean<-group_by(Checklist_IndianOcean, scientificName)
Checklist_IndianOcean<-Checklist_IndianOcean %>% summarise(
  records = sum(records)
)
write_csv(Checklist_IndianOcean, "Checklist_IndianOcean.csv")


#Arctic
Geometry_ArcticOcean<-World_Oceans$coor[22:24]
Checklist_ArcticOcean=NULL
for(x in Geometry_ArcticOcean){
  ARC<-checklist(geometry = x)
  ARC<- as.data.frame(ARC)
  ARC<-subset(ARC, select = c(scientificName,records))
  Checklist_ArcticOcean<-rbind(Checklist_ArcticOcean,ARC)
}

Checklist_ArcticOcean<-group_by(Checklist_ArcticOcean, scientificName)
Checklist_ArcticOcean<-Checklist_ArcticOcean %>% summarise(
  records = sum(records)
)
write_csv(Checklist_ArcticOcean, "Checklist_ArcticOcean.csv")


#Southern
Geometry_SouthernOcean<-World_Oceans$coor[33]

Checklist_SouthernOcean<-checklist(geometry = Geometry_SouthernOcean)
Checklist_SouthernOcean<- as.data.frame(Checklist_SouthernOcean)
Checklist_SouthernOcean<-subset(Checklist_SouthernOcean, select = c(scientificName,records))

Checklist_SouthernOcean<-group_by(Checklist_SouthernOcean, scientificName)
Checklist_SouthernOcean<-Checklist_SouthernOcean %>% summarise(
  records = sum(records)
)
write_csv(Checklist_SouthernOcean, "Checklist_SouthernOcean.csv")

