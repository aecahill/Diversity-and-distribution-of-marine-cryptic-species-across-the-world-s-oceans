rm(list=ls())

library(dplyr)
library(robis)
library(readr)
library(rinat)
library(leaflet)
library(ggplot2)
library(vegan)
library(abdiv)


setwd("~/Desktop/Summer research/2021/Github")
specieslist <- read_csv("specieslist.csv")

Checklist_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_AtlanticOcean.csv")
Checklist_IndianOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_IndianOcean.csv")
Checklist_PacificOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_PacificOcean.csv")
Checklist_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_ArcticOcean.csv")
Checklist_SouthernOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_SouthernOcean.csv")

Checklist_ArcticOcean$Ocean <-"Artic"
Checklist_IndianOcean$Ocean <-"Indian"
Checklist_AtlanticOcean$Ocean <-"Atlantic"
Checklist_PacificOcean$Ocean <-"Pacific"
Checklist_SouthernOcean$Ocean <-"Southern"
Oceans_checklist<-rbind(Checklist_ArcticOcean, Checklist_IndianOcean, Checklist_AtlanticOcean, Checklist_PacificOcean,
                        Checklist_SouthernOcean)



#Calculate the Shannon diversity
Checklist<-table(Oceans_checklist$Ocean, Oceans_checklist$scientificName)
diversity<- diversity(Checklist)
nspecies<-ncol(Checklist)

#calculate the margalef richness
a<-matrix(t(Checklist), ncol= 5)
collumnname<-c("Arctic","Atlantic","Indian","Pacific","Southern")
colnames(a)<-collumnname
d<-as.data.frame(a)

richness=NULL
for(x in collumnname){
  richness<-rbind(richness,margalef(select(d,x)))
}

#Make a table
table_Oceans<- data.frame(diversity, richness)
write.csv(table_Oceans,"~/Desktop/Summer research/2021/table_Oceans.csv")



#Get the cryptic species only
Checklist_ArcticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_ArcticOcean_cryptic<-rbind(Checklist_ArcticOcean_cryptic,filter(Checklist_ArcticOcean, scientificName==x))
}

Checklist_IndianOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_IndianOcean_cryptic<-rbind(Checklist_IndianOcean_cryptic,filter(Checklist_IndianOcean, scientificName==x))
}

Checklist_AtlanticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_AtlanticOcean_cryptic<-rbind(Checklist_AtlanticOcean_cryptic,filter(Checklist_AtlanticOcean, scientificName==x))
}

Checklist_PacificOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_PacificOcean_cryptic<-rbind(Checklist_PacificOcean_cryptic,filter(Checklist_PacificOcean, scientificName==x))
}

Checklist_SouthernOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_SouthernOcean_cryptic<-rbind(Checklist_SouthernOcean_cryptic,filter(Checklist_SouthernOcean, scientificName==x))
}

Oceans_checklist_cryptic<-rbind(Checklist_ArcticOcean_cryptic,Checklist_IndianOcean_cryptic,Checklist_AtlanticOcean_cryptic,
                                Checklist_PacificOcean_cryptic,Checklist_SouthernOcean_cryptic)

#Calculate the Shannon diversity
Checklist_cryptic<-table(Oceans_checklist_cryptic$Ocean, Oceans_checklist_cryptic$scientificName)
diversity<- diversity(Checklist_cryptic)
nspecies<-ncol(Checklist_cryptic)

#calculate the margalef richness
b<-matrix(t(Checklist_cryptic), ncol= 5)
collumnname<-c("Arctic","Atlantic","Indian","Pacific","Southern")
colnames(b)<-collumnname
e<-as.data.frame(b)

richness=NULL
for(x in collumnname){
  richness<-rbind(richness,margalef(select(e,x)))
}

#Make a table
table_Oceans_cryptic<- data.frame(diversity, richness)
write.csv(table_Oceans_cryptic,"~/Desktop/Summer research/2021/table_Oceans_cryptic.csv")

#Calculate beta diversity over cells approximately 500 miles across
#example English channel
colname <- c("scientificName", "decimalLatitude", "decimalLongitude")
English_Channel<- occurrence( geometry = World_Oceans$coor[16], fields = colname )

dggs<- dgconstruct(spacing=500, metric=FALSE, resround='down')
English_Channel$cell <- dgtransform(dggs,English_Channel$decimalLatitude,English_Channel$decimalLongitude)
English_Channel_Checklist<-table(English_Channel$cell, English_Channel$scientificName)
betadiv_English_Channel<-betadiver(English_Channel_Checklist,"w")
a<-as.numeric(betadiv_English_Channel)
betadiv_English_Channe_mean<-mean(a)


#Arctic Ocean beta diversity
Occ_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/Occ_ArcticOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_ArcticOcean$cell <- dgtransform(dggs,Occ_ArcticOcean$decimalLatitude,Occ_ArcticOcean$decimalLongitude)

Occ_ArcticOcean_Checklist<-table(Occ_ArcticOcean$cell, Occ_ArcticOcean$scientificName)

betadiv_ArcticOcean<-betadiver(Occ_ArcticOcean_Checklist,"w")
betadiv_ArcticOcean_mean<-mean(as.numeric(betadiv_ArcticOcean))

#Indian Ocean beta diversity
Occ_IndianOcean <- read_csv("World's Oceans/Occ Oceans/Occ_IndianOcean.csv")

Occ_IndianOcean$cell <- dgtransform(dggs,Occ_ArcticOcean$decimalLatitude,Occ_ArcticOcean$decimalLongitude)

Occ_IndianOcean_Checklist<-table(Occ_IndianOcean$cell, Occ_IndianOcean$scientificName)

betadiv_IndianOcean<-betadiver(Occ_IndianOcean_Checklist,"w")
betadiv_IndianOcean_mean<-mean(as.numeric(betadiv_IndianOcean))

#Southern Ocean beta diversity
Occ_SouthOcean <- read_csv("World's Oceans/Occ Oceans/Occ_SouthOcean.csv")

Occ_SouthOcean$cell <- dgtransform(dggs,Occ_SouthOcean$decimalLatitude,Occ_SouthOcean$decimalLongitude)

Occ_SouthOcean_Checklist<-table(Occ_SouthOcean$cell, Occ_SouthOcean$scientificName)

betadiv_SouthOcean<-betadiver(Occ_SouthOcean_Checklist,"w")
betadiv_SouthOcean_mean<-mean(as.numeric(betadiv_SouthOcean))


#Atlantic Ocean beta diversity
Occ_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/Occ_AtlanticOcean.csv")

Occ_AtlanticOcean$cell <- dgtransform(dggs,Occ_AtlanticOcean$decimalLatitude,Occ_AtlanticOcean$decimalLongitude)

Occ_AtlanticOcean_Checklist<-table(Occ_AtlanticOcean$cell, Occ_AtlanticOcean$scientificName)

betadiv_AtlanticOcean<-betadiver(Occ_AtlanticOcean_Checklist,"w")
betadiv_AtlanticOcean_mean<-mean(as.numeric(betadiv_AtlanticOcean))

#Pacific Ocean beta diversity
Occ_PacificOcean <- read_csv("World's Oceans/Occ Oceans/Occ_PacificOcean.csv")

Occ_PacificOcean$cell <- dgtransform(dggs,Occ_PacificOcean$decimalLatitude,Occ_PacificOcean$decimalLongitude)

Occ_PacificOcean_Checklist<-table(Occ_PacificOcean$cell, Occ_PacificOcean$scientificName)

betadiv_PacificOcean<-betadiver(Occ_PacificOcean_Checklist,"w")
betadiv_PacificOcean_mean<-mean(as.numeric(betadiv_PacificOcean))
