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

