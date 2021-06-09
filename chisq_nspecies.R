rm(list=ls())

library(dplyr)
library(robis)
library(readr)
library(rinat)
library(leaflet)
library(ggplot2)
library(vegan)
library(abdiv)
library(sf)


setwd("~/Desktop/Summer research/2021/Github")
specieslist <- read_csv("specieslist.csv")


#Import checklists
Checklist_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_AtlanticOcean.csv")
Checklist_IndianOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_IndianOcean.csv")
Checklist_PacificOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_PacificOcean.csv")
Checklist_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_ArcticOcean.csv")
Checklist_SouthernOcean <- read_csv("World's Oceans/Occ Oceans/Checklist_SouthernOcean.csv")

Checklist_ArcticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_ArcticOcean_cryptic<-rbind(Checklist_ArcticOcean_cryptic,filter(Checklist_ArcticOcean, scientificName==x))}

Checklist_IndianOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_IndianOcean_cryptic<-rbind(Checklist_IndianOcean_cryptic,filter(Checklist_IndianOcean, scientificName==x))}

Checklist_AtlanticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_AtlanticOcean_cryptic<-rbind(Checklist_AtlanticOcean_cryptic,filter(Checklist_AtlanticOcean, scientificName==x))}

Checklist_PacificOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_PacificOcean_cryptic<-rbind(Checklist_PacificOcean_cryptic,filter(Checklist_PacificOcean, scientificName==x))}

Checklist_SouthernOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Checklist_SouthernOcean_cryptic<-rbind(Checklist_SouthernOcean_cryptic,filter(Checklist_SouthernOcean, scientificName==x))}


#Chi-square test
nspecies<-c(nrow(Checklist_ArcticOcean),nrow(Checklist_AtlanticOcean),nrow(Checklist_IndianOcean),
            nrow(Checklist_PacificOcean),nrow(Checklist_SouthernOcean))

nspecies_cryptic<-c(nrow(Checklist_ArcticOcean_cryptic),nrow(Checklist_AtlanticOcean_cryptic),nrow(Checklist_IndianOcean_cryptic),
            nrow(Checklist_PacificOcean_cryptic),nrow(Checklist_SouthernOcean_cryptic))

table_nspecies<-t(data.frame(nspecies,nspecies_cryptic))
collumnname<-c("Arctic","Atlantic","Indian","Pacific","Southern")
colnames(table_nspecies)<-collumnname

chisq.test(table_nspecies)


#Making tables of percentages to compare
table_nspecies_percent<-table_nspecies
table_nspecies_percent[1,] <- table_nspecies[1,] / sum(table_nspecies[1,])
table_nspecies_percent[2,] <- table_nspecies[2,] / sum(table_nspecies[2,])

write.csv(table_nspecies,"~/Desktop/Summer research/2021/table_nspecies.csv")
write.csv(table_nspecies_percent,"~/Desktop/Summer research/2021/table_nspecies_percent.csv")

           