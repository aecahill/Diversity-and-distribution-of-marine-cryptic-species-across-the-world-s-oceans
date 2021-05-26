rm(list=ls())

library(dplyr)
library(robis)
library(readr)
library(rinat)
library(leaflet)
library(ggplot2)
library(vegan)
library(abdiv)

setwd("~/Desktop/Summer research/2021")
specieslist <- read_csv("specieslist.csv")

#Occurrence in ABNJ
a<-area()
ABNJid<- c(1, 6, 114, 176, 178, 222, 223, 225)
ABNJ <- occurrence(areaid= ABNJid)
map_leaflet(ABNJ)

Occ_ABNJ_Arctic= NULL
Occ_ABNJ_Indian= NULL
Occ_ABNJ_NAtlantic= NULL
Occ_ABNJ_SAtlantic= NULL
Occ_ABNJ_NPacific= NULL
Occ_ABNJ_SPacific= NULL
Occ_ABNJ_Southern= NULL



#Extract the occurrence in the ABNJ
#Only extract "scientificName", "decimalLatitude", "decimalLongitude" for the sake of time. 
#(Can add more if needed)
colname <- c("scientificName", "decimalLatitude", "decimalLongitude")

for (x in specieslist$SpeciesName) {
  Occ_ABNJ_Arctic<-rbind(Occ_ABNJ_Arctic, occurrence(areaid= 6,scientificname = x, fields = colname ))
  
  Occ_ABNJ_Indian<-rbind(Occ_ABNJ_Indian, occurrence(areaid= 114,scientificname = x, fields = colname ))
  
  Occ_ABNJ_NAtlantic<-rbind(Occ_ABNJ_NAtlantic, occurrence(areaid= 176,scientificname = x, fields = colname ))
  
  Occ_ABNJ_SAtlantic<-rbind(Occ_ABNJ_SAtlantic, occurrence(areaid= 222,scientificname = x, fields = colname ))
  
  Occ_ABNJ_NPacific<-rbind(Occ_ABNJ_NPacific, occurrence(areaid= 178,scientificname = x, fields = colname ))
  
  Occ_ABNJ_SPacific<-rbind(Occ_ABNJ_SPacific, occurrence(areaid= 225,scientificname = x, fields = colname ))
  
  Occ_ABNJ_Southern<-rbind(Occ_ABNJ_Southern, occurrence(areaid= 223,scientificname = x, fields = colname ))
}

Occ_ABNJ_Arctic$Ocean <-"Artic"
Occ_ABNJ_Indian$Ocean <-"Indian"
Occ_ABNJ_NAtlantic$Ocean <- "North Atlantic"
Occ_ABNJ_SAtlantic$Ocean <-"South Atlantic"
Occ_ABNJ_NPacific$Ocean <-"North Pacific"
Occ_ABNJ_SPacific$Ocean <-"South Pacific"
Occ_ABNJ_Southern$Ocean <-"Southern"

Occ_ABNJ<-rbind(Occ_ABNJ_Arctic, Occ_ABNJ_Indian,Occ_ABNJ_NAtlantic, 
      Occ_ABNJ_SAtlantic, Occ_ABNJ_NPacific, Occ_ABNJ_SPacific, Occ_ABNJ_Southern)
write.csv(Occ_ABNJ,"~/Desktop/Summer research/2021/OccurrenceABNJ.csv") 

#Calculate the Shannon diversity
nspecies<-ncol(checklist_ABNJ)
checklist_ABNJ<-table(Occ_ABNJ$Ocean, Occ_ABNJ$scientificName)
diversity<- diversity(checklist_ABNJ)

#calculate the margalef richness
a<-matrix(t(checklist_ABNJ), ncol= 7)
collumnname<-c("Arctic", "Indian", "NAtlantic", "SAtlantic", "NPacific", "SPacific", "Southern")
colnames(a)<-collumnname
d<-as.data.frame(a)

richness=NULL
for(x in collumnname){
richness<-rbind(richness,margalef(select(d,x)))
}

#Make a table
table_ABNJ<- data.frame(diversity, richness)
class(table_ABNJ)
write.csv(table_ABNJ,"~/Desktop/Summer research/2021/ABNJtable.csv")
