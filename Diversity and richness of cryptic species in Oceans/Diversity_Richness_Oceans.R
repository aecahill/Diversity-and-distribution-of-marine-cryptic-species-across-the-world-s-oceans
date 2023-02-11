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
library(dggridR)


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
table_Oceans<- data.frame( collumnname, diversity, richness)
colnames(table_Oceans)<- c("Ocean","diversity","richness")
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
table_Oceans_cryptic<- data.frame( collumnname, diversity, richness)
colnames(table_Oceans_cryptic)<- c("Ocean","diversity","richness")
write.csv(table_Oceans_cryptic,"~/Desktop/Summer research/2021/table_Oceans_cryptic.csv")

#Make charts
ggplot(table_Oceans, aes(x=Ocean, y= diversity)) +
  geom_col(width = 0.5) +
  theme_bw()

ggplot(table_Oceans, aes(x=Ocean, y= richness)) +
  geom_col(width = 0.5) +
  theme_bw()

ggplot(table_Oceans_cryptic, aes(x=Ocean, y= diversity)) +
  geom_col(width = 0.5) +
  theme_bw()

ggplot(table_Oceans_cryptic, aes(x=Ocean, y= richness)) +
  geom_col(width = 0.5) +
  theme_bw()


#Calculate beta diversity over cells approximately 1000 miles across
#Arctic Ocean beta diversity
Occ_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/Occ_ArcticOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_ArcticOcean$cell <- dgtransform(dggs,Occ_ArcticOcean$decimalLatitude,Occ_ArcticOcean$decimalLongitude)

Occ_ArcticOcean_Checklist<-table(Occ_ArcticOcean$cell, Occ_ArcticOcean$scientificName)

betadiv_ArcticOcean<-betadiver(Occ_ArcticOcean_Checklist,"w")
betadiv_ArcticOcean<-as.numeric(betadiv_ArcticOcean)
betadiv_ArcticOcean<-data.frame(betadiv_ArcticOcean)

write_csv(betadiv_ArcticOcean, "~/Desktop/Summer research/2021/betadiv_ArcticOcean.csv")



#Indian Ocean beta diversity
Occ_IndianOcean <- read_csv("World's Oceans/Occ Oceans/Occ_IndianOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_IndianOcean$cell <- dgtransform(dggs,Occ_IndianOcean$decimalLatitude,Occ_IndianOcean$decimalLongitude)

Occ_IndianOcean_Checklist<-table(Occ_IndianOcean$cell, Occ_IndianOcean$scientificName)

betadiv_IndianOcean<-betadiver(Occ_IndianOcean_Checklist,"w")
betadiv_IndianOcean<-as.numeric(betadiv_IndianOcean)
betadiv_IndianOcean<-data.frame(betadiv_IndianOcean)

write_csv(betadiv_IndianOcean, "~/Desktop/Summer research/2021/betadiv_IndianOcean.csv")


#Southern Ocean beta diversity
Occ_SouthOcean <- read_csv("World's Oceans/Occ Oceans/Occ_SouthOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_SouthOcean$cell <- dgtransform(dggs,Occ_SouthOcean$decimalLatitude,Occ_SouthOcean$decimalLongitude)

Occ_SouthOcean_Checklist<-table(Occ_SouthOcean$cell, Occ_SouthOcean$scientificName)

betadiv_SouthOcean<-betadiver(Occ_SouthOcean_Checklist,"w")
betadiv_SouthOcean<-as.numeric(betadiv_SouthOcean)
betadiv_SouthOcean<-data.frame(betadiv_SouthOcean)

write_csv(betadiv_SouthOcean, "~/Desktop/Summer research/2021/betadiv_SouthOcean.csv")


#Atlantic Ocean beta diversity
Occ_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/Occ_AtlanticOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_AtlanticOcean$cell <- dgtransform(dggs,Occ_AtlanticOcean$decimalLatitude,Occ_AtlanticOcean$decimalLongitude)

Occ_AtlanticOcean_Checklist<-table(Occ_AtlanticOcean$cell, Occ_AtlanticOcean$scientificName)

betadiv_AtlanticOcean<-betadiver(Occ_AtlanticOcean_Checklist,"w")
betadiv_AtlanticOcean<-as.numeric(betadiv_AtlanticOcean)
betadiv_AtlanticOcean<-data.frame(betadiv_AtlanticOcean)

write_csv(betadiv_AtlanticOcean, "~/Desktop/Summer research/2021/betadiv_AtlanticOcean.csv")

#Pacific Ocean beta diversity
Occ_PacificOcean <- read_csv("World's Oceans/Occ Oceans/Occ_PacificOcean.csv")

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_PacificOcean$cell <- dgtransform(dggs,Occ_PacificOcean$decimalLatitude,Occ_PacificOcean$decimalLongitude)

Occ_PacificOcean_Checklist<-table(Occ_PacificOcean$cell, Occ_PacificOcean$scientificName)

betadiv_PacificOcean<-betadiver(Occ_PacificOcean_Checklist,"w")
betadiv_PacificOcean<-as.numeric(betadiv_PacificOcean)
betadiv_PacificOcean<-data.frame(betadiv_PacificOcean)

write_csv(betadiv_PacificOcean, "~/Desktop/Summer research/2021/betadiv_PacificOcean.csv")

#Arctic Ocean cryptic species beta diversity
Occ_ArcticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Occ_ArcticOcean_cryptic<-rbind(Occ_ArcticOcean_cryptic,filter(Occ_ArcticOcean, scientificName==x))
}

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_ArcticOcean_cryptic$cell <- dgtransform(dggs,Occ_ArcticOcean_cryptic$decimalLatitude,Occ_ArcticOcean_cryptic$decimalLongitude)

Occ_ArcticOcean_cryptic_Checklist<-table(Occ_ArcticOcean_cryptic$cell, Occ_ArcticOcean_cryptic$scientificName)

betadiv_ArcticOcean_cryptic<-betadiver(Occ_ArcticOcean_cryptic_Checklist,"w")
betadiv_ArcticOcean_cryptic<-as.numeric(betadiv_ArcticOcean_cryptic)
betadiv_ArcticOcean_cryptic<-data.frame(betadiv_ArcticOcean_cryptic)

write_csv(betadiv_ArcticOcean_cryptic, "~/Desktop/Summer research/2021/betadiv_ArcticOcean_cryptic.csv")

#South Ocean cryptic species beta diversity
Occ_SouthOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Occ_SouthOcean_cryptic<-rbind(Occ_SouthOcean_cryptic,filter(Occ_SouthOcean, scientificName==x))
}

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_SouthOcean_cryptic$cell <- dgtransform(dggs,Occ_SouthOcean_cryptic$decimalLatitude,Occ_SouthOcean_cryptic$decimalLongitude)

Occ_SouthOcean_cryptic_Checklist<-table(Occ_SouthOcean_cryptic$cell, Occ_SouthOcean_cryptic$scientificName)

betadiv_SouthOcean_cryptic<-betadiver(Occ_SouthOcean_cryptic_Checklist,"w")
betadiv_SouthOcean_cryptic<-as.numeric(betadiv_SouthOcean_cryptic)
betadiv_SouthOcean_cryptic<-data.frame(betadiv_SouthOcean_cryptic)

write_csv(betadiv_SouthOcean_cryptic, "~/Desktop/Summer research/2021/betadiv_SouthOcean_cryptic.csv")

#Indian Ocean cryptic species beta diversity
Occ_IndianOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Occ_IndianOcean_cryptic<-rbind(Occ_IndianOcean_cryptic,filter(Occ_IndianOcean, scientificName==x))
}

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_IndianOcean_cryptic$cell <- dgtransform(dggs,Occ_IndianOcean_cryptic$decimalLatitude,Occ_IndianOcean_cryptic$decimalLongitude)

Occ_IndianOcean_cryptic_Checklist<-table(Occ_IndianOcean_cryptic$cell, Occ_IndianOcean_cryptic$scientificName)

betadiv_IndianOcean_cryptic<-betadiver(Occ_IndianOcean_cryptic_Checklist,"w")
betadiv_IndianOcean_cryptic<-as.numeric(betadiv_IndianOcean_cryptic)
betadiv_IndianOcean_cryptic<-data.frame(betadiv_IndianOcean_cryptic)

write_csv(betadiv_IndianOcean_cryptic, "~/Desktop/Summer research/2021/betadiv_IndianOcean_cryptic.csv")

#Pacific Ocean cryptic species beta diversity
Occ_PacificOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Occ_PacificOcean_cryptic<-rbind(Occ_PacificOcean_cryptic,filter(Occ_PacificOcean, scientificName==x))
}

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_PacificOcean_cryptic$cell <- dgtransform(dggs,Occ_PacificOcean_cryptic$decimalLatitude,Occ_PacificOcean_cryptic$decimalLongitude)

Occ_PacificOcean_cryptic_Checklist<-table(Occ_PacificOcean_cryptic$cell, Occ_PacificOcean_cryptic$scientificName)

betadiv_PacificOcean_cryptic<-betadiver(Occ_PacificOcean_cryptic_Checklist,"w")
betadiv_PacificOcean_cryptic<-as.numeric(betadiv_PacificOcean_cryptic)
betadiv_PacificOcean_cryptic<-data.frame(betadiv_PacificOcean_cryptic)

write_csv(betadiv_PacificOcean_cryptic, "~/Desktop/Summer research/2021/betadiv_PacificOcean_cryptic.csv")

#Atlantic Ocean cryptic species beta diversity
Occ_AtlanticOcean_cryptic=NULL
for(x in specieslist$SpeciesName){
  Occ_AtlanticOcean_cryptic<-rbind(Occ_AtlanticOcean_cryptic,filter(Occ_AtlanticOcean, scientificName==x))
}

dggs<- dgconstruct(spacing=1000, metric=FALSE, resround='down')

Occ_AtlanticOcean_cryptic$cell <- dgtransform(dggs,Occ_AtlanticOcean_cryptic$decimalLatitude,Occ_AtlanticOcean_cryptic$decimalLongitude)

Occ_AtlanticOcean_cryptic_Checklist<-table(Occ_AtlanticOcean_cryptic$cell, Occ_AtlanticOcean_cryptic$scientificName)

betadiv_AtlanticOcean_cryptic<-betadiver(Occ_AtlanticOcean_cryptic_Checklist,"w")
betadiv_AtlanticOcean_cryptic<-as.numeric(betadiv_AtlanticOcean_cryptic)
betadiv_AtlanticOcean_cryptic<-data.frame(betadiv_AtlanticOcean_cryptic)

write_csv(betadiv_AtlanticOcean_cryptic, "~/Desktop/Summer research/2021/betadiv_AtlanticOcean_cryptic.csv")

#Make a table of beta diversity and standard deviation
betadiv_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_ArcticOcean.csv")
betadiv_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_AtlanticOcean.csv")
betadiv_IndianOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_IndianOcean.csv")
betadiv_PacificOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_PacificOcean.csv")
betadiv_SouthOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_SouthOcean.csv")

betadiversity<- c(mean(as.matrix(betadiv_ArcticOcean)),mean(as.matrix(betadiv_AtlanticOcean)),mean(as.matrix(betadiv_IndianOcean)),
                  mean(as.matrix(betadiv_PacificOcean)),mean(as.matrix(betadiv_SouthOcean)))
stdev<- c(sd(as.matrix(betadiv_ArcticOcean)),sd(as.matrix(betadiv_AtlanticOcean)),sd(as.matrix(betadiv_IndianOcean)),
         sd(as.matrix(betadiv_PacificOcean)),sd(as.matrix(betadiv_SouthOcean)))


betadiv_ArcticOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_ArcticOcean_cryptic.csv")
betadiv_AtlanticOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_AtlanticOcean_cryptic.csv")
betadiv_IndianOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_IndianOcean_cryptic.csv")
betadiv_PacificOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_PacificOcean_cryptic.csv")
betadiv_SouthOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_SouthOcean_cryptic.csv")

betadiversity_cryptic<- c(mean(as.matrix(betadiv_ArcticOcean_cryptic)),mean(as.matrix(betadiv_AtlanticOcean_cryptic)),mean(as.matrix(betadiv_IndianOcean_cryptic)),
                  mean(as.matrix(betadiv_PacificOcean_cryptic)),mean(as.matrix(betadiv_SouthOcean_cryptic)))
stdev_cryptic<- c(sd(as.matrix(betadiv_ArcticOcean_cryptic)),sd(as.matrix(betadiv_AtlanticOcean_cryptic)),sd(as.matrix(betadiv_IndianOcean_cryptic)),
          sd(as.matrix(betadiv_PacificOcean_cryptic)),sd(as.matrix(betadiv_SouthOcean_cryptic)))

table_betadiv<- data.frame( collumnname, betadiversity, stdev, betadiversity_cryptic, stdev_cryptic)
colnames(table_betadiv)<- c("Ocean","beta diversity","standard deviation","beta diversity cryptic","standard deviation cryptic")
write.csv(table_betadiv,"~/Desktop/Summer research/2021/table_betadiv.csv")


#Make charts 

ggplot(table_betadiv) +
  geom_bar( aes(x=collumnname, y=betadiversity), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=collumnname, ymin=betadiversity-stdev, ymax=betadiversity+stdev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  xlab("Oceans")+
  ylab("beta diversity")

ggplot(table_betadiv) +
  geom_bar( aes(x=collumnname, y=betadiversity_cryptic), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=collumnname, ymin=betadiversity_cryptic-stdev_cryptic, ymax=betadiversity_cryptic+stdev_cryptic), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  xlab("Oceans")+
  ylab("beta diversity just cryptic")

#T-test on beta diversity
#Just cryptic between oceans
betadiv_ArcticOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_ArcticOcean_cryptic.csv")
betadiv_AtlanticOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_AtlanticOcean_cryptic.csv")
betadiv_IndianOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_IndianOcean_cryptic.csv")
betadiv_PacificOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_PacificOcean_cryptic.csv")
betadiv_SouthOcean_cryptic <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_SouthOcean_cryptic.csv")

betadiv_ArcticOcean_cryptic$Ocean<-"Arctic"
betadiv_AtlanticOcean_cryptic$Ocean<-"Atlantic"
betadiv_IndianOcean_cryptic$Ocean<-"Indian"
betadiv_PacificOcean_cryptic$Ocean<-"Pacific"
betadiv_SouthOcean_cryptic$Ocean<-"South"
collumnname<-c("betadiv","Ocean")
colnames(betadiv_ArcticOcean_cryptic)<-collumnname
colnames(betadiv_AtlanticOcean_cryptic)<-collumnname
colnames(betadiv_IndianOcean_cryptic)<-collumnname
colnames(betadiv_PacificOcean_cryptic)<-collumnname
colnames(betadiv_SouthOcean_cryptic)<-collumnname


betadiv_Ocean_cryptic<-rbind(betadiv_ArcticOcean_cryptic,betadiv_AtlanticOcean_cryptic,betadiv_IndianOcean_cryptic,
                             betadiv_PacificOcean_cryptic,betadiv_SouthOcean_cryptic)


TukeyHSD(aov(betadiv_Ocean_cryptic$betadiv~betadiv_Ocean_cryptic$Ocean))

#All species between oceans
betadiv_ArcticOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_ArcticOcean.csv")
betadiv_AtlanticOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_AtlanticOcean.csv")
betadiv_IndianOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_IndianOcean.csv")
betadiv_PacificOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_PacificOcean.csv")
betadiv_SouthOcean <- read_csv("World's Oceans/Occ Oceans/beta diversity/betadiv_SouthOcean.csv")

betadiv_ArcticOcean$Ocean<-"Arctic"
betadiv_AtlanticOcean$Ocean<-"Atlantic"
betadiv_IndianOcean$Ocean<-"Indian"
betadiv_PacificOcean$Ocean<-"Pacific"
betadiv_SouthOcean$Ocean<-"South"
collumnname<-c("betadiv","Ocean")
colnames(betadiv_ArcticOcean)<-collumnname
colnames(betadiv_AtlanticOcean)<-collumnname
colnames(betadiv_IndianOcean)<-collumnname
colnames(betadiv_PacificOcean)<-collumnname
colnames(betadiv_SouthOcean)<-collumnname


betadiv_Ocean<-rbind(betadiv_ArcticOcean,betadiv_AtlanticOcean,betadiv_IndianOcean,
                             betadiv_PacificOcean,betadiv_SouthOcean)


TukeyHSD(aov(betadiv_Ocean$betadiv~betadiv_Ocean$Ocean))

#Compare all species to just cryptic
#Arctic
betadiv_ArcticOcean_cryptic$data<-"cryptic"
betadiv_ArcticOcean$data<-"all"
betadiv_Arctic<-rbind(betadiv_ArcticOcean_cryptic,betadiv_ArcticOcean)
summary(aov(betadiv_Arctic$betadiv~betadiv_Arctic$data))

#Atlantic
betadiv_AtlanticOcean_cryptic$data<-"cryptic"
betadiv_AtlanticOcean$data<-"all"
betadiv_Atlantic<-rbind(betadiv_AtlanticOcean_cryptic,betadiv_AtlanticOcean)
summary(aov(betadiv_Atlantic$betadiv~betadiv_Atlantic$data))

#Indian
betadiv_IndianOcean_cryptic$data<-"cryptic"
betadiv_IndianOcean$data<-"all"
betadiv_Indian<-rbind(betadiv_IndianOcean_cryptic,betadiv_IndianOcean)
summary(aov(betadiv_Indian$betadiv~betadiv_Indian$data))

#Pacific
betadiv_PacificOcean_cryptic$data<-"cryptic"
betadiv_PacificOcean$data<-"all"
betadiv_Pacific<-rbind(betadiv_PacificOcean_cryptic,betadiv_PacificOcean)
summary(aov(betadiv_Pacific$betadiv~betadiv_Pacific$data))

#South
betadiv_SouthOcean_cryptic$data<-"cryptic"
betadiv_SouthOcean$data<-"all"
betadiv_South<-rbind(betadiv_SouthOcean_cryptic,betadiv_SouthOcean)
summary(aov(betadiv_South$betadiv~betadiv_South$data))


