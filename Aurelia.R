rm(list=ls())

library(dplyr)
library(robis)
library(readr)
library(ggplot2)
library(ape)
library(phangorn)
library(bold)
library(rgeos)
library(rgdal)
library(seqinr)

#OBIS Occ of Aurelia aurita
Occ_Aureliaaurita<-occurrence(scientificname = "Aurelia aurita")



#BOLD- extract gene seq and location
Seq_Aureliaaurita <- bold_seqspec(taxon="Aurelia Aurita", marker="COI-5P")
Seq_Aureliaaurita<-Seq_Aureliaaurita%>%filter(!is.na(lat))

#Put occurrences into polygons

World_Oceans <- read_csv("~/Desktop/Summer research/2021/Github/World's Oceans/Coordinates/World_Oceans.csv")

point.x=Seq_Aureliaaurita$lat
point.y=Seq_Aureliaaurita$lon

Seq_Aureliaaurita$poly<-"blank"
for (x in World_Oceans$X1){
  Polygon= readWKT(paste(World_Oceans$coor[x]))
  Poly <- ggplot2::fortify(Polygon)
  inpoly<-point.in.polygon(point.x, point.y, Poly$lat,
                           Poly$long, mode.checked=FALSE)
  #if the element in >0 in poly assign poly name to poly column 
  Seq_Aureliaaurita$poly[inpoly>0]<-World_Oceans$Regions[x]
}

Seq_Aureliaaurita$poly[c(9,10,12,17,47,51,57,82,88,114,118,128)]<-"North Sea"
Seq_Aureliaaurita$poly[c(1:4,30:33,70:73, 107,108, 137:139)]<-"Pacific, Northeast"
Seq_Aureliaaurita$poly[c(21:25,58:66,94:100,130:133,154)]<-"Bay of Fundy"
Seq_Aureliaaurita$poly[c(26:28,67,101)]<-"Labrador Sea"
Seq_Aureliaaurita$poly[135]<-"Skagerrak"
Seq_Aureliaaurita$poly[29]<-"Hudson Bay"

map_leaflet(Seq_Aureliaaurita)

Seq_Aureliaaurita$Area<-NULL
Seq_Aureliaaurita$Area[which(Seq_Aureliaaurita$poly == "Labrador Sea" | Seq_Aureliaaurita$poly == "Bay of Fundy"|
                               Seq_Aureliaaurita$poly == "Gulf of St. Lawrence")] <- "Northwest Atlantic"
Seq_Aureliaaurita$Area[which(Seq_Aureliaaurita$poly == "North Sea" | Seq_Aureliaaurita$poly == "Skagerrak")] <- "Northeast Atlantic"
Seq_Aureliaaurita$Area[which(Seq_Aureliaaurita$poly == "Pacific, Northeast" | Seq_Aureliaaurita$poly == "Pacific, Eastern Central")] <- "East Pacific"
table(Seq_Aureliaaurita$Area)

write.csv(Seq_Aureliaaurita,"~/Desktop/Summer research/2021/Seq_Aureliaaurita.csv")

#APE
Trial<-read.dna("~/Desktop/sequence.aln", format = "clustal")
Trial_phydat<-phyDat(Trial, type = "DNA", levels = NULL)
mt<-modelTest(Trial_phydat)
print(mt)
dna_dist <- dist.ml(Trial_phydat)
Trial_UPGMA <- upgma(dna_dist)
plot(Trial_UPGMA, main="UPGMA")
Trial_NJ  <- NJ(dna_dist)
plot(Trial_NJ, main = "Neighbor Joining")

#Phylogenies
write.fasta(sequences = as.list(Seq_Aureliaaurita$nucleotides), names = as.list(Seq_Aureliaaurita$processid),file.out = "Seq_Aureliaaurita.fasta" ) 
Aureliaaurita_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Aureliaaurita_Aligned.aln", format = "clustal")
Aureliaaurita_phydat<-phyDat(Aureliaaurita_Aligned, type = "DNA", levels = NULL)

Aureliaaurita_dist <- dist.ml(Aureliaaurita_phydat)

#UPGMA
Aureliaaurita_UPGMA <- upgma(Aureliaaurita_dist)
tipcol <- rep('black', length(Aureliaaurita_UPGMA$tip.label))
Area<-c("East Pacific" ,"Northeast Atlantic" ,"Northwest Atlantic")
Color<-c("blue","red","green")
for (i in 1:3){
  tipcol[grep(Area[i], Seq_Aureliaaurita$Area)] <- Color[i]
}

plot(Aureliaaurita_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.25)
nodelabels(cex=0.2)

#Neighbor joining
Aureliaaurita_NJ  <- NJ(Aureliaaurita_dist)
plot(Aureliaaurita_NJ, main = "Neighbor Joining","u",tip.color=tipcol,use.edge.length = FALSE, cex = 0.25)

#Maximum likelihood
fit_UPGMA <- pml(Aureliaaurita_UPGMA, Aureliaaurita_phydat)
fit_NJ <- pml(Aureliaaurita_NJ, Aureliaaurita_phydat)
print(fit_NJ)
fitJC  <- optim.pml(fit_NJ, TRUE)
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol, cex = 0.5,use.edge.length = FALSE)



#Edited tree
Aureliaaurita_Aligned_edited<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/aurelia_edited.fas", format = "fasta")
Aureliaaurita_edited_phydat<-phyDat(Aureliaaurita_Aligned_edited, type = "DNA", levels = NULL)

Aureliaaurita_dist_edited <- dist.ml(Aureliaaurita_edited_phydat)

#UPGMA
Aureliaaurita_edited_UPGMA <- upgma(Aureliaaurita_dist_edited)
tipcol <- rep('black', length(Aureliaaurita_edited_UPGMA$tip.label))
Area<-c("East Pacific" ,"Northeast Atlantic" ,"Northwest Atlantic")
Color<-c("blue","red","green")
for (i in 1:3){
  tipcol[grep(Area[i], Seq_Aureliaaurita$Area)] <- Color[i]
}

plot(Aureliaaurita_edited_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.25)
nodelabels(cex=0.2)

#Neighbor joining
Aureliaaurita_NJ  <- NJ(Aureliaaurita_dist_edited)
plot(Aureliaaurita_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 0.25)

#Maximum likelihood
fit_UPGMA <- pml(Aureliaaurita_UPGMA, Aureliaaurita_phydat)
fit_NJ <- pml(Aureliaaurita_NJ, Aureliaaurita_phydat)
print(fit_NJ)
fitJC  <- optim.pml(fit_NJ, TRUE)
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol, cex = 0.5,use.edge.length = FALSE)



#ITS1 gene
Seq_Aureliaaurita_ITS1_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Seq_Aureliaaurita_ITS1.aln", format = "clustal")
Aureliaaurita_ITS1_Aligned_phydat<-phyDat(Seq_Aureliaaurita_ITS1_Aligned, type = "DNA", levels = NULL)

Aureliaaurita_ITS1_dist <- dist.ml(Aureliaaurita_ITS1_Aligned_phydat)


#UPGMA
Aureliaaurita_ITS1_UPGMA <- upgma(Aureliaaurita_ITS1_dist)
tipcol <- rep('black', length(Aureliaaurita_ITS1_UPGMA$tip.label))
Area<-c("West Pacific","East Pacific" ,"Northeast Atlantic" ,"Northwest Atlantic", "Southwest Atlantic")
Color<-c("yellow","blue","red","green","purple")
Aureliaaurita_ITS1_area<-c("Northeast Atlantic", "Arctic","Northeast Atlantic","Northwest Atlantic","Northeast Atlantic","Northeast Atlantic",
                           "Northwest Atlantic","West Pacific","Southwest Atlantic","Southwest Atlantic","East Pacific",
                           "East Pacific","East Pacific","East Pacific","East Pacific","East Pacific","Northeast Atlantic","Northeast Atlantic",
                           "West Pacific","West Pacific","West Pacific")
for (i in 1:5){
  tipcol[grep(Area[i], Aureliaaurita_ITS1_area)] <- Color[i]
}
plot(Aureliaaurita_ITS1_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 1)
nodelabels(cex=0.2)

#Neighbor joining
Aureliaaurita_ITS1_NJ  <- NJ(Aureliaaurita_ITS1_dist)
plot(Aureliaaurita_ITS1_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 0.75)

#Maximum likelihood
fit_UPGMA <- pml(Aureliaaurita_ITS1_UPGMA, Aureliaaurita_ITS1_Aligned_phydat)
fit_NJ <- pml(Aureliaaurita_ITS1_NJ, Aureliaaurita_ITS1_Aligned_phydat)
print(fit_UPGMA)
fitJC  <- optim.pml(fit_NJ, TRUE)
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol, cex = 1,use.edge.length = FALSE)



#ITS1 gene according to the paper
Seq_Aurelia_ITS1_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Seq_Aurelia_ITS1_aligned.aln", format = "clustal")
Aurelia_ITS1_Aligned_phydat<-phyDat(Seq_Aurelia_ITS1_Aligned, type = "DNA", levels = NULL)

Aurelia_ITS1_dist <- dist.ml(Aurelia_ITS1_Aligned_phydat)


#UPGMA
Aurelia_ITS1_UPGMA <- upgma(Aurelia_ITS1_dist)
tipcol <- rep('black', length(Aurelia_ITS1_UPGMA$tip.label))
Area<-c("Narragansett Bay","White Sea" ,"North Sea" ,"Marmara Sea", "Black Sea","Caspian Sea")
Color<-c("yellow","blue","red","green","purple","brown")
Aurelia_ITS1_area<-c("Narragansett Bay", "Narragansett Bay","Narragansett Bay","White Sea","White Sea","White Sea",
                           "North Sea","North Sea","North Sea","North Sea","Marmara Sea",
                           "Marmara Sea","Marmara Sea","Black Sea","Black Sea","Black Sea","Caspian Sea","Caspian Sea",
                           "Caspian Sea","Caspian Sea","Caspian Sea","Caspian Sea")
for (i in 1:6){
  tipcol[grep(Area[i], Aurelia_ITS1_area)] <- Color[i]
}
plot(Aurelia_ITS1_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 1)nodelabels(cex=0.2)

#Neighbor joining
Aurelia_ITS1_NJ  <- NJ(Aurelia_ITS1_dist)
plot(Aurelia_ITS1_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 0.75)

#Maximum likelihood
fit_UPGMA <- pml(Aureliaaurita_ITS1_UPGMA, Aureliaaurita_ITS1_Aligned_phydat)
fit_NJ <- pml(Aurelia_ITS1_NJ, Aurelia_ITS1_Aligned_phydat)
print(fit_UPGMA)
fitJC  <- optim.pml(fit_NJ, TRUE)
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol, cex = 1,use.edge.length = FALSE)


#All Aurelia COI
Seq_Aurelia_COI <- bold_seqspec(taxon="Aurelia", marker="COI-5P")
Seq_Aurelia_COI<-Seq_Aurelia_COI%>%filter(!is.na(lat))

point.x=Seq_Aurelia_COI$lat
point.y=Seq_Aurelia_COI$lon


Seq_Aurelia_COI$poly<-"blank"
for (x in World_Oceans$X1){
  Polygon= readWKT(paste(World_Oceans$coor[x]))
  Poly <- ggplot2::fortify(Polygon)
  inpoly<-point.in.polygon(point.x, point.y, Poly$lat,
                           Poly$long, mode.checked=FALSE)
  #if the element in >0 in poly assign poly name to poly column 
  Seq_Aurelia_COI$poly[inpoly>0]<-World_Oceans$Regions[x]
}
Seq_Aurelia_COI<-Seq_Aurelia_COI%>%filter(poly!="blank")

#Phylogenies
write.fasta(sequences = as.list(Seq_Aurelia_COI$nucleotides), names = as.list(Seq_Aurelia_COI$processid),file.out = "Seq_Aurelia_COI.fasta" ) 
Aurelia_COI_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Seq_Aurelia_COI_Aligned.aln", format = "clustal")
Aurelia_COI_Aligned_phydat<-phyDat(Aurelia_COI_Aligned, type = "DNA", levels = NULL)

Aurelia_COI_Aligned_dist <- dist.ml(Aurelia_COI_Aligned_phydat)

#UPGMA
Aurelia_COI_Aligned_UPGMA <- upgma(Aurelia_COI_Aligned_dist)
tipcol <- rep('black', length(Aurelia_COI_Aligned_UPGMA$tip.label))
Area<-c("Arctic, Central West","Gulf of St. Lawrence","Labrador Sea","North Sea","Pacific, Eastern Central",
        "Pacific, Northeast","Pacific, Nothwest","Skagerrak")
Color<-c("yellow","blue","red","green","purple","brown","black","orange")
for (i in 1:8){
  tipcol[grep(Area[i], Seq_Aurelia_COI$poly)] <- Color[i]
}

plot(Aurelia_COI_Aligned_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.5)
nodelabels(cex=0.2)

#Neighbor joining
Aurelia_COI_Aligned_NJ  <- NJ(Aurelia_COI_Aligned_dist)
plot(Aurelia_COI_Aligned_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 0.5)






