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



#BOLD- extract gene seq and location COI
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

write.csv(Seq_Aureliaaurita,"~/Desktop/Summer research/2021/Seq_Aurelia_COI.csv")

#Phylogenies
write.fasta(sequences = as.list(Seq_Aureliaaurita$nucleotides), names = as.list(Seq_Aureliaaurita$processid),file.out = "Seq_Aurelia_COI.fasta" ) 
Aureliaaurita_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Seq_Aurelia_COI_Aligned.aln", format = "clustal")
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



#All Aurelia cryptic species COI_BOLD
Seq_Aurelia_COI <- bold_seqspec(taxon="Aurelia", marker="COI-5P")
#Seq_Aurelia_COI<-Seq_Aurelia_COI%>%filter(!is.na(lat))
Seq_Aurelia_COI<-Seq_Aurelia_COI%>%filter(species_name=="Aurelia aurita"|species_name=="Aurelia coerulea"|
                                            species_name=="Aurelia limbata"|species_name=="Aurelia labiata")

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
write.fasta(sequences = as.list(Seq_Aurelia_COI$nucleotides), names = as.list(Seq_Aurelia_COI$processid),file.out = "Seq_Aurelia_COI_species.fasta" ) 
Aurelia_COI_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Sequences/Seq_Aurelia_COI_Aligned.aln", format = "clustal")
Aurelia_COI_Aligned_phydat<-phyDat(Aurelia_COI_Aligned, type = "DNA", levels = NULL)

Aurelia_COI_Aligned_dist <- dist.ml(Aurelia_COI_Aligned_phydat)

#UPGMA
Aurelia_COI_Aligned_UPGMA <- upgma(Aurelia_COI_Aligned_dist)
tipcol <- rep('black', length(Aurelia_COI_Aligned_UPGMA$tip.label))
Species<-c("Aurelia aurita","Aurelia coerulea","Aurelia labiata","Aurelia limbata")
Color<-c("yellow","blue","red","green")
for (i in 1:4){
  tipcol[grep(Species[i], Seq_Aurelia_COI$species_name)] <- Color[i]
}

plot(Aurelia_COI_Aligned_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.25)
nodelabels(cex=0.2)

#Neighbor joining
Aurelia_COI_Aligned_NJ  <- NJ(Aurelia_COI_Aligned_dist)
plot(Aurelia_COI_Aligned_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 0.5)


#ITS1 gene according to the paper(Korsun, 2012)
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


#Phylogenies cryptic ITS1 from genbank
Aurelia_cryptic_ITS1<-read.dna("~/Desktop/Summer research/2021/Github/Aurelia Aurita/Sequences/Aurelia_ITS1_Genbank.aln", format = "clustal")
Aurelia_cryptic_ITS1_phydat<-phyDat(Aurelia_cryptic_ITS1, type = "DNA", levels = NULL)

Aurelia_cryptic_ITS1_dist <- dist.ml(Aurelia_cryptic_ITS1_phydat)

#UPGMA
Aurelia_cryptic_ITS1_UPGMA <- upgma(Aurelia_cryptic_ITS1_dist)
tipcol <- rep('black', length(Aurelia_cryptic_ITS1_UPGMA$tip.label))
Species<-c("aurita","labiata","limbata")
Color<-c("blue","red","green")
for (i in 1:3){
  tipcol[grep(Species[i],Aurelia_cryptic_ITS1_UPGMA$tip.label )] <- Color[i]
}

plot(Aurelia_cryptic_ITS1_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.5)
nodelabels(cex=0.2)

#Neighbor joining
Aurelia_cryptic_ITS1_NJ  <- NJ(Aurelia_cryptic_ITS1_dist)
plot(Aurelia_cryptic_ITS1_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 1)

#Calculating their differences
x<-as.matrix(Aurelia_cryptic_ITS1_dist)
m <- data.frame(t(combn(rownames(x),2)), as.numeric(Aurelia_cryptic_ITS1_dist))
names(m) <- c("c1", "c2", "distance")
n<-as.data.frame(m)

aurita_labiata<-c(intersect(grep("aurita",n$c1), grep("labiata",n$c2)),intersect(grep("labiata",n$c1), grep("aurita",n$c2)))
aurita_limbata<-c(intersect(grep("aurita",n$c1), grep("limbata",n$c2)),intersect(grep("limbata",n$c1), grep("aurita",n$c2)))
labiata_limbata<-c(intersect(grep("labiata",n$c1), grep("limbata",n$c2)),intersect(grep("limbata",n$c1), grep("labiata",n$c2)))

dif_aurita_labiata<-mean(n$distance[aurita_labiata])
dif_aurita_limbata<-mean(n$distance[aurita_limbata])
dif_labiata_limbata<-mean(n$distance[labiata_limbata])

aurita1_c1<-c(grep("KC767902.1",n$c1),grep("AY935206.1",n$c1),grep("AY935205.1",n$c1) , grep("AY319849.1",n$c1))
aurita1_c2<-c(grep("KC767902.1",n$c2),grep("AY935206.1",n$c2),grep("AY935205.1",n$c2) , grep("AY319849.1",n$c2))

aurita2_c1<-c(grep("AF461405.1",n$c1),grep("EF010536.1",n$c1),grep("KC767900.1",n$c1))
aurita2_c2<-c(grep("AF461405.1",n$c2),grep("EF010536.1",n$c2),grep("KC767900.1",n$c2))

aurita1_aurita2<- c(intersect(aurita1_c1,aurita2_c2), intersect(aurita1_c2,aurita2_c1))
dif_aurita1_aurita2<- mean(n$distance[aurita1_aurita2])

aurita3_c1<-c(grep("AY319848.1",n$c1),grep("KC767901.1",n$c1))
aurita3_c2<-c(grep("AY319848.1",n$c2),grep("KC767901.1",n$c2))

aurita1_aurita3<- c(intersect(aurita1_c1,aurita3_c2), intersect(aurita1_c2,aurita3_c1))
dif_aurita1_aurita3<- mean(n$distance[aurita1_aurita3])

aurita2_aurita3<- c(intersect(aurita2_c1,aurita3_c2), intersect(aurita2_c2,aurita3_c1))
dif_aurita2_aurita3<- mean(n$distance[aurita2_aurita3])

dif_ITS1<-data.frame(dif_aurita_labiata,dif_aurita_limbata,dif_labiata_limbata)
names(dif_ITS1) <- c("aurita and labiata", "aurita and limbata", "labiata and limbata")
write_csv(dif_ITS1, "~/Desktop/Summer research/2021/Github/Aurelia Aurita/dif_Aurelia_ITS1_Genbank.csv")

#Heatmap
heatmap(x,scale = "none",cexRow = 0.5,cexCol = 0.40,symm = TRUE)



