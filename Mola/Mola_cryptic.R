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


#Extract COI sequences and location
Seq_Mola <- bold_seqspec(taxon="Mola", marker="COI-5P")
Seq_Mola<-Seq_Mola%>%filter(!is.na(lat))

#Putting sequences into regions
#Specify the regions
Seq_Mola$country[5]<-"Northwest_Pacific"
Seq_Mola$country[which(Seq_Mola$country == "Pacific Ocean")] <- "New_Zealand"
Seq_Mola$country[which(Seq_Mola$country == "United States")] <- "US_West_Coast"

#Making a map to visualize
world <- map_data("world")
#Region map
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_point(data = filter(Seq_Mola,country == "New Zealand"), aes(x = lon, y = lat), color='Red')+
  geom_point(data = filter(Seq_Mola,country == "US West Coast"), aes(x = lon, y = lat), color='Blue')+
  geom_point(data = filter(Seq_Mola,country == "Australia"), aes(x = lon, y = lat), color='Orange')+
  geom_point(data = filter(Seq_Mola,country == "Sweden"), aes(x = lon, y = lat), color='Green')+
  geom_point(data = filter(Seq_Mola,country == "Portugal"), aes(x = lon, y = lat), color='Purple')+
  geom_point(data = filter(Seq_Mola,country == "Northwest Pacific"), aes(x = lon, y = lat), color='Black')

#Species map
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_point(data = filter(Seq_Mola,species_name == "Mola mola"), aes(x = lon, y = lat), color='Blue')+
  geom_point(data = filter(Seq_Mola,species_name == "Mola sp. A MN-2017"), aes(x = lon, y = lat), color='Red')+
  geom_point(data = filter(Seq_Mola,species_name == "Mola sp. B"), aes(x = lon, y = lat), color='Green')+
  geom_point(data = filter(Seq_Mola,species_name == "Mola tecta"), aes(x = lon, y = lat), color='Purple')+
  geom_point(data = filter(Seq_Mola,species_name == ""), aes(x = lon, y = lat), color='Black')


write.csv(Seq_Mola,"~/Desktop/Summer research/2021/Github/Mola/Seq_Mola_COI.csv")
write.fasta(sequences = as.list(Seq_Mola$nucleotides), names = as.list(paste(Seq_Mola$processid,Seq_Mola$country,sep = "_")),file.out = "Seq_Mola.fasta" )

#Making phyologenies
Mola_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Mola/Seq_Mola_COI_aligned.aln", format = "clustal")
Mola_dist <- dist.dna(Mola_Aligned)

#UPGMA
Mola_UPGMA <- upgma(Mola_dist)

#Color tip based on location
tipcol <- rep('black', length(Mola_UPGMA$tip.label))
Area<-c("Australia" ,"New_Zealand" ,"Northwest_Pacifi","Portugal","Sweden","US_West_Coast")
Color<-c("blue","red","green","purple","orange","black")
for (i in 1:6){
  tipcol[grep(Area[i], Mola_UPGMA$tip.label)] <- Color[i]
}
plot(Mola_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 1)

#color tip based on identified species
tipcol_species <- rep('black', length(Mola_UPGMA$tip.label))
Mola_mola<-Seq_Mola$processid[which(Seq_Mola$species_name == "Mola mola")]
Mola_spA<-Seq_Mola$processid[which(Seq_Mola$species_name == "Mola sp. A MN-2017")]
Mola_spB<-Seq_Mola$processid[which(Seq_Mola$species_name == "Mola sp. B")]
Mola_tecta<-Seq_Mola$processid[which(Seq_Mola$species_name == "Mola tecta")]

for (i in 1:8){tipcol_species[grep(Mola_mola[i], Mola_UPGMA$tip.label)] <- "blue"}
for (i in 1:3){tipcol_species[grep(Mola_spA[i], Mola_UPGMA$tip.label)] <- "red"}
tipcol_species[grep(Mola_spB, Mola_UPGMA$tip.label)] <- "green"
tipcol_species[grep(Mola_tecta, Mola_UPGMA$tip.label)] <- "purple"

plot(Mola_UPGMA, main="UPGMA",tip.color=tipcol_species,use.edge.length = FALSE,cex = 1)

#Neighbor joining
Mola_NJ  <- NJ(Mola_dist)
plot(Mola_NJ, main = "Neighbor Joining",tip.color=tipcol,use.edge.length = FALSE, cex = 1)
plot(Mola_NJ, main = "Neighbor Joining",tip.color=tipcol_species,use.edge.length = FALSE, cex = 1)

#Maximum likelihood
Mola_phydat<-phyDat(Mola_Aligned, type = "DNA", levels = NULL)
fit_UPGMA <- pml(Mola_UPGMA,Mola_phydat)
fit_NJ <- pml(Mola_NJ, Mola_phydat)
print(fit_UPGMA)
print(fit_NJ)
fitJC  <- optim.pml(fit_UPGMA, TRUE)

tipcol_species_bs <- rep('black', length(fitJC$tree$tip.label))

for (i in 1:8){tipcol_species_bs[grep(Mola_mola[i], fitJC$tree$tip.label)] <- "blue"}
for (i in 1:3){tipcol_species_bs[grep(Mola_spA[i], fitJC$tree$tip.label)] <- "red"}
tipcol_species_bs[grep(Mola_spB, fitJC$tree$tip.label)] <- "green"
tipcol_species_bs[grep(Mola_tecta, fitJC$tree$tip.label)] <- "purple"

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol_species_bs, cex = 1,use.edge.length = FALSE)


#Find the differences between clades in the Phylogeny
x<-as.matrix(Mola_dist)
l <- data.frame(t(combn(rownames(x),2)), as.numeric(Mola_dist))
names(l) <- c("c1", "c2", "distance")
mola_dif<-as.data.frame(l)
write_csv(mola_dif, "~/Desktop/Summer research/2021/Github/Mola/Mola_COI_dist.csv")

mola_spA_c1<-c(grep("GBMIN97867",mola_dif$c1),grep("GBMIN97866",mola_dif$c1))
mola_spA_c2<-c(grep("GBMIN97867",mola_dif$c2),grep("GBMIN97866",mola_dif$c2))

mola_spB_c1<-c(grep("GBMIN97865",mola_dif$c1),grep("GBMIN97864",mola_dif$c1),grep("GBMIN127327",mola_dif$c1),grep("GBMIN133155",mola_dif$c1),grep("GBMIN127325",mola_dif$c1),
               grep("GBMIN127326",mola_dif$c1),grep("ANGBF29549",mola_dif$c1),grep("AMS174",mola_dif$c1),grep("FOAO2277",mola_dif$c1))
mola_spB_c2<-c(grep("GBMIN97865",mola_dif$c2),grep("GBMIN97864",mola_dif$c2),grep("GBMIN127327",mola_dif$c2),grep("GBMIN133155",mola_dif$c2),grep("GBMIN127325",mola_dif$c2),
               grep("GBMIN127326",mola_dif$c2),grep("ANGBF29549",mola_dif$c2),grep("AMS174",mola_dif$c2),grep("FOAO2277",mola_dif$c2))

mola_spC_c1<-c(grep("ANGBF46642",mola_dif$c1),grep("ANGBF46643",mola_dif$c1),grep("GBMIN127324",mola_dif$c1),grep("GBMIN122614",mola_dif$c1))
mola_spC_c2<-c(grep("ANGBF46642",mola_dif$c2),grep("ANGBF46643",mola_dif$c2),grep("GBMIN127324",mola_dif$c2),grep("GBMIN122614",mola_dif$c2))

mola_spD_c1<-c(grep("US_",mola_dif$c1),grep("ANGBF46644",mola_dif$c1),grep("FMVIC396",mola_dif$c1))
mola_spD_c2<-c(grep("US_",mola_dif$c2),grep("ANGBF46644",mola_dif$c2),grep("FMVIC396",mola_dif$c2))

mola_spE_c1<-c(grep("Portugal",mola_dif$c1),grep("Sweden",mola_dif$c1))
mola_spE_c2<-c(grep("Portugal",mola_dif$c2),grep("Sweden",mola_dif$c2))

dif_spA_B<-mean(mola_dif$distance[c(intersect(mola_spA_c1,mola_spB_c2), intersect(mola_spA_c2,mola_spB_c1))])
dif_spA_C<-mean(mola_dif$distance[c(intersect(mola_spA_c1,mola_spC_c2), intersect(mola_spA_c2,mola_spC_c1))])
dif_spA_D<-mean(mola_dif$distance[c(intersect(mola_spA_c1,mola_spD_c2), intersect(mola_spA_c2,mola_spD_c1))])
dif_spA_E<-mean(mola_dif$distance[c(intersect(mola_spA_c1,mola_spE_c2), intersect(mola_spA_c2,mola_spE_c1))])
dif_spB_C<-mean(mola_dif$distance[c(intersect(mola_spB_c1,mola_spC_c2), intersect(mola_spB_c2,mola_spC_c1))])
dif_spB_D<-mean(mola_dif$distance[c(intersect(mola_spB_c1,mola_spD_c2), intersect(mola_spB_c2,mola_spD_c1))])
dif_spB_E<-mean(mola_dif$distance[c(intersect(mola_spB_c1,mola_spE_c2), intersect(mola_spB_c2,mola_spE_c1))])
dif_spC_D<-mean(mola_dif$distance[c(intersect(mola_spC_c1,mola_spD_c2), intersect(mola_spC_c2,mola_spD_c1))])
dif_spC_E<-mean(mola_dif$distance[c(intersect(mola_spC_c1,mola_spE_c2), intersect(mola_spC_c2,mola_spE_c1))])
dif_spD_E<-mean(mola_dif$distance[c(intersect(mola_spD_c1,mola_spE_c2), intersect(mola_spD_c2,mola_spE_c1))])

dif_Mola<-data.frame(dif_spA_B,dif_spA_C,dif_spA_D,dif_spA_E,dif_spB_C,dif_spB_D,dif_spB_E,dif_spC_D,dif_spC_E,dif_spD_E)
write_csv(dif_Mola, "~/Desktop/Summer research/2021/Github/Mola/dif_Mola_COI.csv")


#Mola species D-loop gene from Ichthyol Res (2009)
Mola_Dloop <- read_csv("~/Desktop/Summer research/2021/Github/Mola/Mola_Dloop.csv")
seq_Mola_Dloop<-read.GenBank(Mola_Dloop$`Accession number`)
write.fasta(sequences = as.character.DNAbin(seq_Mola_Dloop), names = as.list(paste(Mola_Dloop$`Accession number`,Mola_Dloop$Location,sep = "_")),file.out = "Seq_Mola_Dloop.fasta" )

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_point(data = filter(Mola_Dloop,Location == "Pacific_northeastern_Japan"), aes(x = lon, y = lat), color='Red')+
  geom_point(data = filter(Mola_Dloop,Location == "Ogasawara (Bonin) Islands, Japan"), aes(x = lon, y = lat), color='Blue')+
  geom_point(data = filter(Mola_Dloop,Location == "Pacific_southern_Japan"), aes(x = lon, y = lat), color='Orange')+
  geom_point(data = filter(Mola_Dloop,Location == "Ryukyu_Islands, Japan"), aes(x = lon, y = lat), color='Green')+
  geom_point(data = filter(Mola_Dloop,Location == "Sea_of_Japan, Japan"), aes(x = lon, y = lat), color='Purple')+
  geom_point(data = filter(Mola_Dloop,Location == "North_Pacific Ocean (outside Japan)"), aes(x = lon, y = lat), color='Black')+
  geom_point(data = filter(Mola_Dloop,Location == "Taiwan"), aes(x = lon, y = lat), color='Pink')+
  geom_point(data = filter(Mola_Dloop,Location == "Denmark"), aes(x = lon, y = lat), color='Yellow')+
  geom_point(data = filter(Mola_Dloop,Location == "Australia (Tasman Sea)"), aes(x = lon, y = lat), color='Brown')+
  geom_point(data = filter(Mola_Dloop,Location == "UK"), aes(x = lon, y = lat), color='White')


ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "#dddddd") +
  geom_point(data = filter(Mola_Dloop,Group == "A"), aes(x = lon, y = lat), color='Red')+
  geom_point(data = filter(Mola_Dloop,Group == "B"), aes(x = lon, y = lat), color='Blue')

#Making phyologenies
Mola_Dloop_Aligned<-read.dna("~/Desktop/Summer research/2021/Github/Mola/Seq_Mola_Dloop_Aligned.aln", format = "clustal")
Mola_Dloop_dist <- dist.dna(Mola_Dloop_Aligned)

#UPGMA
Mola_Dloop_UPGMA <- upgma(Mola_Dloop_dist)

#Color tip based on location
tipcol <- rep('black', length(Mola_Dloop_UPGMA$tip.label))
Area<-c("Pacific_northeastern" ,"Ogasawara" ,"Pacific_southern","Ryukyu","Sea_of_Japan",
        "North_Pacific","Taiwan","Denmark","Australia","UK")
Color<-c("red","blue","orange","green","purple","black","pink","yellow","brown","white")
for (i in 1:10){
  tipcol[grep(Area[i], Mola_Dloop_UPGMA$tip.label)] <- Color[i]
}

#UPGMA
plot(Mola_Dloop_UPGMA, main="UPGMA",tip.color=tipcol,use.edge.length = FALSE,cex = 0.3)

#Neighbor joining
Mola_Dloop_NJ  <- NJ(Mola_Dloop_dist)
plot(Mola_Dloop_NJ, main = "Neighbor Joining",tip.color=tipcol, cex = 0.3)


#Maximum likelihood
Mola_Dloop_phydat<-phyDat(Mola_Dloop_Aligned, type = "DNA", levels = NULL)
fit_NJ <- pml(Mola_Dloop_NJ, Mola_Dloop_phydat)
print(fit_NJ)
fitJC  <- optim.pml(fit_NJ, TRUE)

tipcol_bs <- rep('black', length(fitJC$tree$tip.label))

for (i in 1:10){
  tipcol_bs[grep(Area[i], fitJC$tree$tip.label)] <- Color[i]
}

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p",tip.color=tipcol_bs, cex = 0.3)


#Color tip based on groups
tipcol_group <- rep('black', length(Mola_Dloop_UPGMA$tip.label))
groupA<-Mola_Dloop$`Accession number`[which(Mola_Dloop$Group == "A")]
for (i in 1:20){
  tipcol_group[grep(groupA[i], Mola_Dloop_UPGMA$tip.label)] <- "red"
}

groupB<-Mola_Dloop$`Accession number`[which(Mola_Dloop$Group == "B")]
for (i in 1:101){
  tipcol_group[grep(groupB[i], Mola_Dloop_UPGMA$tip.label)] <- "blue"
}


#UPGMA
plot(Mola_Dloop_UPGMA, main="UPGMA",tip.color=tipcol_group,cex = 0.3)


#Neighbor joining
plot(Mola_Dloop_NJ, main = "Neighbor Joining",tip.color=tipcol_group, cex = 0.3)

