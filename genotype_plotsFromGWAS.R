#2023-11-26 gwas genotype plots

library(ggplot2)
library(dplyr)
library(ggstatsplot)
library(RColorBrewer)
library(pals)
library(ggforce)
library(ggpubr)
setwd("../Manhattan_plots")


#load significant SNP GTs
sigsnpGT <- list.files("Manhattan_plots",pattern="GTs.txt",full.names = T)

GT_DF <- lapply(sigsnpGT, function(x){
  read.delim(x,header=F)})

#load sample order
vcf_sample_order <- read.table("vcf_sample_order.txt", quote="\"", comment.char="")

significant_snps_metadata <- read.delim("significant_snps_metadata.txt")
significant_snps_metadata$snpids <- paste(significant_snps_metadata$chr, significant_snps_metadata$ps,sep="_")

for (i in 1:length(GT_DF)){
  colnames(GT_DF[[i]]) <- c("chr","pos",vcf_sample_order$V1)
  snpids <- paste(GT_DF[[i]]$chr,GT_DF[[i]]$pos,sep="_")
  GT_DF[[i]] <- as.data.frame(t(GT_DF[[i]][,-c(1,2)]))
  colnames(GT_DF[[i]]) <- snpids
  GT_DF[[i]]$CRC <- rownames(GT_DF[[i]])
  GT_DF[[i]][GT_DF[[i]]=="0|1" | GT_DF[[i]]=="1|0"] <- "0/1"
  GT_DF[[i]][GT_DF[[i]]=="1|1"] <- "1/1"
  
  
  
  }
names(GT_DF) <-  gsub("Manhattan_plots/","",sigsnpGT)
names(GT_DF)<- gsub("_GTs.txt","",names(GT_DF) )



#load raw phenotype data and meta data 
metadata_CVC <- read.delim("2024-03-29_meta_data.txt")
metadata_CVC$Location = gsub("_","-",metadata_CVC$Location)


#load packline data
outrep <- read.delim("../PhenotypeDataProcessing/two_year_packline_phenotypes_mean_per_tree.txt")
packline <- unique(merge(outrep, metadata_CVC, by="CRC"))


dest <- read.delim("../PhenotypeDataProcessing/TWO_years_destructive_fruit_quality_phenotypes_per_tree.txt")

dest <- merge(dest, metadata_CVC, by="CRC")

#mark parents of fairchild
packline$CRC <- ifelse(packline$CRC =="279","Clementine",packline$CRC)
packline$CRC <- ifelse(packline$CRC =="3874","Orlando",packline$CRC)

dest$CRC <- ifelse(dest$CRC =="279","Clementine",dest$CRC)
dest$CRC <- ifelse(dest$CRC =="3874","Orlando",dest$CRC)


#add GT data
packtrait_and_GT <- list()
desttrait_andGT <- list()
traits <- gsub("Manhattan_plots/","",sigsnpGT)
traits <- gsub("_GTs.txt","",traits)

for (i in 1:length(traits)){
  packtrait_and_GT[[i]] <- merge(packline,GT_DF[[i]], GT_DF[[i]], by.x="CRC",by.y="CRC")
  desttrait_andGT[[i]] <- merge(dest,GT_DF[[i]], GT_DF[[i]], by.x="CRC",by.y="CRC")
}

names(packtrait_and_GT) <- traits
names(desttrait_andGT) <- traits



pdf("MajDiamgenotype_plot.pdf")
ggplot(packtrait_and_GT[[11]] %>% filter(chr1_2820321 != "./."), aes(x=chr1_2820321,y=MajorDiameter.mean, color = chr1_2820321))+geom_boxplot(show.legend = F) +
  geom_jitter(width=0.1,alpha=0.5,size=2)+

  theme_classic() + 
  xlab("Genotype") + facet_grid(~Year)+ scale_color_manual(values=c("slategrey","slateblue1","slateblue4"))



#clean up market types by collapsing them
packtrait_and_GT[[11]]$Type.x <- gsub("-Bittersweet","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Blood oranges","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Valencia","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Common","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Navel","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Pigmented","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-White","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("-Valencia","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub(" \\(Graft\\)","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("GFT-Hybrid","GFT",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("GFT-Hybrid","GFT",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("Citrus species","Mand-Hybrid",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("TRIF-","",packtrait_and_GT[[11]]$Type.x)
packtrait_and_GT[[11]]$Type.x <- gsub("PUMM-Hybrid","GFT",packtrait_and_GT[[11]]$Type.x)


ggplot(packtrait_and_GT[[11]] %>% filter(chr1_2820321 != "./."), aes(x=chr1_2820321,y=MajorDiameter.mean))+geom_boxplot(show.legend = F,outlier.alpha = 0) +
  geom_jitter(width=0.1,alpha=1,size=2,aes(color=Type.x))+
  theme_classic() +
  xlab("Genotype") + facet_grid(~Year) + scale_colour_brewer(palette="Paired")

dev.off()


