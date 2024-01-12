library(ggplot2)
library(dplyr)
library(regioneR)
library(ggforce)
detach("package:plyr", unload = TRUE)

TEs <- read.delim("TEs.gff", header=FALSE, comment.char="#")

TEs$V9 <- gsub("\\?","",TEs$V9)
TEs$V9 <- gsub("MuLE","MULE",TEs$V9)
TEs$subClass <- gsub(".*?Class=(.*?)\\;.*", "\\1", TEs$V9)
TEs$Class <- gsub("(.*?)\\/.*", "\\1", TEs$subClass)


TEs$size <- TEs$V5 - TEs$V4


summary(TEs$size)

TEs <- TEs[(grep("chr", TEs$V1)),]

mean_size <- mean(TEs$size)
TE_size <- ggplot(TEs, aes(x = size)) + geom_histogram(binwidth = 10) + xlim(0,10000) + xlab("TE size (bp)") + geom_vline(aes(xintercept = mean_size, col = 'red')) + theme_bw() + theme(legend.position = "none") +facet_wrap_paginate(~Class, nrow =3 , ncol =3,page=1)

DNA_TE <- TEs %>% filter(Class =="DNA")
Retro <- TEs %>% filter(Class == "LINE" | Class == "LTR")

#do permutations to see if AS ACRs overlap TEs more than expected
print("Importing AS marks")

ATAC_AS <- read.delim("ATAC_AS_wID.bed", header=FALSE)

#load in AS and all Histone mark peaks
H3K27me3_AS <- read.delim("H3K27me3_AS_wID.bed", header=FALSE)

H3K56ac_AS <- read.delim("H3K56ac_AS_wID.bed", header=FALSE)

H3K4me3_AS <- read.delim("H3K4me3_AS_wID.bed", header=FALSE)

H3K36me3_AS <- read.delim("H3K36me3_AS_wID.bed", header=FALSE)


#all peaks
print("Importing peaks and ACRs")
ACRs <- read.delim("ACRs_final.bed", header=FALSE)

H3K4me3_peaks <- read.delim("H3K4me3_peaks.txt.bed", header=FALSE)

H3K36me3_peaks <- read.delim("H3K36me3_peaks.txt.bed", header=FALSE)

H3K56ac_peaks <- read.delim("H3K56ac_peaks.txt.bed", header=FALSE)

H3K27me3_peaks <- read.delim("H3K27me3_peaks.txt.bed", header=FALSE)


print("Importing Promoter files")

All_promoters <- read.delim("All_promoters.sort.bed", header=FALSE)

Non_ASE_geneiase_promoters <- read.delim("Non_ASE_geneiase_promoters.bed", header=FALSE)

Geneiase_onlyPromoters<- read.delim("Geneiase_onlyPromoters.bed", header=FALSE)

ASE_promoter_1<- read.delim("ASE_promoter_1.sort.bed", header=FALSE)

ASE_promoter_2<- read.delim("ASE_promoter_2.sort.bed", header=FALSE)

ASE_promoters_both <- read.delim("ASE_promoters_both.bed", header=FALSE)

####
print("Making gRanges objects")

#make gRange objects of peaks and genes and promoters

H3K4me3_AS <- filterChromosomes(toGRanges(H3K4me3_AS), chr.type = "canonical")
H3K36me3_AS <- filterChromosomes(toGRanges(H3K36me3_AS), chr.type = "canonical")
H3K56ac_AS <- filterChromosomes(toGRanges(H3K56ac_AS), chr.type = "canonical")
H3K27me3_AS <- filterChromosomes(toGRanges(H3K27me3_AS), chr.type = "canonical")
ATAC_AS <- filterChromosomes(toGRanges(ATAC_AS), chr.type = "canonical")

promoters <- filterChromosomes(toGRanges(All_promoters, chr.type = "canonical"))


geneiase_promoters <- filterChromosomes(toGRanges(Geneiase_onlyPromoters), chr.type = "canonical")
Non_ASE_geneiase_promoters <- filterChromosomes(toGRanges(Non_ASE_geneiase_promoters), chr.type = "canonical")

#make promoter object with promoters that show ase in both reps
ASE_prom_both <- merge(ASE_promoter_1, ASE_promoter_2, by=c("V5"))[,c(1)]
onlyASE <- ASE_promoters_both$V5 %in% ASE_prom_both
ASE_promoters_both <- unique(ASE_promoters_both[which(onlyASE),])



#Make Granges object of promoters
ASE_promoters_1 <- filterChromosomes(toGRanges(ASE_promoter_1), chr.type = "canonical")
ASE_promoters_2 <- filterChromosomes(toGRanges(ASE_promoter_2), chr.type = "canonical")
ASE_bothreps <- filterChromosomes(toGRanges(ASE_promoters_both), chr.type = "canonical")

#Make Granges objects of peaks
H3K4me3_peaks <- filterChromosomes(toGRanges(H3K4me3_peaks), chr.type = "canonical")

H3K27me3_peaks <- filterChromosomes(toGRanges(H3K27me3_peaks), chr.type = "canonical")

H3K56ac_peaks <- filterChromosomes(toGRanges(H3K56ac_peaks), chr.type = "canonical")

H3K36me3_peaks <- filterChromosomes(toGRanges(H3K36me3_peaks), chr.type = "canonical")

ATAC_peaks <- filterChromosomes(toGRanges(ACRs), chr.type = "canonical")

TEs_ <- filterChromosomes(toGRanges(TEs[,c(1,4,5,11)], "canonical"))

DNA_TEs <- filterChromosomes(toGRanges(DNA_TE[,c(1,4,5,11)], "canonical"))

Retro_TEs <- filterChromosomes(toGRanges(Retro[,c(1,4,5,11)], "canonical"))


#Do my phased SNPs overlap TEs?

Linked_read_SNPs_sort <- read.delim("Linked_read_SNPs_sort.bed", header=FALSE)

SNPs <- filterChromosomes(toGRanges(Linked_read_SNPs_sort))

numOverlaps(SNPs, TEs_, count.once = TRUE)
TEs_test <- permTest(A=ATAC_AS, B = TEs_, randomize.function = resampleRegions, universe = ATAC_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parallel = FALSE, ntimes = 1000)
plot(TEs_test)

numOverlaps(ATAC_AS, TEs_, count.once=TRUE)


H3K27me3_te_test <- permTest(A=H3K27me3_AS, B = TEs_, randomize.function = resampleRegions, universe = H3K27me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K27me3_te_test)

H3K4me3_te_test <- permTest(A=H3K4me3_AS, B = TEs_, randomize.function = resampleRegions, universe =H3K4me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K4me3_te_test)

H3K36me3_te_test <- permTest(A=H3K36me3_AS, B = TEs_, randomize.function = resampleRegions, universe = H3K36me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K36me3_te_test)

H3K56ac_te_test <- permTest(A=H3K56ac_AS, B = TEs_, randomize.function = resampleRegions, universe = H3K56ac_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K56ac_te_test)


###

TE_AS_perm_objects <- list(TEs_test,H3K4me3_te_test,H3K36me3_te_test,H3K56ac_te_test,H3K27me3_te_test)
Marks <- c("ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3")

TE_AS_perm_objects[[1]]$Mark <- "ATAC"
TE_AS_perm_objects[[2]]$Mark <- "H3K4me3"
TE_AS_perm_objects[[3]]$Mark <- "H3K36me3"
TE_AS_perm_objects[[4]]$Mark <- "H3K56ac"
TE_AS_perm_objects[[5]]$Mark <- "H3K27me3"

TE_AS_perm_values <- list()
TE_AS_observed_numOverlap <- list()

for (i in 1:5){
 TE_AS_perm_values[[i]] <- TE_AS_perm_objects[[i]][["numOverlaps"]][["permuted"]]
 TE_AS_perm_values[[i]] <- as.data.frame(TE_AS_perm_values[[i]])
 TE_AS_perm_values[[i]]$Type <- "Permuted"
 TE_AS_perm_values[[i]]$Mark <- Marks[i]
 colnames(TE_AS_perm_values[[i]]) <- c("Overlaps","Type","Mark")
 TE_AS_observed_numOverlap[[i]] <- TE_AS_perm_objects[[i]][["numOverlaps"]][["observed"]]
 TE_AS_observed_numOverlap[[i]] <- as.data.frame(TE_AS_observed_numOverlap[[i]])
 TE_AS_observed_numOverlap[[i]]$Type <- "Observed"
 TE_AS_observed_numOverlap[[i]]$Mark <- Marks[i]
 colnames(TE_AS_observed_numOverlap[[i]]) <- c("Overlaps","Type","Mark")

}

All_TE_AS_overlaps <- rbind(TE_AS_perm_values[[1]],TE_AS_perm_values[[2]],TE_AS_perm_values[[3]],TE_AS_perm_values[[4]],TE_AS_perm_values[[5]],TE_AS_observed_numOverlap[[1]],TE_AS_observed_numOverlap[[2]],TE_AS_observed_numOverlap[[3]],TE_AS_observed_numOverlap[[4]],TE_AS_observed_numOverlap[[5]])

ggplot(All_TE_AS_overlaps, aes(x = Mark, fill = Type, y = Overlaps)) + geom_boxplot(outlier.shape = NA) + theme_bw() + scale_fill_brewer(palette ="Set2" ) + ylab("Number of Unique Overlaps with TEs")


DNA_TEs_test <- permTest(A=ATAC_AS, B = DNA_TEs, randomize.function = resampleRegions, universe = ATAC_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parallel = FALSE, ntimes = 1000)
plot(DNA_TEs_test)

numOverlaps(ATAC_AS, DNA_TEs, count.once=TRUE)


DNA_H3K27me3_te_test <- permTest(A=H3K27me3_AS, B = DNA_TEs, randomize.function = resampleRegions, universe = H3K27me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(DNA_H3K27me3_te_test)

DNA_H3K4me3_te_test <- permTest(A=H3K4me3_AS, B = DNA_TEs, randomize.function = resampleRegions, universe =H3K4me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(DNA_H3K4me3_te_test)

DNA_H3K36me3_te_test <- permTest(A=H3K36me3_AS, B = DNA_TEs, randomize.function = resampleRegions, universe = H3K36me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(DNA_H3K36me3_te_test)

DNA_H3K56ac_te_test <- permTest(A=H3K56ac_AS, B = DNA_TEs, randomize.function = resampleRegions, universe = H3K56ac_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(DNA_H3K56ac_te_test)


###

DNA_TE_AS_perm_objects <- list(DNA_TEs_test,DNA_H3K4me3_te_test,DNA_H3K36me3_te_test,DNA_H3K56ac_te_test,DNA_H3K27me3_te_test)
Marks <- c("ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3")

DNA_TE_AS_perm_objects[[1]]$Mark <- "ATAC"
DNA_TE_AS_perm_objects[[2]]$Mark <- "H3K4me3"
DNA_TE_AS_perm_objects[[3]]$Mark <- "H3K36me3"
DNA_TE_AS_perm_objects[[4]]$Mark <- "H3K56ac"
DNA_TE_AS_perm_objects[[5]]$Mark <- "H3K27me3"

DNA_TE_AS_perm_values <- list()
DNA_TE_AS_observed_numOverlap <- list()

for (i in 1:5) {
 DNA_TE_AS_perm_values[[i]] <- DNA_TE_AS_perm_objects[[i]][["numOverlaps"]][["permuted"]]
 DNA_TE_AS_perm_values[[i]] <- as.data.frame(DNA_TE_AS_perm_values[[i]])
 DNA_TE_AS_perm_values[[i]]$Type <- "Permuted"
 DNA_TE_AS_perm_values[[i]]$Mark <- Marks[i]
 colnames(DNA_TE_AS_perm_values[[i]]) <- c("Overlaps","Type","Mark")
 DNA_TE_AS_observed_numOverlap[[i]] <- DNA_TE_AS_perm_objects[[i]][["numOverlaps"]][["observed"]]
 DNA_TE_AS_observed_numOverlap[[i]] <- as.data.frame(DNA_TE_AS_observed_numOverlap[[i]])
 DNA_TE_AS_observed_numOverlap[[i]]$Type <- "Observed"
 DNA_TE_AS_observed_numOverlap[[i]]$Mark <- Marks[i]
 colnames(DNA_TE_AS_observed_numOverlap[[i]]) <- c("Overlaps","Type","Mark")
}

DNA_All_TE_AS_overlaps <- rbind(DNA_TE_AS_perm_values[[1]],DNA_TE_AS_perm_values[[2]],DNA_TE_AS_perm_values[[3]],DNA_TE_AS_perm_values[[4]],DNA_TE_AS_perm_values[[5]],DNA_TE_AS_observed_numOverlap[[1]],DNA_TE_AS_observed_numOverlap[[2]],DNA_TE_AS_observed_numOverlap[[3]],DNA_TE_AS_observed_numOverlap[[4]],DNA_TE_AS_observed_numOverlap[[5]])

ggplot(DNA_All_TE_AS_overlaps %>% filter(Mark == "ATAC"), aes(x = Mark, fill = Type, y = Overlaps)) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(text = element_text(size=20)) + scale_fill_manual(values=c("slategrey","slateblue1")) + ylab("Number of Unique Overlaps with DNA_TEs")



Retro_TEs_test <- permTest(A=ATAC_AS, B = Retro_TEs, randomize.function = resampleRegions, universe = ATAC_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parallel = FALSE, ntimes = 1000)
plot(Retro_TEs_test)

numOverlaps(ATAC_AS, Retro_TEs, count.once=TRUE)


RetrosH3K27me3_te_test <- permTest(A=H3K27me3_AS, B = Retro_TEs, randomize.function = resampleRegions, universe = H3K27me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(RetrosH3K27me3_te_test)

RetrosH3K4me3_te_test <- permTest(A=H3K4me3_AS, B = Retro_TEs, randomize.function = resampleRegions, universe =H3K4me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(RetrosH3K4me3_te_test)

RetrosH3K36me3_te_test <- permTest(A=H3K36me3_AS, B = Retro_TEs, randomize.function = resampleRegions, universe = H3K36me3_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(RetrosH3K36me3_te_test)

RetrosH3K56ac_te_test <- permTest(A=H3K56ac_AS, B = Retro_TEs, randomize.function = resampleRegions, universe = H3K56ac_peaks, evaluate.function = numOverlaps, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(RetrosH3K56ac_te_test)


###

RetrosTE_AS_perm_objects <- list(Retro_TEs_test,RetrosH3K4me3_te_test,RetrosH3K36me3_te_test,RetrosH3K56ac_te_test,RetrosH3K27me3_te_test)
Marks <- c("ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3")

RetrosTE_AS_perm_objects[[1]]$Mark <- "ATAC"
RetrosTE_AS_perm_objects[[2]]$Mark <- "H3K4me3"
RetrosTE_AS_perm_objects[[3]]$Mark <- "H3K36me3"
RetrosTE_AS_perm_objects[[4]]$Mark <- "H3K56ac"
RetrosTE_AS_perm_objects[[5]]$Mark <- "H3K27me3"

RetrosTE_AS_perm_values <- list()
RetrosTE_AS_observed_numOverlap <- list()

for (i in 1:5) {
 RetrosTE_AS_perm_values[[i]] <- RetrosTE_AS_perm_objects[[i]][["numOverlaps"]][["permuted"]]
 RetrosTE_AS_perm_values[[i]] <- as.data.frame(RetrosTE_AS_perm_values[[i]])
 RetrosTE_AS_perm_values[[i]]$Type <- "Permuted"
 RetrosTE_AS_perm_values[[i]]$Mark <- Marks[i]
 colnames(RetrosTE_AS_perm_values[[i]]) <- c("Overlaps","Type","Mark")
 RetrosTE_AS_observed_numOverlap[[i]] <- RetrosTE_AS_perm_objects[[i]][["numOverlaps"]][["observed"]]
 RetrosTE_AS_observed_numOverlap[[i]] <- as.data.frame(RetrosTE_AS_observed_numOverlap[[i]])
 RetrosTE_AS_observed_numOverlap[[i]]$Type <- "Observed"
 RetrosTE_AS_observed_numOverlap[[i]]$Mark <- Marks[i]
 colnames(RetrosTE_AS_observed_numOverlap[[i]]) <- c("Overlaps","Type","Mark")
}

RetrosAll_TE_AS_overlaps <- rbind(RetrosTE_AS_perm_values[[1]],RetrosTE_AS_perm_values[[2]],RetrosTE_AS_perm_values[[3]],RetrosTE_AS_perm_values[[4]],RetrosTE_AS_perm_values[[5]],RetrosTE_AS_observed_numOverlap[[1]],RetrosTE_AS_observed_numOverlap[[2]],RetrosTE_AS_observed_numOverlap[[3]],RetrosTE_AS_observed_numOverlap[[4]],RetrosTE_AS_observed_numOverlap[[5]])




#Test overlap with All deletions small and large

#2023-04-19_add ancestry assigned dels
dels.vcf <- read.delim("delsAncestryAssigned.txt", header=FALSE)

dels.vcf <- dels.vcf %>% filter(V3 - V2 > 0)
summary(dels.vcf$V3 - dels.vcf$V2)

##Plot size

ggplot(dels.vcf, aes(x = V3 - V2)) + geom_density() + theme_classic() + xlab("Deletion Size (bp)")

dels <- filterChromosomes(toGRanges(dels.vcf[,c(1,2,3)], chr.type= "canonical"))

SVs_test <- permTest(A=ATAC_AS, B = dels, randomize.function = resampleRegions, universe = ATAC_peaks, evaluate.function = meanDistance, force.parallel = FALSE, ntimes = 1000)

plot(SVs_test)

print("permutations of deletions and as marks")
H3K27me3_SV_test <- permTest(A=H3K27me3_AS, B = dels, randomize.function = resampleRegions, universe = H3K27me3_peaks, evaluate.function = meanDistance, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K27me3_SV_test)

H3K4me3_SV_test <- permTest(A=H3K4me3_AS, B = dels, randomize.function = resampleRegions, universe =H3K4me3_peaks, evaluate.function = meanDistance, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K4me3_SV_test)

H3K36me3_SV_test <- permTest(A=H3K36me3_AS, B = dels, randomize.function = resampleRegions, universe = H3K36me3_peaks, evaluate.function = meanDistance, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K36me3_SV_test)

H3K56ac_SV_test <- permTest(A=H3K56ac_AS, B = dels, randomize.function = resampleRegions, universe = H3K56ac_peaks, evaluate.function = meanDistance, count.once = TRUE, force.parralel = FALSE, ntimes = 1000)
plot(H3K56ac_SV_test)


## ran bedtools intersect -wo -a "$SVS" -b ASE_promoter_1.sort.bed ASE_promoter_2.sort.bed ASE_promoters_both.sort.bed -names Rep1 Rep2 Both > Deletions_int_ASE_promoters.txt

Promoter_type_list <- list(ASE_promoters_1,ASE_promoters_2,ASE_bothreps)
ASE_non_geneiase_universe_results <- list()
ASE_all_geneiase_universe_res <- list()
ASE_perm_values <- list()
ASE_obs_values <- list()

print("running SV promoter permutations")

for (i in 1:3){
 ASE_non_geneiase_universe_results[[i]] <- permTest(A=Promoter_type_list[[i]], B = dels, randomize.function = resampleRegions, universe = Non_ASE_geneiase_promoters, evaluate.function = numOverlaps, count.once=TRUE, force.parallel = FALSE, ntimes = 1000)

 ASE_all_geneiase_universe_res[[i]] <- permTest(A=Promoter_type_list[[i]], B = dels, randomize.function = resampleRegions, universe = geneiase_promoters, evaluate.function = numOverlaps, count.once=TRUE, force.parallel = FALSE, ntimes = 1000)
 ASE_perm_values[[i]] <- ASE_all_geneiase_universe_res[[i]][["numOverlaps"]][["permuted"]]
 ASE_perm_values[[i]] <- as.data.frame(ASE_perm_values[[i]])
 ASE_perm_values[[i]]$Type <- "Permuted"
 colnames(ASE_perm_values[[i]]) <- c("Overlaps","Type")
 ASE_obs_values[[i]] <- ASE_all_geneiase_universe_res[[i]][["numOverlaps"]][["observed"]]
 ASE_obs_values[[i]] <- as.data.frame(ASE_obs_values[[i]])
 ASE_obs_values[[i]]$Type <- "Observed"
 colnames(ASE_obs_values[[i]]) <- c("Overlaps","Type")
}

plot(ASE_all_geneiase_universe_res[[1]])
plot(ASE_non_geneiase_universe_results[[3]])

ASE_perm_values[[1]]$Rep <- "1"
ASE_perm_values[[2]]$Rep <- "2"
ASE_perm_values[[3]]$Rep <- "ASE in both reps"

ASE_obs_values[[1]]$Rep <- "1"
ASE_obs_values[[2]]$Rep <- "2"
ASE_obs_values[[3]]$Rep <- "ASE in both reps"

ASEgenes_deletions <- rbind(ASE_perm_values[[1]],ASE_perm_values[[2]],ASE_perm_values[[3]],ASE_obs_values[[1]],ASE_obs_values[[2]],ASE_obs_values[[3]])

ASEgenes_deletions_bothrep <- rbind(ASE_perm_values[[3]],ASE_obs_values[[3]])

ggplot(ASEgenes_deletions_bothrep, aes(x = Type, fill = Type, y = Overlaps)) + geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values=c("slategrey","slateblue1")) + ylab("Number of Unique Overlaps with Deletions")+ theme(text = element_text(size=22),legend.title=element_blank()) + ylim(0,160)



Deletions_int_ASE_prom_wphaseblock_REP1 <- read.delim("Deletions_int_ASE_prom_wphaseblock_REP1.bed", header=FALSE)

Deletions_int_ASE_prom_wphaseblock_REP2 <- read.delim("Deletions_int_ASE_prom_wphaseblock_REP2.bed", header=FALSE)

Deletions_int_ASE_prom_wphaseblock_REPboth <- read.delim("Deletions_int_ASE_prom_wphaseblock_REPboth.bed", header=FALSE)

# Filter to include only genes that show ASE in both reps, rather than merging reps 1 and 2

print("Filtering deletions with nearest gene in the same phase block")
del_list <- list(Deletions_int_ASE_prom_wphaseblock_REP1, Deletions_int_ASE_prom_wphaseblock_REP2,Deletions_int_ASE_prom_wphaseblock_REPboth)

for (i in 1:3){
  colnames(del_list[[i]]) <- c("prom_chr","prom_start","prom_end","strand","geneid","del_chr","del_start","del_end","del_GT","del_gene_distance","phs_block_chr","phs_block_start","phs_block_end")
  del_list[[i]] <- del_list[[i]] %>% filter(del_start > phs_block_start & del_end< phs_block_end & prom_start > phs_block_start & prom_end < phs_block_end & del_GT != "1|1")
}


RNA1_with_pos <- read.delim("RNA1_with_pos.tab")

RNA2_with_pos <- read.delim("RNA2_with_pos.tab")

#Linked_read_SNP_GT <- read.delim("Linked_read_SNP_GT.bed", header=FALSE)

Linked_read_SNP_GT <- read.delim("SNPs_with_ancestry_assignment.txt", header=FALSE)


RNA1_counts_wSNP <- merge(RNA1_with_pos, Linked_read_SNP_GT, by.x=c("chr","pos"),by.y=c("V1","V2"))
RNA2_counts_wSNP <- merge(RNA2_with_pos, Linked_read_SNP_GT, by.x=c("chr","pos"),by.y=c("V1","V2"))


print("Calculating haplotype counts")
rna <- list(RNA1_counts_wSNP, RNA2_counts_wSNP)
sep_rep <- list()
processed <- list()
final <- list()

for (i in 1:2) {
  sep_rep[[i]] <- rna[[i]]
  sep_rep[[i]]$HAP1 <- ifelse(sep_rep[[i]]$V4 == "0|1", sep_rep[[i]]$ref.dp, sep_rep[[i]]$alt.dp)
  sep_rep[[i]]$HAP2 <- ifelse(sep_rep[[i]]$V4 == "0|1", sep_rep[[i]]$alt.dp, sep_rep[[i]]$ref.dp)
  processed[[i]] <- as.data.frame(sep_rep[[i]] %>% group_by(gene) %>% summarize(HAP1 = sum(HAP1), HAP2 =sum(HAP2)))
  processed[[i]]$HAP1_log2FC <- log2(processed[[i]]$HAP1 + 1) - log2(processed[[i]]$HAP2+1)
  processed[[i]]$HAP2_log2FC <- log2(processed[[i]]$HAP2 + 1) - log2(processed[[i]]$HAP1+1)
  final[[i]] <- merge(processed[[i]], del_list[[3]], by.x=c("gene"), by.y=c("geneid"))
  final[[i]]$HAPD_FC <- ifelse(final[[i]]$del_GT == "0|1", final[[i]]$HAP2_log2FC, final[[i]]$HAP1_log2FC)
  final[[i]]$Del_in_promoter <- ifelse(final[[i]]$del_gene_distance == 0 , "Yes","No")
}

final[[1]]$Rep <- "1"
final[[2]]$Rep <- "2"
final_both <- rbind(final[[1]],final[[2]])

ase_both <- final_both$gene %in% ASE_promoters_both$V5
final_both <- final_both[which(ase_both),]

final_both_del_1 <- final[[1]] %>% filter(Del_in_promoter == "Yes")
final_both_del_2 <- final[[2]] %>% filter(Del_in_promoter == "Yes")

a <-dim(final_both_del_1)[1]
b <-dim(final_both_del_2)[1]

final_both_avg <- final_both %>% group_by(gene) %>% summarise(HAPD_FC = mean(HAPD_FC), Del_in_promoter = Del_in_promoter)

final_both_avg <- unique(final_both_avg)

final_both_avg$HAPD_FC <- final_both_avg$HAPD_FC




deletion_ASE_genes1 <- final_both_avg %>% filter(Del_in_promoter == "Yes")
deletion_ASE_genes1 <- deletion_ASE_genes1$HAPD_FC

nodeletion_ASE_genes1 <- final_both_avg %>% filter(Del_in_promoter == "No")
nodeletion_ASE_genes1 <- nodeletion_ASE_genes1$HAPD_FC


#test stat is difference in means
test_stat <- abs(mean(nodeletion_ASE_genes1) - mean(deletion_ASE_genes1))

#combine everything
pooleddat1 <- c(deletion_ASE_genes1, nodeletion_ASE_genes1)

n.sim = 10000
nulldist1 =rep(NA, n.sim)

no <- table(final_both_avg$Del_in_promoter)[1]
total <- table(final_both_avg$Del_in_promoter)[1] + table(final_both_avg$Del_in_promoter)[2]

print("Running 10000 permutations for Haplotype D fold change for Rep1")
print(table(final[[1]]$Del_in_promoter))
set.seed(1)
for (i in 1:n.sim) {
  #randomly split data into two groups one group has 1972 samples, other has 154
  these <- sample(1:total, no)
  fakeno_del <- pooleddat1[these]
  fake_del <- pooleddat1[-these]

  # store difference in means
  nulldist1[i] <- abs(mean(fakeno_del) - mean(fake_del))

}

hist(nulldist1,main = "Difference in mean HAPD of randomly permuted groups",xlim=c(0,0.7))
abline(v = test_stat, col="red", lwd=3, lty=2)
print("Permuted test statistic (difference in mean) is greater than or equal to observed test statistic")
#calc pval
table(nulldist1 >= test_stat)

p_value <- table(nulldist1 >= test_stat)[[2]] / n.sim

##calculate effect size
library("esvis")

effect_size <- hedg_g(final_both_avg, HAPD_FC ~ Del_in_promoter)
eff_Size <- abs(effect_size$hedg_g[1])

ggplot(final_both_avg) +
  aes(x=Del_in_promoter, y=(HAPD_FC), fill = Del_in_promoter) +
  geom_boxplot(width = 0.5,outlier.size = 0.3,outlier.alpha = 0.2)+
  geom_jitter(width = 0.2,size=0.3,alpha=0.2)+
  theme(axis.text.x = element_text(angle = 20)) +
  ylab("HAPD_log2FC") + xlab("Promoter SV") + theme_bw() + labs(fill="Promoter SV")+ scale_fill_manual(values=c("slategrey","slateblue1")) + theme_classic() + theme(text = element_text(size=20)) + ggtitle(eff_Size,p_value)



print("Do SVs overlap TEs?")

Deletions_int_TEs <- read.delim("Deletions_int_TEs.bed", header=FALSE)

SVs_in_ASE_promoter_TEannotated <-read.delim("SVs_in_ASE_promoter_TEannotated.bed", header=FALSE)

SVs_in_GENEIASEpromoters_TEannotated <- read.delim("SVs_in_GENEIASEpromoters_TEannotated.bed", header=FALSE)

AS_ACRs_int_TEs <- read.delim("AS_ACRs_int_TEs.bed", header=FALSE)

TEs <- read.delim("TEs_wClass.bed", header=FALSE)


TEs$size <- TEs$V3-TEs$V2

te <- TEs %>% group_by(V6) %>% summarise(sum(size))
te <- as.data.frame(te)
te$Type <- "All_TEs"
te$Percentage <- te$`sum(size)` / sum(te$`sum(size)`) * 100
colnames(te) <- c("Var1","Freq","Type","Percentage")

all_sv <- as.data.frame(table(Deletions_int_TEs$V10))
ase_prom_sv <- as.data.frame(table(SVs_in_ASE_promoter_TEannotated$V15))
geneiase_prom_sv <- as.data.frame(table(SVs_in_GENEIASEpromoters_TEannotated$V15))

all_sv$Type <- "All_SVs"
all_sv$Percentage <- all_sv$Freq / sum(all_sv$Freq) * 100
ase_prom_sv$Type <- "ASE_promoter_SVs"
ase_prom_sv$Percentage <- ase_prom_sv$Freq / sum(ase_prom_sv$Freq) * 100
geneiase_prom_sv$Type <- "All_promoter_SVs"
geneiase_prom_sv$Percentage <- geneiase_prom_sv$Freq / sum(geneiase_prom_sv$Freq) * 100

SV_ann <- rbind(all_sv,ase_prom_sv,te,geneiase_prom_sv)
colnames(SV_ann) <- c("TE_type","Freq","Type","Percentage")

Types_detail <- unique(TEs[,c(6,7)])


Full <- merge(SV_ann, Types_detail, by.x=c("TE_type"),by.y=c("V6"))


ggplot(Full %>% filter(V7 == "DNA"),aes(fill=Type, y=Percentage, x = TE_type)) + geom_bar(position="dodge",stat="identity") + theme_classic() + scale_fill_manual(values=c("slategrey","slateblue1","slateblue4","purple")) + theme(text = element_text(size=15)) + coord_flip() + ggtitle("DNA TEs")


ggplot(Full %>% filter(V7 == "LTR"),aes(fill=Type, y=Percentage, x = TE_type)) + geom_bar(position="dodge",stat="identity") + theme_classic() + scale_fill_manual(values=c("slategrey","slateblue1","slateblue4","purple")) + theme(text = element_text(size=15)) + coord_flip() + ggtitle("LTR")




print("Separate TEs into groups")

class_2_tes <- c("MULE-MuDR","hAT","Helitron","CMC-EnSpm","PIF-Harbinger")
class_1_tes <- c("Gypsy","Copia","L1")

sub_class2_tes <- list()
regions <- list()
SV_test <- list()
SV_perm_values <-list()
SV_obs_values <- list()

ASE_prom_sv <- toGRanges(filterChromosomes(SVs_in_ASE_promoter_TEannotated[,c(1,2,3,15)]))



for (i in 1:5){
 sub_class2_tes[[i]] <- TEs %>% filter(grepl(class_2_tes[[i]],TEs$V6))
 SV_test[[i]] <- permTest(A=ASE_prom_sv, B = toGRanges(filterChromosomes(sub_class2_tes[[i]][,c(1,2,3,6)],chr.type="canonical")), randomize.function = resampleRegions, universe = dels, evaluate.function = numOverlaps, count.once=TRUE, force.parallel = FALSE, ntimes = 1000)
 SV_perm_values[[i]] <- mean(SV_test[[i]][["numOverlaps"]][["permuted"]])
 SV_perm_values[[i]] <- as.data.frame(SV_perm_values[[i]])
 SV_perm_values[[i]]$Type <- class_2_tes[i]
 SV_perm_values[[i]]$P_val <- SV_test[[i]][["numOverlaps"]][["pval"]]
 colnames(SV_perm_values[[i]]) <- c("Overlaps","Type","P_val")
 SV_obs_values[[i]] <- SV_test[[i]][["numOverlaps"]][["observed"]]
 SV_obs_values[[i]] <- as.data.frame(SV_obs_values[[i]])
 SV_obs_values[[i]]$Type <- class_2_tes[i]
 SV_obs_values[[i]]$P_val <- SV_test[[i]][["numOverlaps"]][["pval"]]
 colnames(SV_obs_values[[i]]) <- c("Overlaps","Type","P_val")
 SV_obs_values[[i]]$Log_enrichment <- log2(SV_obs_values[[i]]$Overlaps + 1) - log2(SV_perm_values[[i]]$Overlaps)
}

print("Permutations for overlaps between ASE promoter SVs and class1 TEs")

sub_class1_tes <- list()
regions <- list()
SV_test <- list()
SV_perm_values1 <-list()
SV_obs_values1 <- list()

for (i in 1:3){
 sub_class1_tes[[i]] <- TEs %>% filter(grepl(class_1_tes[[i]],TEs$V6))
 SV_test[[i]] <- permTest(A=ASE_prom_sv, B = toGRanges(filterChromosomes(sub_class1_tes[[i]][,c(1,2,3,6)],chr.type="canonical")), randomize.function = resampleRegions, universe = dels, evaluate.function = numOverlaps, count.once=TRUE, force.parallel = FALSE, ntimes = 1000)
 SV_perm_values1[[i]] <- mean(SV_test[[i]][["numOverlaps"]][["permuted"]])
 SV_perm_values1[[i]] <- as.data.frame(SV_perm_values1[[i]])
 SV_perm_values1[[i]]$Type <- class_1_tes[i]
 SV_perm_values1[[i]]$P_val <- SV_test[[i]][["numOverlaps"]][["pval"]]
 colnames(SV_perm_values1[[i]]) <- c("Overlaps","Type","P_val")
 SV_obs_values1[[i]] <- SV_test[[i]][["numOverlaps"]][["observed"]]
 SV_obs_values1[[i]] <- as.data.frame(SV_obs_values1[[i]])
 SV_obs_values1[[i]]$Type <- class_1_tes[i]
 SV_obs_values1[[i]]$P_val <- SV_test[[i]][["numOverlaps"]][["pval"]]
 colnames(SV_obs_values1[[i]]) <- c("Overlaps","Type","P_val")
 SV_obs_values1[[i]]$Log_enrichment <- log2(SV_obs_values1[[i]]$Overlaps + 1) - log2(SV_perm_values1[[i]]$Overlaps)
}

SV_obs_values1[[3]]
SV_perm_values1[[3]]

ALL_class1_obs <- do.call(rbind, SV_obs_values1)
ALLclass2_obs <-do.call(rbind, SV_obs_values)

All_obs <- rbind(ALL_class1_obs,ALLclass2_obs)

AS_ACRs_int_TEs <- read.delim("AS_ACRs_int_TEs.bed", header=FALSE)

AS_ACRs_int_TEs$V10 <- gsub("DNA/", "",AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("LTR/", "",AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("RC/", "", AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("RC/", "", AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("LINE/", "", AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("hAT.*", "hAT", AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("CMC.*", "CMC", AS_ACRs_int_TEs$V10)
AS_ACRs_int_TEs$V10 <- gsub("unknown", "Unknown_LTR", AS_ACRs_int_TEs$V10)

as_acr_te <- as.data.frame(table(AS_ACRs_int_TEs$V10))
as_acr_te$Percentage <- as_acr_te$Freq / sum(as_acr_te$Freq) * 100
colnames(as_acr_te) <- c("TE_type","Freq","Percentage")

Full$TE_type <- gsub("DNA/", "",Full$TE_type)
Full$TE_type <- gsub("LTR/", "",Full$TE_type)
Full$TE_type <- gsub("RC/", "", Full$TE_type)
Full$TE_type <- gsub("RC/", "", Full$TE_type)
Full$TE_type <- gsub("LINE/", "", Full$TE_type)
Full$TE_type <- gsub("hAT.*", "hAT", Full$TE_type)
Full$TE_type <- gsub("CMC.*", "CMC", Full$TE_type)


Full$TE_type <- gsub("unknown", "Unknown_LTR", Full$TE_type)
Full <- Full %>% filter(TE_type != "DNA" & TE_type != "LTR")
promsvs <- Full %>% filter(Type == "ASE_promoter_SVs" | Type == "All_promoter_SVs")


ggplot(promsvs,aes(fill=Type, y=Percentage, x = TE_type)) + geom_bar(position="dodge",stat="identity") + theme_classic() + scale_fill_manual(values=c("slategrey","slateblue1","slateblue4","purple")) + theme(text = element_text(size=15)) + coord_flip() + scale_x_discrete(limits=c("Unknown","Unknown_LTR","Gypsy","Copia","L1","L2","CMC","MULE-MuDR","Merlin","PIF_Harbinger","Crypton-H","Simple_repeat","Helitron"))
dev.off()

ggplot(as_acr_te,aes(y=Percentage, x = TE_type)) + geom_bar(position="dodge",stat="identity") + theme_classic() + theme(text = element_text(size=15)) + coord_flip() + scale_x_discrete(limits=c("Unknown","Unknown_LTR","Gypsy","Copia","L1","L2","CMC","hAT","MULE-MuDR","Merlin","PIF-Harbinger","TcMar-Pogo","Helitron"))
