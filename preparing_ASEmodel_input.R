library(dplyr)
library(glmnet)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(stringr)


gene_body_countsAllmods_INT <- read.delim("gene_body_countsAllmods_INT.txt", header=FALSE)

promoter_counts_INT <- read.delim("promoter_counts_INT.txt", header=FALSE)

ACRs_inwindow_UP_counts_INT <- read.delim("ACRs_inwindow_UP_counts_INT.txt", header=FALSE)

ACRs_inwindow_DOWN_counts_INT <- read.delim("ACRs_inwindow_DOWN_counts_INT.txt", header=FALSE)

#methylation, needs some processing

promoter_total_meth_int <- read.delim("promoter_total_meth_int.txt", header =F)

promoter_meth_int <- read.delim("promoter_meth_int.txt", header =F)

gene_body_total_meth <- read.delim("gene_body_meth_total_int.txt", header =F)

gene_body_meth_int <- read.delim("gene_body_meth_int.txt", header=F)

up_meth_total_int <- read.delim("ACRs_inwindow_UP_total_meth_INT.txt", header=F)

up_meth_int <- read.delim("ACRs_inwindow_UP_meth_INT.txt", header=F)

down_meth_total_int <- read.delim("ACRs_inwindow_DOWN_total_meth_INT.txt", header=F)

down_meth_int <- read.delim("ACRs_inwindow_DOWN_meth_INT.txt", header=F)

print("Processing methylation data")

total_meth_list <- list(gene_body_total_meth,promoter_total_meth_int, up_meth_total_int, down_meth_total_int)

ASE_meth_list <- list(gene_body_meth_int,promoter_meth_int, up_meth_int, down_meth_int)

regions <- c("Gene_body","Promoter","Upstream","Downstream")


total_meth_data <- list()

ASE_meth_data <- list()

meth_hap_summary <- list()

for (i in 1:2){
  colnames(ASE_meth_list[[i]]) <- c("chr","start","end","geneid","null","strand","context","chr_dup","meth_pos","meth_pos_dup","hap1_meth","hap1_unmeth","hap2_meth","hap2_unmeth")

  ASE_meth_data[[i]] <- ASE_meth_list[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarize(HAP1_meth = sum(hap1_meth), HAP2_meth = sum(hap2_meth))
  ASE_meth_data[[i]]$Region <- regions[i]

  ASE_meth_data[[i]]$Methylation_log2fc <- log2(ASE_meth_data[[i]]$HAP1_meth + 1) - log2(ASE_meth_data[[i]]$HAP2_meth + 1)

  ASE_meth_data[[i]] <- ASE_meth_data[[i]][,c(1,5)]

  colnames(total_meth_list[[i]]) <- c("chr","start","end","geneid","null","strand","context","chr_dup","meth_pos","meth_pos_dup","context2","meth")

  total_meth_data[[i]] <- total_meth_list[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarize(Meth = log2(mean(meth+0.001)))
  total_meth_data[[i]]$region <- regions[i]
  total_meth_data[[i]] <- total_meth_data[[i]][,c(1,2)]

}

for (i in 3:4){
  colnames(ASE_meth_list[[i]]) <- c("chr","start","end","geneid","null","strand","peakid","peak_start","peak_end","context","chr_dup","meth_pos","meth_pos_dup","hap1_meth","hap1_unmeth","hap2_meth","hap2_unmeth")

  ASE_meth_data[[i]] <- ASE_meth_list[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarize(HAP1_meth = sum(hap1_meth), HAP2_meth = sum(hap2_meth))
  ASE_meth_data[[i]]$Region <- regions[i]

  ASE_meth_data[[i]]$Methylation_log2fc <- log2(ASE_meth_data[[i]]$HAP1_meth + 1) - log2(ASE_meth_data[[i]]$HAP2_meth + 1)

  ASE_meth_data[[i]] <- ASE_meth_data[[i]][,c(1,5)]

  colnames(total_meth_list[[i]]) <- c("chr","start","end","geneid","null","strand","peakid","peak_start","peak_end","context","chr_dup","meth_pos","meth_pos_dup","context2","meth")

  total_meth_data[[i]] <- total_meth_list[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarize(Meth = log2(mean(meth+0.001)))
  total_meth_data[[i]]$region <- regions[i]
  total_meth_data[[i]] <- total_meth_data[[i]][,c(1,2)]

}

rm(total_meth_list,ASE_meth_list,gene_body_meth_int,promoter_meth_int,up_meth_int,down_meth_int,up_meth_total_int,down_meth_total_int,gene_body_total_meth,promoter_total_meth_int)


colnames(ASE_meth_data[[1]]) <- c("geneid","gene_body_meth")
colnames(ASE_meth_data[[2]]) <- c("geneid","promoter_meth")
colnames(ASE_meth_data[[3]]) <- c("geneid","UP_meth")
colnames(ASE_meth_data[[4]]) <- c("geneid","DOWN_meth")

colnames(total_meth_data[[1]]) <- c("geneid","gene_body_meth_tot")
colnames(total_meth_data[[2]]) <- c("geneid","promoter_meth_tot")
colnames(total_meth_data[[3]]) <- c("geneid","UP_meth_tot")
colnames(total_meth_data[[4]]) <- c("geneid","DOWN_meth_tot")

#deletion
promoter_del_INT <- read.delim("promoter_del_INT.txt", header=FALSE)

promoter_del_NO_int <- read.delim("promoter_del_NO_int.txt", header=FALSE)


##2024-04-13 INCORPORATE ancestry assigned dels
dels_ancestry <- read.delim("2023-05-11_delsAncestryAssigned.txt", header=FALSE)


promoter_delINT_anc <- merge(promoter_del_INT, dels_ancestry, by.x=c("V7","V8","V9"),by.y=c("V1","V2","V3"))

promoter_delINT_anc$HAP1_promoter_del <- ifelse(promoter_delINT_anc$V4.y == "1|0", 1, 0)
promoter_del_NO_int$HAP1_promoter_del <- 0

promoter_delINT_anc$HAP2_promoter_del <- ifelse(promoter_delINT_anc$V4.y == "1|0", 0, 1)
promoter_del_NO_int$HAP2_promoter_del <- 0

promoter_delINT_anc <- promoter_delINT_anc[,c(7,12,13)]
colnames(promoter_delINT_anc) <- c("geneid","HAP1_promoter_del","HAP2_promoter_del")
promoter_del_NO_int <- promoter_del_NO_int[,c(4,11,12)]
colnames(promoter_del_NO_int) <- c("geneid","HAP1_promoter_del","HAP2_promoter_del")


promoter_Del <- rbind(promoter_delINT_anc, promoter_del_NO_int)
colnames(promoter_Del) <- c("geneid","HAP1_promoter_del","HAP2_promoter_del")




All_genebodies <- read.delim("All_genebodies.bed", header=FALSE)
gene_strand <- All_genebodies[,c(6,4)]

##

Promoter_SNP_counts <- read.delim("Promoter_SNP_counts.txt", header=FALSE)
Promoter_SNP_counts <- Promoter_SNP_counts[,c(4,7)]
Promoter_SNP_counts$Promoter_SNPs <- log2(Promoter_SNP_counts$V7 +1)
Promoter_SNP_counts <- Promoter_SNP_counts[,c(1,3)]
colnames(Promoter_SNP_counts) <- c("geneid","Promoter_SNPs")

########################################################

gene_body_countsAllmods_INT$Region <- "Gene_body"
promoter_counts_INT$Region <- "Promoter"

#inserted snps with parent assignment

phased_SNPs_gts <- read.delim("SNPs_with_ancestry_assignment.txt", header=FALSE)

SNPEFF <- read.delim("2022-11-21_SNPeff_impacts.txt")

SNPEFF$IMPACT <- gsub(",.*", "",SNPEFF$ANN....IMPACT)
SNPEFF$Annotation <- gsub(",.*", "",SNPEFF$ANN....EFFECT)
SNPEFF <- SNPEFF[,c(1,2,3,6,7)]

head(phased_SNPs_gts)
colnames(phased_SNPs_gts) <- c("chr","pos","snp_pos_dup","GT")

phased_SNPs_gts <- merge(phased_SNPs_gts, SNPEFF, by.x=c("chr","pos"),by.y=c("CHROM","POS"))

INT_dfs <- list(gene_body_countsAllmods_INT,promoter_counts_INT)
INT_dfs_with_SNPs <- list()

for (i in 1:2) {
  print(dim(INT_dfs[[i]]))
  colnames(INT_dfs[[i]]) <- c("chr","start","end","geneid","NULL","strand","DataType","Snp_chr","snp_pos","snp_pos_dup","ref_count","alt_count","Region")
  INT_dfs_with_SNPs[[i]] <- merge(INT_dfs[[i]], phased_SNPs_gts, by.x=c("Snp_chr","snp_pos"), by.y=c("chr","pos"))
  print(dim(INT_dfs_with_SNPs[[i]]))
}

gene_body <- INT_dfs_with_SNPs[[1]]
promoter <- INT_dfs_with_SNPs[[2]]



###### process UP STREAM window AND DOWNSTREAM WINDOW separately

print("Processing up window")
UP_window <- merge(ACRs_inwindow_UP_counts_INT, phased_SNPs_gts, by.x=c("V11","V12"), by.y = c("chr","pos"))
window_UP <- UP_window[,c(3,4,5,6,8,9,10,11,12,1,2,13,14,15,16,17,18,19)]
colnames(window_UP) <- c("chr","ACR_start","ACR_end","geneid","strand","ACR","prom_start","prom_end","DataType","snp_chr","snp_pos","snp_pos_dup","ref_count","alt_count","pos_dup","GT","IMPACT","Variant_type")


DOWN_window <- merge(ACRs_inwindow_DOWN_counts_INT, phased_SNPs_gts, by.x=c("V11","V12"), by.y = c("chr","pos"))
#for some reason this rearranges things alot
window_DOWN <- DOWN_window[,c(3,4,5,6,8,9,10,11,12,1,2,13,14,15,16,17,18,19)]
colnames(window_DOWN) <- c("chr","ACR_start","ACR_end","geneid","strand","ACR","gene_start","gene_end","DataType","snp_chr","snp_pos","snp_pos_dup","ref_count","alt_count","pos_dup","GT","IMPACT","Variant_type")

#Split gene body by data type
gb_genom <- gene_body %>% filter(DataType == "Genomic")
gb_RNA1 <- gene_body %>% filter(DataType == "RNA1")
gb_RNA2 <- gene_body %>% filter(DataType == "RNA2")
gb_H3K4me3 <- gene_body %>% filter(DataType == "H3K4me3")
gb_H3K36me3 <- gene_body %>% filter(DataType == "H3K36me3")
gb_H3K56ac <- gene_body %>% filter(DataType == "H3K56ac")
gb_H3K27me3 <- gene_body %>% filter(DataType == "H3K27me3")
gb_H3K9me2 <- gene_body %>% filter(DataType == "H3K9me2")
gb_ATAC <- gene_body %>% filter(DataType == "ATAC")



prom_ATAC <- promoter %>% filter(DataType == "ATAC")
prom_H3K4me3 <- promoter %>% filter(DataType == "H3K4me3")
prom_H3K36me3 <- promoter %>% filter(DataType == "H3K36me3")
prom_H3K56ac <- promoter %>% filter(DataType == "H3K56ac")
prom_H3K27me3 <- promoter %>% filter(DataType == "H3K27me3")
prom_H3K9me2 <- promoter %>% filter(DataType == "H3K9me2")



wind_ATAC_UP <- window_UP %>% filter(DataType == "ATAC")
wind_H3K4me3_UP <- window_UP %>% filter(DataType == "H3K4me3")
wind_H3K36me3_UP <- window_UP %>% filter(DataType == "H3K36me3")
wind_H3K56ac_UP <- window_UP %>% filter(DataType == "H3K56ac")
wind_H3K27me3_UP <- window_UP %>% filter(DataType == "H3K27me3")
wind_H3K9me2_UP <- window_UP %>% filter(DataType == "H3K9me2")


wind_ATAC_DOWN <- window_DOWN %>% filter(DataType == "ATAC")
wind_H3K4me3_DOWN <- window_DOWN %>% filter(DataType == "H3K4me3")
wind_H3K36me3_DOWN <- window_DOWN %>% filter(DataType == "H3K36me3")
wind_H3K56ac_DOWN <- window_DOWN %>% filter(DataType == "H3K56ac")
wind_H3K27me3_DOWN <- window_DOWN %>% filter(DataType == "H3K27me3")
wind_H3K9me2_DOWN <- window_DOWN %>% filter(DataType == "H3K9me2")

cov_dfs_names <- list.files("/mosdepth/",full.names = TRUE)

cov_dfs <- lapply(cov_dfs_names, function(i){
  read.delim(i,head=FALSE)})

all_data <- list(gb_genom, gb_RNA1, gb_RNA2, gb_H3K4me3, gb_H3K36me3,gb_H3K56ac, gb_H3K27me3,gb_H3K9me2,gb_ATAC,prom_ATAC,prom_H3K4me3,prom_H3K36me3,prom_H3K56ac,prom_H3K27me3,prom_H3K9me2,wind_ATAC_UP,wind_H3K4me3_UP,wind_H3K36me3_UP,wind_H3K56ac_UP,wind_H3K27me3_UP,wind_H3K9me2_UP,wind_ATAC_DOWN,wind_H3K4me3_DOWN,wind_H3K36me3_DOWN,wind_H3K56ac_DOWN,wind_H3K27me3_DOWN,wind_H3K9me2_DOWN)

count_sum <- list()

data_type_vector <- c("Genomic","RNA1","RNA2","H3K4me3","H3K36me3","H3K56ac","H3K27me3","H3K9me2","ATAC","ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3","H3K9me2","ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3","H3K9me2","ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3","H3K9me2")

region_type_vector <- c("Gene_body","Gene_body","Gene_body","Gene_body","Gene_body","Gene_body","Gene_body","Gene_body","Gene_body","Promoter","Promoter","Promoter","Promoter","Promoter","Promoter","Upstream","Upstream","Upstream","Upstream","Upstream","Upstream","Downstream","Downstream","Downstream","Downstream","Downstream","Downstream")



for (i in 1:27) {
  all_data[[i]]$HAP1 <- ifelse(all_data[[i]]$GT == "0|1", all_data[[i]]$ref_count, all_data[[i]]$alt_count)
  all_data[[i]]$HAP2 <- ifelse(all_data[[i]]$GT == "1|0", all_data[[i]]$ref_count, all_data[[i]]$alt_count)

  count_sum[[i]] <- as.data.frame(all_data[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarise(HAP1 = sum(HAP1), HAP2 = sum(HAP2)))
  count_sum[[i]]$HAP1_log2FC <- log2(count_sum[[i]]$HAP1+1) - log2(count_sum[[i]]$HAP2+1)
  count_sum[[i]]$Coverage <- log2(count_sum[[i]]$HAP1 + count_sum[[i]]$HAP2)
  print(data_type_vector[i])
  print(region_type_vector[i])
  count_sum[[i]] <- count_sum[[i]][,c(1,4,5)]
  print("Dim before filtering")
  print(dim(count_sum[[i]]))
  print(quantile(count_sum[[i]]$Coverage,probs=seq(0,1,0.1)))
  count_sum[[i]] <- count_sum[[i]] %>% filter(Coverage > quantile(count_sum[[i]]$Coverage,probs=seq(0,1,0.1))[[2]])
  count_sum[[i]] <- count_sum[[i]][,c(1,2)]
  print("Dim after filtering based on feature coverage")
  print(dim(count_sum[[i]]))
}


for (i in 1:34) {
  print(head(cov_dfs_names[[i]]),1)
  print(head(cov_dfs[[i]],3))
}

data_type_vector2 <- c("ATAC","ATAC","ATAC","ATAC","Meth","Meth","Meth","Meth","Genomic","Genomic","Genomic","Genomic","H3K27me3","H3K27me3","H3K27me3","H3K27me3","H3K36me3","H3K36me3","H3K36me3","H3K36me3","H3K4me3","H3K4me3","H3K4me3",
                       "H3K4me3","H3K56ac","H3K56ac","H3K56ac","H3K56ac","H3K9me2","H3K9me2","H3K9me2","H3K9me2","RNA1","RNA2")
region_type_vector2 <- c("DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","DOWN","UP","genebody","promoter","genebody","genebody")



Coverage_data <- list()
for (i in 1:34){

  cov_dfs[[i]]$geneid <- cov_dfs[[i]]$V4
  cov_dfs[[i]]$Coverage_raw <- cov_dfs[[i]]$V5

  hist(cov_dfs[[i]]$Coverage_raw , main = gsub("/mosdepth//","",cov_dfs_names[i]))
  cov_dfs[[i]]$Coverage <- log2(cov_dfs[[i]]$Coverage_raw + 1)

  Coverage_data[[i]] <- as.data.frame(cov_dfs[[i]] %>% dplyr::group_by(geneid) %>% dplyr::summarize(Coverage = mean(Coverage)))

  Coverage_data[[i]]$DataType <- data_type_vector2[i]
  Coverage_data[[i]]$region <- region_type_vector2[i]
  Coverage_data[[i]] <- Coverage_data[[i]] %>% dplyr::select(geneid,DataType,region,Coverage)
  Coverage_data[[i]] <- Coverage_data[[i]] %>% filter(Coverage != -Inf)

  print(ggplot(Coverage_data[[i]], aes(x = Coverage))+geom_histogram()+ggtitle(data_type_vector2[i],region_type_vector2[i]) +theme_classic() + xlab("log(Coverage)"))
  Coverage_data[[i]] <- Coverage_data[[i]][,c(1,4)]

  }


print(ggplot(Coverage_data[[33]], aes(x = Coverage))+geom_histogram()+ggtitle(data_type_vector2[33],region_type_vector2[33]) +theme_classic() + xlab("log(Coverage)"))
print(ggplot(Coverage_data[[34]], aes(x = Coverage))+geom_histogram()+ggtitle(data_type_vector2[34],region_type_vector2[34]) +theme_classic() + xlab("log(Coverage)"))



RNA_TPMs_bothreps <- read.delim("RNA_TPMs_bothreps.txt")

RNA_TPM_R1 <- RNA_TPMs_bothreps[,c(1,2)]
RNA_TPM_R2 <- RNA_TPMs_bothreps[,c(1,3)]

Coverage_data[[33]] <- RNA_TPM_R1
Coverage_data[[34]] <- RNA_TPM_R2

colnames(Coverage_data[[33]]) <- c("geneid","Coverage")
colnames(Coverage_data[[34]]) <- c("geneid","Coverage")

Coverage_data[[33]]$Coverage <- log2(Coverage_data[[33]]$Coverage +1)
Coverage_data[[34]]$Coverage <- log2(Coverage_data[[34]]$Coverage+1)

Coverage_data[[33]] <- Coverage_data[[33]] %>% filter(Coverage != -Inf)
Coverage_data[[34]] <- Coverage_data[[34]] %>% filter(Coverage != -Inf)


Coverage_data[[5]] <- total_meth_data[[4]] #downstream
Coverage_data[[6]] <- total_meth_data[[3]] #upstream
Coverage_data[[7]] <- total_meth_data[[2]] #promoter
Coverage_data[[8]] <- total_meth_data[[1]] #genebody




ggplot(all_data[[1]], aes(x = factor(Annotation), fill = IMPACT)) +
  geom_bar(stat = "count", position = position_dodge()) +
  coord_flip() + facet_wrap(~IMPACT, scales = "free") + theme_classic()

synon <- all_data[[1]] %>% filter(Annotation == "synonymous_variant")

nonsynon <- all_data[[1]] %>% filter(Annotation != "splice_donor_variant&intron_variant" &
                                       Annotation != "splice_acceptor_variant&intron_variant" & Annotation != "stop_retained_variant" & Annotation != "splice_region_variant" & IMPACT != "LOW" & IMPACT != "MODIFIER")

corrected <- rbind(synon, nonsynon)

ggplot(corrected, aes(x = factor(Annotation), fill = IMPACT)) +
  geom_bar(stat = "count", position = position_dodge()) +
  coord_flip() + facet_wrap(~IMPACT, scales = "free") + theme_classic()

#Include only SNPs where the alternate allele from DPI matches alternate allele in Fairchild
dpi_SNPs_final_GTs <- read.table("dpi_SNPs_final_GTs.txt", quote="\"", comment.char="")

fairchild_snps_for_polarization <- read.delim("LinkRead_Pos_ForPolarization.txt", header=FALSE)

combined <- merge(fairchild_snps_for_polarization, dpi_SNPs_final_GTs, by=c("V1","V2"))
combine_match <- combined %>% filter(V3.x == V3.y)
rm(combined)
combine_match <- combine_match[,c(1,2,5)]
colnames(combine_match) <- c("Snp_chr","snp_pos","pol_snp_trif_GT")
head(combine_match)


corrected_for_polarization <- merge(corrected, combine_match, by=c("Snp_chr","snp_pos"))

ggplot(corrected_for_polarization, aes(x = factor(Annotation), fill = IMPACT)) +
  geom_bar(stat = "count", position = position_dodge()) +
  coord_flip() + facet_wrap(~IMPACT, scales = "free") + theme_classic()
#make uniform GT field
corrected_for_polarization$pol_snp_trif_GT <- gsub("\\|","/",corrected_for_polarization$pol_snp_trif_GT)

## Summarize SNP effect data for SNPs that are present in WGS read counts overlapping gene bodies
corrected_summary <- corrected_for_polarization %>% dplyr::group_by(geneid,GT,pol_snp_trif_GT,IMPACT) %>% dplyr::summarize(N=n()) %>%
  mutate(freq = N / sum(N)) %>%
  ungroup() %>%
  complete(geneid,GT,pol_snp_trif_GT,IMPACT, fill = list(N=0,freq=0))


####
#a means GT = 1|0
#b means GT = 0|1
synon_hap1_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "LOW" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, synon_hap1_a = 5)

synon_hap1_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "LOW" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, synon_hap1_b = 5)

synon_hap2_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "LOW" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, synon_hap2_a = 5)

synon_hap2_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "LOW" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, synon_hap2_b = 5)

non_synon_mod_hap1_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "MODERATE" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, non_synon_mod_hap1_a = 5)

non_synon_mod_hap1_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "MODERATE" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, non_synon_mod_hap1_b = 5)

non_synon_mod_hap2_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "MODERATE" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, non_synon_mod_hap2_a = 5)

non_synon_mod_hap2_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "MODERATE" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, non_synon_mod_hap2_b = 5)

non_synon_high_hap1_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "HIGH" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, non_synon_high_hap1_a = 5)

non_synon_high_hap1_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "HIGH" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, non_synon_high_hap1_b = 5)

non_synon_high_hap2_a <- corrected_summary %>% filter(GT == "1|0" & IMPACT == "HIGH" & pol_snp_trif_GT == "1/1") %>% dplyr::select(geneid = 1, non_synon_high_hap2_a = 5)

non_synon_high_hap2_b <- corrected_summary %>% filter(GT == "0|1" & IMPACT == "HIGH" & pol_snp_trif_GT == "0/0") %>% dplyr::select(geneid = 1, non_synon_high_hap2_b = 5)


corrected_SNP_EFFECT_LIST <- list(synon_hap1_a,synon_hap1_b,synon_hap2_a,synon_hap2_b,non_synon_mod_hap1_a,non_synon_mod_hap1_b, non_synon_mod_hap2_a,non_synon_mod_hap2_b,non_synon_high_hap1_a,non_synon_high_hap1_b,non_synon_high_hap2_a,non_synon_high_hap2_b)

corrected_SNP_EFFECT_comp <- Reduce(function(x, y) merge(x, y, all=TRUE), corrected_SNP_EFFECT_LIST)

corrected_SNP_EFFECT_comp$non_synon_hap1 <- corrected_SNP_EFFECT_comp$non_synon_mod_hap1_a + corrected_SNP_EFFECT_comp$non_synon_mod_hap1_b +
  corrected_SNP_EFFECT_comp$non_synon_high_hap1_a + corrected_SNP_EFFECT_comp$non_synon_high_hap1_b

corrected_SNP_EFFECT_comp$non_synon_hap2 <- corrected_SNP_EFFECT_comp$non_synon_mod_hap2_a + corrected_SNP_EFFECT_comp$non_synon_mod_hap2_b +
  corrected_SNP_EFFECT_comp$non_synon_high_hap2_a + corrected_SNP_EFFECT_comp$non_synon_high_hap2_b

corrected_SNP_EFFECT_comp$synon_hap1 <- corrected_SNP_EFFECT_comp$synon_hap1_a + corrected_SNP_EFFECT_comp$synon_hap1_b

corrected_SNP_EFFECT_comp$synon_hap2 <- corrected_SNP_EFFECT_comp$synon_hap2_a + corrected_SNP_EFFECT_comp$synon_hap2_b

plot(corrected_SNP_EFFECT_comp$non_synon_hap1, corrected_SNP_EFFECT_comp$non_synon_hap2)

plot(corrected_SNP_EFFECT_comp$synon_hap1, corrected_SNP_EFFECT_comp$synon_hap2)

full_snpeff <- corrected_SNP_EFFECT_comp[,c(1,14,15,16,17)]

full_snpeff$dn_ds_hap1 <- log2(full_snpeff$non_synon_hap1 +1) - log2(full_snpeff$synon_hap1+1)
full_snpeff$dn_ds_hap2 <- log2(full_snpeff$non_synon_hap2 +1) - log2(full_snpeff$synon_hap2+1)

full_snpeff <- full_snpeff[,c(1,6,7)]

summary(full_snpeff$dn_ds_hap1)


hist(full_snpeff$dn_ds_hap1,breaks=8)
hist(full_snpeff$dn_ds_hap2,breaks=8)
dev.off()

for(i in 1:27){
  print(i)
  print(data_type_vector[i])
  print(region_type_vector[i])
}


#Name columns appropriately
colnames(count_sum[[1]]) <- c("geneid","Genomic")
colnames(count_sum[[2]]) <- c("geneid","RNA1")
colnames(count_sum[[3]]) <- c("geneid","RNA2")
colnames(count_sum[[4]]) <- c("geneid","H3K4me3")
colnames(count_sum[[5]]) <- c("geneid","H3K36me3")
colnames(count_sum[[6]]) <- c("geneid","H3K56ac")
colnames(count_sum[[7]]) <- c("geneid","H3K27me3")
colnames(count_sum[[8]]) <- c("geneid","H3K9me2")
colnames(count_sum[[9]]) <- c("geneid","ATAC")
colnames(count_sum[[10]]) <- c("geneid","Promoter_ATAC")
colnames(count_sum[[11]]) <- c("geneid","Promoter_H3K4me3")
colnames(count_sum[[12]]) <- c("geneid","Promoter_H3K36me3")
colnames(count_sum[[13]]) <- c("geneid","Promoter_H3K56ac")
colnames(count_sum[[14]]) <- c("geneid","Promoter_H3K27me3")
colnames(count_sum[[15]]) <- c("geneid","Promoter_H3K9me2")
colnames(count_sum[[16]]) <- c("geneid","UP_ATAC")
colnames(count_sum[[17]]) <- c("geneid","UP_H3K4me3")
colnames(count_sum[[18]]) <- c("geneid","UP_H3K36me3")
colnames(count_sum[[19]]) <- c("geneid","UP_H3K56ac")
colnames(count_sum[[20]]) <- c("geneid","UP_H3K27me3")
colnames(count_sum[[21]]) <- c("geneid","UP_H3K9me2")
colnames(count_sum[[22]]) <- c("geneid","DOWN_ATAC")
colnames(count_sum[[23]]) <- c("geneid","DOWN_H3K4me3")
colnames(count_sum[[24]]) <- c("geneid","DOWN_H3K36me3")
colnames(count_sum[[25]]) <- c("geneid","DOWN_H3K56ac")
colnames(count_sum[[26]]) <- c("geneid","DOWN_H3K27me3")
colnames(count_sum[[27]]) <- c("geneid","DOWN_H3K9me2")

#Name columns for Coverage
for(i in 1:34){
  print(i)
  print(data_type_vector2[i])
  print(region_type_vector2[i])
}


colnames(Coverage_data[[1]]) <- c("geneid","DOWN_ATAC_cov")
colnames(Coverage_data[[2]]) <- c("geneid","UP_ATAC_cov")
colnames(Coverage_data[[3]]) <- c("geneid","ATAC_cov")
colnames(Coverage_data[[4]]) <- c("geneid","Promoter_ATAC_cov")
colnames(Coverage_data[[5]]) <- c("geneid","DOWN_Meth_cov")
colnames(Coverage_data[[6]]) <- c("geneid","UP_Meth_cov")
colnames(Coverage_data[[7]]) <- c("geneid","Promoter_Meth_cov")
colnames(Coverage_data[[8]]) <- c("geneid","Meth_cov")

colnames(Coverage_data[[9]]) <- c("geneid","DOWN_Genomic_cov")
colnames(Coverage_data[[10]]) <- c("geneid","UP_Genomic_cov")
colnames(Coverage_data[[11]]) <- c("geneid","Genomic_cov")
colnames(Coverage_data[[12]]) <- c("geneid","Promoter_Genomic_cov")
colnames(Coverage_data[[13]]) <- c("geneid","DOWN_H3K27me3_cov")
colnames(Coverage_data[[14]]) <- c("geneid","UP_H3K27me3_cov")
colnames(Coverage_data[[15]]) <- c("geneid","H3K27me3_cov")
colnames(Coverage_data[[16]]) <- c("geneid","Promoter_H3K27me3_cov")
colnames(Coverage_data[[17]]) <- c("geneid","DOWN_H3K36me3_cov")
colnames(Coverage_data[[18]]) <- c("geneid","UP_H3K36me3_cov")
colnames(Coverage_data[[19]]) <- c("geneid","H3K36me3_cov")
colnames(Coverage_data[[20]]) <- c("geneid","Promoter_H3K36me3_cov")
colnames(Coverage_data[[21]]) <- c("geneid","DOWN_H3K4me3_cov")
colnames(Coverage_data[[22]]) <- c("geneid","UP_H3K4me3_cov")
colnames(Coverage_data[[23]]) <- c("geneid","H3K4me3_cov")
colnames(Coverage_data[[24]]) <- c("geneid","Promoter_H3K4me3_cov")
colnames(Coverage_data[[25]]) <- c("geneid","DOWN_H3K56ac_cov")
colnames(Coverage_data[[26]]) <- c("geneid","UP_H3K56ac_cov")
colnames(Coverage_data[[27]]) <- c("geneid","H3K56ac_cov")
colnames(Coverage_data[[28]]) <- c("geneid","Promoter_H3K56ac_cov")
colnames(Coverage_data[[29]]) <- c("geneid","DOWN_H3K9me2_cov")
colnames(Coverage_data[[30]]) <- c("geneid","UP_H3K9me2_cov")
colnames(Coverage_data[[31]]) <- c("geneid","H3K9me2_cov")
colnames(Coverage_data[[32]]) <- c("geneid","Promoter_H3K9me2_cov")
colnames(Coverage_data[[33]]) <- c("geneid","RNA1_cov")
colnames(Coverage_data[[34]]) <- c("geneid","RNA2_cov")




count_sum_full <- count_sum %>% reduce(full_join,by=c("geneid"))

count_sum_list <- list()

rm(count_sum_list)



cov_sum_full <- Coverage_data %>% reduce(full_join, by=c("geneid"))



print("Merge all modeling data together")
All_step1 <- merge(count_sum_full,cov_sum_full, by=c("geneid"),all=TRUE)
All_step2 <- merge(All_step1,ASE_meth_data[[1]],by=c("geneid"),all=TRUE)
All_step3 <- merge(All_step2,ASE_meth_data[[2]], by =c("geneid"),all= TRUE)
All_step4 <- merge(All_step3,ASE_meth_data[[3]], by =c("geneid"),all= TRUE)
All_step5 <- merge(All_step4,ASE_meth_data[[4]], by =c("geneid"),all= TRUE)

All_step6 <- merge(All_step5,promoter_Del, by =c("geneid"),all= TRUE)
All_step7 <- merge(All_step6, full_snpeff, by=c("geneid"),all=TRUE )


All_step7$RNA_BOTH <- (All_step7$RNA1 + All_step7$RNA2) / 2
All_step7$RNA_cov_both <- (All_step7$RNA1_cov + All_step7$RNA2_cov) / 2

All_complete_response <- All_step7 %>% filter(RNA1 != "NA" & RNA2 != "NA")


##Filter to keep only ASE genes
ASE_promoter1_info <- read.delim("ASE_promoter_1.bed", header =FALSE)

ASE_promoter2_info <- read.delim("ASE_promoter_2.bed", header =FALSE)

ASE_prom_both <- merge(ASE_promoter1_info, ASE_promoter2_info, by=c("V5"))[,c(1)]

ASE_promoters_in_model <- merge(ASE_promoter1_info, ASE_promoter2_info, by=c("V5"))[,c(2:5,1)]


print("Filtering genes with no ASE model response")

no_response <- All_step7 %>% filter(is.na(RNA1) | is.na(RNA2))

no_response_ase <- no_response$geneid %in% ASE_prom_both
noRESPONSE_ASE <- unique(no_response[which(no_response_ase),])


onlyASE <- All_complete_response$geneid %in% ASE_prom_both
All_complete_ASE <- unique(All_complete_response[which(onlyASE),])

only_ASE_total <- All_step7$geneid %in% ASE_prom_both
All_ASE <- unique(All_step7[which(only_ASE_total),])

#Filter genes with genomic bias
lowerq = quantile(All_complete_ASE$Genomic, na.rm = TRUE)[2]
upperq = quantile(All_complete_ASE$Genomic, na.rm = TRUE)[4]
iqr = upperq - lowerq #Or use IQR(data

extreme.threshold.upper = (iqr * 1.5) + upperq
extreme.threshold.lower = lowerq - (iqr * 1.5)


#filter ASE genes with outlier WGS coverage

hist(All_complete_ASE$Genomic_cov)

lowerq_cov <- quantile(All_complete_ASE$Genomic_cov, na.rm=TRUE)[2]
upperq_cov <- quantile(All_complete_ASE$Genomic_cov, na.rm =T)[4]

iqr = upperq_cov - lowerq_cov

extreme.threshold.lower_cov = lowerq_cov- (iqr * 3)



All_complete_ASE_filt <- All_complete_ASE %>% filter(Genomic <= extreme.threshold.upper & Genomic >= extreme.threshold.lower & Genomic_cov >= extreme.threshold.lower_cov)

print("Number of WGS read mapping bias")
print(dim(All_complete_ASE))[[1]] - (dim(All_complete_ASE_filt)[[1]])



All_complete_ASE_forOverall_exp <- All_complete_ASE %>% filter(Genomic_cov >= extreme.threshold.lower_cov)


#Filter genes with strong genomic bias from ALL GENES
#Filter genes with strong genomic bias
lowerq_all = quantile(All_step7$Genomic, na.rm = TRUE)[2]
upperq_all = quantile(All_step7$Genomic, na.rm = TRUE)[4]
iqr_all = upperq_all - lowerq_all #Or use IQR(data

extreme.threshold.upper_all = (iqr_all * 1.5) + upperq_all
extreme.threshold.lower_all = lowerq_all - (iqr_all * 1.5)

#filter genes with low/high coverage from all genes
hist(All_step7$Genomic_cov)
hist(All_complete_response$Genomic_cov)


hist(All_step7$Genomic_cov)

lowerq_cov <- quantile(All_step7$Genomic_cov, na.rm=TRUE)[2]
upperq_cov <- quantile(All_step7$Genomic_cov, na.rm =T)[4]

iqr = upperq_cov - lowerq_cov

extreme.threshold.lower_cov = lowerq_cov- (iqr * 1.5)
extreme.threshold.upper_cov = (iqr * 1.5) + upperq_cov


All_complete_response_filt <- All_complete_response %>% filter(Genomic <= extreme.threshold.upper_all & Genomic >= extreme.threshold.lower_all & Genomic_cov >= extreme.threshold.lower_cov )


All_forexp <- All_step7 %>% filter(Genomic_cov >=extreme.threshold.lower_cov & Genomic_cov <= extreme.threshold.upper_cov)


b <- dim(All_ASE)[1]


variables <- data.frame(colnames(All_ASE))
missing_data <- vector()

for (i in 1:72) {
  missing_data[i] <- table(is.na(All_ASE[,c(i)]))[2] / b * 100
}

missing_data_table <- cbind(variables,missing_data)
median(missing_data_table$missing_data,na.rm = T)
write.table(missing_data_table,"/ASE_genes_Missing_data_report_Featlvl.txt", sep = "\t", col.names =TRUE, row.names=FALSE,quote=FALSE)


#Plot summarys of Data
cov <- All_complete_ASE[,c(29:62,72)]
ratio <- All_complete_ASE[,c(2:28,63,64,65,66,71)]


for (i in 1:length(cov)) {
  print(ggplot(cov, aes(x=2^cov[,c(i)])) + geom_histogram()+ geom_vline(xintercept = mean(2^cov[,c(i)]),col = 'red') + ggtitle(colnames(cov)[i]) + theme_classic()+ xlab("Sum of read counts overlapping SNPs"))

}


for (i in 1:length(ratio)) {
  print(ggplot(ratio, aes(x=ratio[,c(i)])) + geom_histogram(bins=100) + ggtitle(colnames(ratio)[i]) + theme_classic()+ xlab("HAP1/HAP2 ratio"))

}



#Missing data for All genes

b <- dim(All_step7)[1]

variables <- data.frame(colnames(All_step7))
missing_data <- vector()

for (i in 1:72) {
  missing_data[i] <- table(is.na(All_step7[,c(i)]))[2] / b * 100
}

missing_data_table <- cbind(variables,missing_data)
median(missing_data_table$missing_data,na.rm = T)
write.table(missing_data_table,"ALL_genes_Missing_data_report_Featlvl.txt", sep = "\t", col.names =TRUE, row.names=FALSE,quote=FALSE)



write.table(All_complete_ASE_forOverall_exp,"ASE_genes_for_overall_expression.txt", col.names = TRUE, row.names=FALSE, sep = "\t", quote=FALSE)

write.table(All_complete_ASE_filt,"ASE_genes_for_ASE.txt", col.names = TRUE, row.names=FALSE, sep = "\t", quote=FALSE)

write.table(All_forexp, "ALL_genes_foroverall_expression.txt", col.names=TRUE, row.names = FALSE, sep = "\t",quote=FALSE)

write.table(All_complete_response_filt, "ALL_genes_for_ASE.txt", col.names=TRUE, row.names = FALSE, sep = "\t",quote=FALSE)
