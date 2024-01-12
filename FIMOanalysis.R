#2023-07-20
#parse and analyse ERF48 fimo data
library(tidyverse)
library(dplyr)
library(qvalue)
library(gghighlight)


fimo_df <- read.delim("ERF48_pat_fimo.txt",head=T)


AS_ACR_chromatin_counts <- read.delim("/bigdata/seymourlab/idiaz026/Results/Fairchild/Annotate_feats/ACR_summits/For_meme/AS_ACR_chromatin_counts.txt", header=FALSE)

phased_SNPs_gts <- read.delim("/bigdata/seymourlab/idiaz026/Results/Fairchild/Assign_parentage_haplotype_blocks/2023-03-21_SNPs_assigned_final/2023-05-04_SNPs_with_ancestry_assignment.txt", header=FALSE)

ACR_countswSNPs <- merge(AS_ACR_chromatin_counts, phased_SNPs_gts, by.x=c("V6","V7"), by.y=c("V1","V2"))

colnames(ACR_countswSNPs) <- c("chr","snp_pos","chrb","peak_start","peak_end","peak_id","DataType","snp_posdup","ref_count","alt_count","snp_posdup2","GT")

Acrs_bed <- ACR_countswSNPs[,c(1,4,5,6)]

fimo_scores <- list()

fimo_df$sequence.name <- gsub("_.*","",fimo_df$sequence.name)
fimo_df $chr <- str_split_i(fimo_df$sequence.name, pattern = ":",1)
fimo_df$peak_start <- str_split_i(str_split_i(fimo_df$sequence.name, pattern = "-",1),pattern=":",2)
fimo_df$peak_end <- str_split_i(fimo_df$sequence.name, pattern = "-",2)

fimo_df$V13 <- fimo_df$q.value
colnames(fimo_df) <- c(colnames(fimo_df[1:12]),fimo_motifs[i])
fimo_df <- merge(fimo_df,Acrs_bed,by=c("chr","peak_start") )

print("subset")
fimo_df <- unique(fimo_df)
print(dim(fimo_df))

fimo_scores <- fimo_df %>% group_by(peak_id) %>% summarise(qvalue= q.value,score = mean(score))
fimo_scores <- as.data.frame(fimo_scores)
fimo_scores$q_value_overall <- fimo_scores$score

fimo_scores <- fimo_scores[,c(1,4)]
fimo_scores <- unique(fimo_scores)
colnames(fimo_scores) <- c("peak_id","ERF48_pat")
print("dim after")
print(dim(fimo_scores))

print("merge motifs data together")
final_motif_df <- fimo_scores %>% reduce(full_join,by='peak_id')


#split_by_data_type

ACR_ATAC <- ACR_countswSNPs %>% filter(DataType == "ATAC")
ACR_H3K4me3 <- ACR_countswSNPs %>% filter(DataType == "H3K4me3")
ACR_H3K36me3 <- ACR_countswSNPs %>% filter(DataType == "H3K36me3")
ACR_H3K56ac <- ACR_countswSNPs %>% filter(DataType == "H3K56ac")
ACR_H3K27me3 <- ACR_countswSNPs %>% filter(DataType == "H3K27me3")


all_data <- list(ACR_ATAC,ACR_H3K4me3,ACR_H3K36me3,ACR_H3K56ac,ACR_H3K27me3)

datatype_vector <- c("ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3")

count_sum <- list()

for (i in 1:5) {
  all_data[[i]]$HAP1 <- ifelse(all_data[[i]]$GT == "0|1", all_data[[i]]$ref_count, all_data[[i]]$alt_count)
  all_data[[i]]$HAP2 <- ifelse(all_data[[i]]$GT == "1|0", all_data[[i]]$ref_count, all_data[[i]]$alt_count)

  count_sum[[i]] <- as.data.frame(all_data[[i]] %>% dplyr::group_by(peak_id) %>% dplyr::summarise(HAP1 = sum(HAP1), HAP2 = sum(HAP2)))
  count_sum[[i]]$HAP1_log2FC <- log2(count_sum[[i]]$HAP1+1) - log2(count_sum[[i]]$HAP2+1)
  count_sum[[i]]$Coverage <- log2(count_sum[[i]]$HAP1 + count_sum[[i]]$HAP2)
  count_sum[[i]]$Type <- datatype_vector[i]
  count_sum[[i]] <- count_sum[[i]][,c(1,4)]
  colnames(count_sum[[i]]) <- c("peak_id",datatype_vector[i])

}


final_chrom_df <- count_sum %>% reduce(full_join,by='peak_id')

datatype_vector <- c("ATAC","H3K4me3","H3K36me3","H3K56ac","H3K27me3")


for (i in 1:5) {
  p <- ggscatter(full_meta, x = ERF48_pat, datatype_vector[i],add = "reg.line", conf.int = F, size = 1) + stat_cor(aes(),method = "spearman",label.sep = sprintf(", n = %s, ", sum(complete.cases(full_meta[c(fimo_motifs[i],datatype_vector[1])]))))
  print(p) + geom_vline(aes(xintercept=1.30103))

}

dev.off()

head(full_meta)



for (i in 1:5) {
  p <- ggscatter(full_meta, x = datatype_vector[1], datatype_vector[i],add = "reg.line", conf.int = F, size = 1) + stat_cor(aes(),method = "spearman",label.sep = sprintf(", n = %s, ", sum(complete.cases(full_meta[c(datatype_vector[1],datatype_vector[i])]))))

  print(p)
  }


###############
# Test fimo scores of motif3 (repressive TF) in ase promoter svs

ase_prom_sv_fimo <- read.delim("AZF1_ase_promSV_fimo.txt", header=TRUE)

ase_prom_sv_fimo$sequence.name <- gsub("_.*","",ase_prom_sv_fimo$sequence.name)
ase_prom_sv_fimo$chr <- str_split_i(ase_prom_sv_fimo$sequence.name, pattern = ":",1)
ase_prom_sv_fimo$sv_start <- str_split_i(str_split_i(ase_prom_sv_fimo$sequence.name, pattern = "-",1),pattern=":",2)
ase_prom_sv_fimo$sv_end <- str_split_i(ase_prom_sv_fimo$sequence.name, pattern = "-",2)

ase_prom_sv_fimo$V13 <- ase_prom_sv_fimo$q.value


sv_int_aseprom <- read.delim("Deletions_int_ASE_RepBothpromoters.bed", header=FALSE)[,c(6,7,8,5,9)]

ASE_DATA <- read.delim("SV_ASE_data.txt")

aspromsvFIMO <- ase_prom_sv_fimo %>% group_by(sequence.name) %>% summarise(chr = chr,sv_start = sv_start, sv_end = sv_end, score = mean(score), q_val = mean(q.value))

aspromsvFIMO <- unique(aspromsvFIMO)
#merge fimo w SVs

fimo_sv <- merge(aspromsvFIMO, sv_int_aseprom, by.x=c("chr","sv_start","sv_end"),by.y=c("V6","V7","V8"))


head(fimo_sv)


plot(fimo_sv$score, fimo_sv$q_val)

#merge_w_asedata

colnames(fimo_sv) <- c("chr","sv_start","sv_end","sequence.name","score","q_val","geneid","GT")

meta_sv <- unique(merge(fimo_sv, ASE_DATA ,by.x=c("geneid"),by.y=c("gene")))


sig_motif3 <- meta_sv$geneid


in_sigmotif <- ASE_DATA$gene %in% sig_motif3


ASE_DATA$AZF1 <- ifelse(ASE_DATA$gene %in% sig_motif3, "Sig", "Not_Sig")


ASE_DATA_filt <- ASE_DATA %>% filter(Del_in_promoter == "Yes")

data <- ASE_DATA_filt[order(ASE_DATA_filt$HAPD_FC), ]



ggplot(data, aes(x=reorder(gene,-(HAPD_FC)), y=HAPD_FC,fill = AZF1)) +
  geom_bar(stat='identity',aes(fill=AZF1),width = 0.6)+ theme_classic()+ theme(axis.text.y = element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +xlab("Gene") +
  ylab('log2(HAPD_FC)')  + gghighlight(AZF1 == "Sig")  + ggtitle("ASE genes with promoter SVs containing AZF1 motif") + scale_fill_manual(values=c("slateblue1"))



AZF1 <- data %>% filter(AZF1 == "Sig")
AZF1 <- AZF1$HAPD_FC

noAZF1 <- data %>% filter(AZF1 == "Not_Sig")
noAZF1 <- noAZF1$HAPD_FC


#test stat is difference in means
test_stat <- abs(mean(AZF1) - mean(noAZF1))

#combine everything
pooleddat1 <- c(AZF1, noAZF1)

n.sim = 10000
nulldist1 =rep(NA, n.sim)

no <- table(data$AZF1)[1]
total <- table(data$AZF1)[1] + table(data$AZF1)[2]

print("Running 10000 permutations for Haplotype D fold change for Rep1")

print(table(data$AZF1))
set.seed(1)
for (i in 1:n.sim) {
  #randomly split data into two groups one group has 73 samples, other has 28
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


ggplot(data) +
  aes(x=AZF1, y=(HAPD_FC), fill = AZF1) +
  geom_boxplot(width = 0.5, outlier.size=.3,outlier.alpha = 0.3)+
  geom_jitter(size=0.3,width = 0.1,alpha=0.3) +
  theme(axis.text.x = element_text(angle = 20)) +
  ylab("HAPD_log2FC") + xlab("AZF1") + theme_bw() + labs(fill="AZF1")+ scale_fill_manual(values=c("slategrey","slateblue1")) + theme_bw() +
  theme(text = element_text(size=20)) + ggtitle(p_value)
