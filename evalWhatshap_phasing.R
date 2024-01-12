#2023-02-28

library(dplyr)
library(ggplot2)
#comparing SNPs phased w Longranger, w SNPs phased with whatshap
phased_wTrioSNPs_GTs <- read.delim("phased_wTrioSNPs_GTs.txt", header=FALSE)

linked_readSNPs_GTs <- read.delim("linked_readSNPs_GTs.txt", header=FALSE)

phased_wLR_whatshap <- read.delim("phased_wLR_whatshap_GTs.txt", header=FALSE)


colnames(phased_wTrioSNPs_GTs) <- c("chr","pos","REF","ALT","PED_GT")
colnames(linked_readSNPs_GTs) <- c("chr","pos","REF","ALT","LR_GT")
colnames(phased_wLR_whatshap) <- c("chr","pos","REF","ALT","LR_wh_GT")


phased_wTrioSNPs_GTs$SNP <- (nchar(phased_wTrioSNPs_GTs$REF) + nchar(phased_wTrioSNPs_GTs$ALT))
phased_wTrioSNPs_gts <- phased_wTrioSNPs_GTs %>% filter(SNP == 2)

linked_readSNPs_GTs$SNP <- (nchar(linked_readSNPs_GTs$REF) + nchar(linked_readSNPs_GTs$ALT))
linked_readSNPs_only <- linked_readSNPs_GTs %>% filter(SNP == 2)



shared <- merge(linked_readSNPs_only,phased_wTrioSNPs_gts, by=c("chr","pos"))



match_gts <- shared %>% filter(LR_GT == PED_GT & ALT.x == ALT.y)
table(shared$LR_GT)


head(match_gts)
table(match_gts$LR_GT)
table(match_gts$PED_GT)

b <- shared %>% filter(LR_GT == "0|1" | LR_GT == "1|0")

#what are the parental GTs of unphased variants
trio_GTs <- read.delim("phased_wTrio_TrioGTs.txt",header=F)[,c(1:7)]

colnames(trio_GTs) <- c("chr","pos","ref","alt","orl_gt","clem_gt","fair_gt")
trio_GTs <- trio_GTs[,c(1,2,5,6,7)]

meta <- merge(trio_GTs,shared,by=c("chr","pos"))
phased <- trio_GTs %>% filter(fair_gt=="0|1" | fair_gt == "1|0")

table(phased$clem_gt)
table(phased$orl_gt)

unphased <- meta %>% filter(fair_gt=="0/1")
table(unphased$clem_gt,unphased$orl_gt)
table(unphased$orl_gt)

head(meta)

FINAL_PHASED <- meta %>% filter(PED_GT == "1|0" | PED_GT == "0|1" & ALT.x == ALT.y)


FINAL_PHASED$orl_gt <- ifelse(FINAL_PHASED$orl_gt == "0|1" | FINAL_PHASED$orl_gt == "1|0", "0/1", FINAL_PHASED$orl_gt)

FINAL_PHASED$clem_gt <- ifelse(FINAL_PHASED$clem_gt == "0|1" | FINAL_PHASED$clem_gt == "1|0", "0/1", FINAL_PHASED$clem_gt)

table(FINAL_PHASED$clem_gt, FINAL_PHASED$orl_gt)

FINAL_PHASED <- FINAL_PHASED[,c(1,2,8,12)]



modeling_counts <- list.files("/allelic_bias_counts/",full.names = T)

modeling_count_dfs <- lapply(modeling_counts,function(i){read.delim(i,head=T)})

dim_before <- vector()
dim_after <- vector()
out <- data.frame()

merged_w_parentage <- list()

for ( i in 1:26){
 print(head(modeling_counts[i],1))
 dim_before[i] <- dim(modeling_count_dfs[[i]])[1]
 merged_w_parentage[[i]] <- merge(modeling_count_dfs[[i]],FINAL_PHASED, by.x=c("chr","snp_pos"), by.y=c("chr","pos"))
 dim_after[i] <- dim(merged_w_parentage[[i]])[1]
 out <- cbind(modeling_counts[i],dim_before[i],dim_after[i])



}

out <- cbind(modeling_counts,dim_before,dim_after)
out <- as.data.frame(out)


out$modeling_counts <- gsub("_allelic_counts.txt","",out$modeling_counts)

out$percent_assigned <- (as.numeric(out$dim_after) / as.numeric(out$dim_before)) *100

ggplot(out, aes(x=modeling_counts, y = percent_assigned)) + geom_bar(position="dodge",stat="identity") + coord_flip() + theme_bw()

summary(out$percent_assigned)

forint_w_blocks <- merge(FINAL_PHASED, linked_readSNPs_GTs, all=T, by = c("chr","pos"))

forint_w_blocks$end <- forint_w_blocks$pos + 1

forint_w_blocks <- forint_w_blocks[,c(1,2,9,3,4,5,6,7,8)]

write.table(forint_w_blocks,"SNPs_forintWphaseblocks.bed", col.names = F, row.names = F, sep = "\t",quote=FALSE)
