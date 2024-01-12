#2023-03-16
#checking parentage assignment within phase blocks

library(dplyr)
library(ggplot2)
library(karyoploteR)
library(tidyr)
library(plyr)

fairchild_lr_dels <- read.delim("dels_intPhaseblocks.txt", header=FALSE)

fairchild_lr_dels$block_id <- paste(fairchild_lr_dels$V5,fairchild_lr_dels$V6, sep= "_")

colnames(fairchild_lr_dels) <- c("del_chr","del_start","del_end","del_gt","phs_block_chr","phs_block_start","phs_block_end","block_id")


trio_GTs <- read.delim("phased_wTrio_TrioGTs.txt",header=F)[,c(1:7)]

colnames(trio_GTs) <- c("chr","pos","ref","alt","orl_gt","clem_gt","fair_gt")

trio_GTs <- trio_GTs %>% filter((clem_gt == "0/0" & orl_gt == "1/1") | (orl_gt == "0/0" & clem_gt == "1/1"))

assigned_snps_plus_blockinfo <- read.delim("SNPs_mergedwALLblocks.txt", header=FALSE)

colnames(assigned_snps_plus_blockinfo) <- c("chr","snp_pos","pos_plus1","LR_GT_match","PED_GT","REF","ALT","LR_GT_full","allele_len","phs_block_chr","phs_block_start","phs_block_end")

assigned_snps_plus_blockinfo$block_id <- paste(assigned_snps_plus_blockinfo$phs_block_chr,assigned_snps_plus_blockinfo$phs_block_start,sep="_")

assigned_snps_plus_blockinfo$unassigned <- ifelse(is.na(assigned_snps_plus_blockinfo$LR_GT_match), 1, 0)

#remove indels
assigned_snps_plus_blockinfo <- assigned_snps_plus_blockinfo %>% filter(allele_len ==2)

assigned <- assigned_snps_plus_blockinfo %>% filter(unassigned == 0)
unassigned <- assigned_snps_plus_blockinfo %>% filter(unassigned == 1)

head(unassigned)


assigned$snp_id <- paste(assigned$chr, assigned$snp_pos, sep="_")
trio_GTs$snp_id <- paste(trio_GTs$chr, trio_GTs$pos, sep="_")

intrio_hom <- assigned$snp_id %in% trio_GTs$snp_id
assigned <- assigned[which(intrio_hom),]
assigned <- assigned[,(-15)]


write.table(unassigned[,c(1,2)], "unassigned_snps.txt",col.names = F, row.names=F,sep="\t",quote=F)

print("Assign ancestry of reference allele")


haplotype_1 <- vector()
haplotype_2 <- vector()

hap1_parent <- vector()
hap2_parent <- vector()

hap1 <- data.frame()
hap2 <- data.frame()


assigned$haplotype1 <- ifelse(assigned$LR_GT_full == "0|1", "0","1")
assigned$haplotype2 <- ifelse(assigned$LR_GT_full == "0|1", "1","0")

assigned$hap1_parent <- ifelse((assigned$haplotype1 == "0" & assigned$PED_GT == "0|1") | (assigned$haplotype1 == "1" & assigned$PED_GT == "1|0"), "Dad","Mom")

assigned$hap2_parent <- ifelse((assigned$haplotype2 == "0" & assigned$PED_GT == "0|1") | (assigned$haplotype2 == "1" & assigned$PED_GT == "1|0"), "Dad","Mom")

error <- assigned %>% filter(hap1_parent == hap2_parent)


phase_blocks <- unique(assigned$block_id)
block_tally <- list()
block_tally2 <- list()

block_tally_transpose <- list()
block_tally2_transpose <- list()

meta_list <- list()

mom_hap1 <- vector()
dad_hap1 <- vector()

mom_hap2 <- vector()
dad_hap2 <- vector()




for (i in 1:length(phase_blocks)){
 meta_list[[i]] <- assigned %>% filter(block_id == phase_blocks[i])
 #Add pseudo Mom and Dad counts to each haplotype, this will prevent issues with downstream filtering

 b <- nrow(meta_list[[i]]) +1
 c <- nrow(meta_list[[i]]) +2
 meta_list[[i]][b,] <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","Mom","Dad")
 meta_list[[i]][c,] <- c("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","Dad","Mom")

 block_tally[[i]] <- meta_list[[i]] %>% dplyr::group_by(hap1_parent) %>% tally()
 block_tally2[[i]] <- meta_list[[i]] %>% dplyr::group_by(hap2_parent) %>% tally()

 mom_hap1[i] <- data.frame(block_tally[[i]] %>% filter(hap1_parent == "Mom"))[[2]]

 mom_hap2[i]<- data.frame(block_tally2[[i]] %>% filter(hap2_parent == "Mom"))[[2]]

 dad_hap1[i] <- data.frame(block_tally[[i]] %>% filter(hap1_parent == "Dad"))[[2]]

 dad_hap2[i] <- data.frame(block_tally2[[i]] %>% filter(hap2_parent == "Dad"))[[2]]

 #remove pseudo count
 meta_list[[i]]$mom_hap1 <- mom_hap1[i] - 1
 meta_list[[i]]$mom_hap2 <- mom_hap2[i] - 1

 meta_list[[i]]$dad_hap1 <- dad_hap1[i] - 1
 meta_list[[i]]$dad_hap2 <- dad_hap2[i] - 1

}

meta_list_collapse <- do.call(rbind,meta_list)
meta_list_collapse <- meta_list_collapse[complete.cases(meta_list_collapse),]

meta <- meta_list_collapse %>% dplyr::group_by(block_id) %>% dplyr::summarise(phs_block_chr = phs_block_chr, phase_block_start = phs_block_start, phase_block_end = phs_block_end , Paternal_HAP1 = mean(dad_hap1), Paternal_HAP2 = mean(dad_hap2), Maternal_HAP1 = mean(mom_hap1), Maternal_HAP2 = mean(mom_hap2))

meta <- meta[complete.cases(meta),]

unassign_meta <- unique(unassigned %>% dplyr::group_by(block_id) %>% dplyr::summarise(phs_block_chr = phs_block_chr, phase_block_start = phs_block_start, phase_block_end = phs_block_end ,unassigned = sum(unassigned)))

meta_full <- merge(meta, unassign_meta, by=c("block_id","phs_block_chr","phase_block_start","phase_block_end"), all=T)


meta_full$block_len <- as.numeric(meta_full$phase_block_end) - as.numeric(meta_full$phase_block_start)

meta_full <- unique(meta_full)

meta_full <- meta_full %>%
 mutate_at(c('Maternal_HAP1','Maternal_HAP2',"Paternal_HAP1","Paternal_HAP2","unassigned"), ~replace_na(.,0))

hist(log10(meta_full$block_len))

plot(meta_full$Maternal_HAP1, meta_full$Paternal_HAP2)

meta_full$hap1_percentMaternal <- meta_full$Maternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1 + meta_full$unassigned)

meta_full$hap1_percentPaternal <- meta_full$Paternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1 + meta_full$unassigned)

meta_full$hap2_percentMaternal <- meta_full$Maternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2 + meta_full$unassigned)

meta_full$hap2_percentPaternal <- meta_full$Paternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2 + meta_full$unassigned)


meta_full <- meta_full %>%
 mutate_at(c('hap1_percentMaternal','hap1_percentPaternal',"hap2_percentMaternal","hap2_percentPaternal"), ~replace_na(.,0))

meta_full <- meta_full[complete.cases(meta_full),]

#replace NAs in percantages with 0
library("dplyr")
library("tidyr")




hist((meta_full$hap1_percentMaternal),breaks = 100,xlim = c(0:1), main=c("Including_unassigned"))

hist(meta_full$Maternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1),breaks = 100,xlim = c(0:1), main = c("Maternal / Paternal"))

hist(meta_full$hap1_percentPaternal,breaks = 100,xlim=c(0:1), main=c("Including_unassigned"))

hist(meta_full$Paternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1),breaks = 100,xlim = c(0:1),main = c("Paternal / Maternal"))

hist(meta_full$hap2_percentMaternal,breaks = 100,xlim = c(0:1), main=c("Including_unassigned"))

hist(meta_full$Maternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2),breaks = 100,xlim = c(0:1),main = c("Maternal / Paternal"))

hist(meta_full$hap2_percentPaternal,breaks = 100,xlim = c(0:1), main=c("Including_unassigned"))

hist(meta_full$Paternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2),breaks = 100,xlim = c(0:1),main = c("Paternal / Maternal"))

meta_full$percent_unassigned <- meta_full$unassigned / (((meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2 + meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1)/2) + meta_full$unassigned)


library(MASS)

library(viridis)

get_density <- function(x, y, ...) {
 dens <- MASS::kde2d(x, y, ...)
 ix <- findInterval(x, dens$x)
 iy <- findInterval(y, dens$y)
 ii <- cbind(ix, iy)
 return(dens$z[ii])
}

meta_full$Number_assigned <- meta_full$Paternal_HAP1 + meta_full$Paternal_HAP2 + meta_full$Maternal_HAP1 + meta_full$Maternal_HAP2

meta_full$Number_SNPs <- (meta_full$Paternal_HAP1 + meta_full$Paternal_HAP2 + meta_full$Maternal_HAP1 + meta_full$Maternal_HAP2) / 2  + meta_full$unassigned

meta_full$density <- get_density(meta_full$Number_assigned, meta_full$percent_unassigned, n = 100)

ggplot(meta_full) + geom_point(aes(Number_SNPs, percent_unassigned, color = density)) + scale_color_viridis()


hist(meta_full$percent_unassigned, breaks=100,main = c("Percentage unassigned SNPs per phaseblock"),xlab = c("Percent unassigned SNPs"),freq=F)


fully_assigned <- meta_full %>% filter(percent_unassigned == 0)
fully_unassigned <- meta_full %>% filter(percent_unassigned == 1)

intermed_assign <- meta_full %>% filter(percent_unassigned > 0 & percent_unassigned < 0.4)

sum(fully_assigned$block_len)
sum(fully_unassigned$block_len)
sum(intermed_assign$block_len)


length(unique(assigned_snps_plus_blockinfo$block_id))

ggplot(meta_full, aes(x=Maternal_HAP1, y=Paternal_HAP1)) + geom_point() + theme_bw()

ggplot(meta_full, aes(x=Maternal_HAP2, y=Paternal_HAP2)) + geom_point() + theme_bw()



print("Replace % maternal or paternal with calculation that doesnt include unassigned")

meta_full$hap1_percentMaternal <- meta_full$Maternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1)

meta_full$hap1_percentPaternal <- meta_full$Paternal_HAP1 / (meta_full$Maternal_HAP1 + meta_full$Paternal_HAP1)

meta_full$hap2_percentMaternal <- meta_full$Maternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2)

meta_full$hap2_percentPaternal <- meta_full$Paternal_HAP2 / (meta_full$Maternal_HAP2 + meta_full$Paternal_HAP2)

meta_full <- meta_full %>%
 mutate_at(c('hap1_percentMaternal','hap1_percentPaternal',"hap2_percentMaternal","hap2_percentPaternal"), ~replace_na(.,0))


print("filter phase blocks based on they are all-or-nothing ")

no_recomb <- meta_full %>% filter(hap1_percentMaternal == 1 | hap1_percentPaternal == 1 | hap2_percentPaternal ==1 | hap2_percentMaternal == 1)


no_recomb$type = "non_recomb"

recomb <- meta_full %>% filter(hap1_percentMaternal < 1 & hap1_percentPaternal < 1 & hap2_percentMaternal < 1 & hap2_percentPaternal < 1 )

recomb <- recomb %>% filter_at(vars(Paternal_HAP1,Maternal_HAP1,Maternal_HAP2,Paternal_HAP2), all_vars(!is.na(.)))

recomb$type <- "recomb"

recomb <- recomb %>% filter(percent_unassigned != 1)

#how many phase blocks have assigned snps
length(unique(assigned_snps_plus_blockinfo$block_id))
length(unique(assigned$block_id))




hist((recomb$Maternal_HAP1 / (recomb$Maternal_HAP1 + recomb$Paternal_HAP1)),breaks = 100,xlim = c(0,1))

hist((recomb$Paternal_HAP1 / (recomb$Maternal_HAP1 + recomb$Paternal_HAP1)),breaks = 100,xlim = c(0,1))


hist((recomb$Maternal_HAP2 / (recomb$Maternal_HAP2 + recomb$Paternal_HAP2)),breaks = 100,xlim = c(0,1))

hist((recomb$Paternal_HAP2 / (recomb$Maternal_HAP2 + recomb$Paternal_HAP2)),breaks = 100,xlim = c(0,1))

hist((no_recomb$percent_unassigned),breaks = 100,xlim = c(0,1),xlab=c("Percent unassigned"), main = c("Percent unassigned SNPs in no_recomb blocks"))

hist((recomb$percent_unassigned),breaks = 50,xlim = c(0,1),xlab=c("Percent unassigned"), main = c("Percent unassigned SNPs in recomb blocks"))




#filter snps to only include those that are in non-recomb blocks
in_non_recomb <- assigned_snps_plus_blockinfo$block_id %in% no_recomb$block_id
in_non_recomb <- assigned_snps_plus_blockinfo[which(in_non_recomb),]


#descibe unassigned blocks

unassigned_blocks <- meta_full %>% filter(percent_unassigned == 1)

unassigned_blocks$type = "unassigned"

blocks_classified <-rbind(unassigned_blocks,no_recomb,recomb)

blocks_classified[is.na(blocks_classified)] <- 0

blocks_classified <- blocks_classified  %>% filter(percent_unassigned != 1)

blocks_classified$SNP_density <- (blocks_classified$Paternal_HAP1 + blocks_classified$Maternal_HAP1 + blocks_classified$Paternal_HAP2 + blocks_classified$Maternal_HAP2 + blocks_classified$unassigned) / blocks_classified$block_len

unassigned_blocks$SNP_density <- (unassigned_blocks$Paternal_HAP1 + unassigned_blocks$Maternal_HAP1 + unassigned_blocks$Paternal_HAP2 + unassigned_blocks$Maternal_HAP2 + unassigned_blocks$unassigned) / unassigned_blocks$block_len



ggplot(blocks_classified, aes(x=type,y=block_len,fill=type)) + geom_boxplot() + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=20, size=6, color="red", fill="red")

ggplot(blocks_classified, aes(x=type,y=SNP_density,fill=type)) + geom_boxplot() + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=20, size=6, color="red", fill="red")



non_recomb_snps <- assigned_snps_plus_blockinfo$block_id %in% no_recomb$block_id
non_recomb_SNPs <- assigned_snps_plus_blockinfo[which(non_recomb_snps),]

recomb_snps <- assigned_snps_plus_blockinfo$block_id %in% recomb$block_id
recomb_SNPs <- assigned_snps_plus_blockinfo[which(recomb_snps),] %>% filter(unassigned == 0)

from_unassigned_blocks <- assigned_snps_plus_blockinfo$block_id %in% unassigned_blocks$block_id

unassigned_block_SNPs <- assigned_snps_plus_blockinfo[which(from_unassigned_blocks),]

can_be_assigned <- rbind(non_recomb_SNPs,recomb_SNPs)


#get block info for all SNPs that were unassigned in whatshap
SNPs_for_geno_assignment <- merge(assigned_snps_plus_blockinfo, blocks_classified,by=c("block_id"))

head(SNPs_for_geno_assignment)

#DESIGNATE GENOTYPES SO MATERNAL ALLELE IS ALWAYS IN HAP1
a_mat1 <- SNPs_for_geno_assignment %>% filter(type == "non_recomb" & (hap1_percentMaternal == 1 | hap2_percentPaternal == 1))




a_mat1$geno_assign <- ifelse(a_mat1$LR_GT_full == "0|1", "0|1","1|0")

a_mat2 <- SNPs_for_geno_assignment %>% filter(type == "non_recomb" & (hap2_percentMaternal == 1 | hap1_percentPaternal == 1))

a_mat2$geno_assign <- ifelse(a_mat2$LR_GT_full == "1|0", "0|1","1|0")

b <- SNPs_for_geno_assignment %>% filter(type == "recomb" & !is.na(PED_GT))
b$geno_assign <- ifelse(b$PED_GT == "0|1", "1|0","0|1")

#filter to remove phasblock chr5_2661956 because of extensive haplotype switching

b <- b %>% filter(block_id != "chr5_26691956")
b <- b %>% filter((hap1_percentMaternal > hap1_percentPaternal & LR_GT_full == "0|1" & PED_GT == "1|0") | (hap1_percentMaternal > hap1_percentPaternal & LR_GT_full == "1|0" & PED_GT == "0|1") | (hap1_percentPaternal > hap1_percentMaternal & LR_GT_full == "1|0" & PED_GT == "1|0") | (hap1_percentPaternal > hap1_percentMaternal & LR_GT_full == "0|1" & PED_GT == "0|1") )


final_assigned <- rbind(a_mat1,a_mat2,b)

head(meta_list_collapse)



recomb_block_ids <- unique(recomb_SNPs$block_id)


trio_GTs <- read.delim("phased_wTrio_TrioGTs.txt",header=F)[,c(1:7)]

colnames(trio_GTs) <- c("chr","pos","ref","alt","orl_gt","clem_gt","fair_gt")

trio_GTs <- trio_GTs %>% filter((clem_gt == "0/0" & orl_gt == "1/1") | (orl_gt == "0/0" & clem_gt == "1/1"))


test <- list()


for ( i in 1:length(recomb_block_ids)) {
 test[[i]] <- meta_list_collapse %>% filter(block_id == recomb_block_ids[i])
 test[[i]]$HAP1_mom <- ifelse(test[[i]]$hap1_parent == "Mom","1","0")

 SNP_plot <- ggplot(test[[i]], aes(x=snp_pos, y=HAP1_mom, group=HAP1_mom, color=HAP1_mom)) +
  geom_point(size=0.01) +ggtitle(recomb_block_ids[i])+ theme(axis.line.x =element_blank(),
                                axis.text.x=element_blank(),
                                axis.title.x=element_blank(),
  )
 print(SNP_plot)
}



test_hom_block_type <- merge(SNPs_for_geno_assignment, trio_GTs, by.x=c("chr","snp_pos"),by.y=c("chr","pos") )


recomb_blocks_hom <- unique(test_hom_block_type[,c(3,32)])

table(recomb_blocks_hom$type)

hom_blocks <- merge(SNPs_for_geno_assignment, recomb_blocks_hom, by=c("block_id"))


sum(unique(hom_blocks$Number_SNPs))

sum(unique(SNPs_for_geno_assignment$Number_SNPs))

test_hom <- merge(meta_list_collapse, trio_GTs, by.x=c("chr","snp_pos"),by.y=c("chr","pos") )

recomb_blocks_hom <- unique(test_hom$block_id)
test2 <- list()


for ( i in 1:length(recomb)) {
 test2[[i]] <- meta_list_collapse %>% filter(block_id == recomb$block_id[i])
 test2[[i]] <- merge(test2[[i]], trio_GTs, by.x=c("chr","snp_pos"),by.y=c("chr","pos"))
 test2[[i]]$HAP1_mom <- ifelse(test2[[i]]$hap1_parent == "Mom","1","0")

 SNP_plot <- ggplot(test2[[i]], aes(x=snp_pos, y=HAP1_mom, group=HAP1_mom, color=HAP1_mom)) +
  geom_point(size=0.01) +ggtitle(recomb_block_ids[i])+ theme(axis.line.x =element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.title.x=element_blank(),
  )
 print(SNP_plot)
}





b$haplotype1 <- ifelse(b$LR_GT_full == "0|1", "0","1")
b$haplotype2 <- ifelse(b$LR_GT_full == "0|1", "1","0")

b$hap1_parent <- ifelse((b$haplotype1 == "0" & b$PED_GT == "0|1") | (b$haplotype1 == "1" & b$PED_GT == "1|0"), "Dad","Mom")

b$HAP1_mom <- ifelse(b$hap1_parent == "Mom","1","0")


test3 <- list()
for ( i in 1:length(recomb)) {
 test3[[i]] <- b %>% filter(block_id == recomb$block_id[i])

 SNP_plot <- ggplot(test3[[i]], aes(x=snp_pos, y=HAP1_mom, group=HAP1_mom, color=HAP1_mom)) +
  geom_point(size=0.01) +ggtitle(recomb_block_ids[i])+ theme(axis.line.x =element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.title.x=element_blank(),
  )
 print(SNP_plot)
}

dev.off()

dels_blocks <- merge(fairchild_lr_dels,blocks_classified,by=c("block_id"), all =T)

table(dels_blocks$type)

#do promoters that intersect promoters fall into recomb or non recomb blocks
promoter_del_INT <- read.delim("promoter_del_INT.txt", header=FALSE)[,c(7,8,9,10)]

promoter_dels_blocks <- merge(dels_blocks, promoter_del_INT, by.x=c("del_chr","del_start","del_end"), by.y=c("V7","V8","V9"))

table(promoter_dels_blocks$type)

#most fall into recomb blocks

#Where are unclassified blocks/snps
Fairchild_genome <- read.delim("Fairchild_genome.bed", header=FALSE)[c(1:9),]

colnames(Fairchild_genome) <- c("chr","start","end")

custom.genome <- toGRanges(Fairchild_genome)
kp <- plotKaryotype(genome = custom.genome)

unclass_snps <- toGRanges(unassigned[,c(1,2,3)])
unclass_blocks <- toGRanges(unassigned_blocks[,c(2,3,4)])
class_blocks <- toGRanges(can_be_assigned[,c(10,11,12)])
recomb_blocks <- toGRanges(recomb[,c(2,3,4)])
non_recomb <- toGRanges(no_recomb[,c(2,3,4)])

kp <- plotKaryotype(genome=custom.genome, plot.type=2)



kpPlotRegions(kp, data=recomb_blocks,col="red",data.panel=1)
kpPlotRegions(kp, data=non_recomb,data.panel = 2)




#plot snp ancestry
test <- meta_list_collapse[,c(13,2,3,17)]

test$HAP1_mom <- ifelse(test$hap1_parent == "Mom", "1","0")
test <- test[complete.cases(test),]

custom.genome_blocks <- toGRanges(meta_full[,c(1,3,4)])



kp <- plotKaryotype(genome = custom.genome_blocks,chromosomes = c("chr1_10920416", "chr1_10920859", "chr1_10921933"),plot.type = 2)

chr1_test <- test %>% filter(block_id == "chr1_10920416" | block_id == "chr1_10920859" | block_id == "chr1_10921933" )

chrom <- chr1_test$chr
x=as.numeric(chr1_test$snp_pos)
y=as.numeric(chr1_test$HAP1_mom)



modeling_counts <- list.files("allelic_bias_counts",full.names = T)

modeling_count_dfs <- lapply(modeling_counts,function(i){read.delim(i,head=T)})

dim_before <- vector()
dim_after <- vector()
out <- data.frame()

merged_w_parentage <- list()

for ( i in 1:26){
 print(head(modeling_counts[i],1))
 dim_before[i] <- dim(modeling_count_dfs[[i]])[1]
 merged_w_parentage[[i]] <- merge(modeling_count_dfs[[i]],final_assigned, by=c("chr","snp_pos"))
 dim_after[i] <- dim(merged_w_parentage[[i]])[1]
 out <- cbind(modeling_counts[i],dim_before[i],dim_after[i])

}

out <- cbind(modeling_counts,dim_before,dim_after)
out <- as.data.frame(out)

out$modeling_counts <- gsub("_allelic_counts.txt","",out$modeling_counts)

out$percent_assigned <- (as.numeric(out$dim_after) / as.numeric(out$dim_before)) *100

ggplot(out, aes(x=modeling_counts, y = percent_assigned)) + geom_bar(position="dodge",stat="identity") + coord_flip() + theme_bw()

summary(out$percent_assigned)

head(final_assigned)



#assign GT to dels
dels_mat1 <- dels_blocks %>% filter(type == "non_recomb" & (hap1_percentMaternal == 1 | hap2_percentPaternal == 1))

dels_mat1 <- dels_mat1[complete.cases(dels_mat1),]

dels_mat1$geno_assign <- ifelse(dels_mat1$del_gt == "0|1", "0|1","1|0")


dels_mat2 <- dels_blocks %>% filter(type == "non_recomb" & (hap2_percentMaternal == 1 | hap1_percentPaternal == 1))
dels_mat2 <- dels_mat2[complete.cases(dels_mat2),]

dels_mat2$geno_assign <- ifelse(dels_mat2$del_gt == "1|0", "0|1","1|0")

dels_b <- dels_blocks %>% filter(type == "recomb")
dels_b <- dels_b <- dels_b[complete.cases(dels_b),]
dels_b$geno_assign <- NA


dels_for_assign <- rbind(dels_mat1,dels_mat2,dels_b)


dels_for_assign <-dels_for_assign[,c(2,3,4,5,28)]
colnames(dels_for_assign) <- c("chr","snp_pos","pos_plus1","LR_GT_full","geno_assign")

assigned_for_dels <- final_assigned[,c(2,3,4,9,34)]

toassigndels <- rbind.fill(assigned_for_dels,dels_for_assign)



missing_idx <- which(is.na(toassigndels[,5]))
missing_store <- missing_idx

# function to find the closest 10 variants on the same chromosome, 5 on each side
#do this assigned genotype, then based on dels GT and the nearest geno assigned we actually assign the genotype

closest_ten <- function(idx) {
 chr <- toassigndels[idx, 1]
 subset_data <- toassigndels[toassigndels[, 1] == chr,]
 distances <- abs(subset_data[idx, 2] - subset_data[, 2])
 closest_idx <- c(order(distances[1:(idx-1)])[1:5], order(distances[(idx+1):nrow(subset_data)])[1:5] + idx)
 closest_genotypes <- subset_data[closest_idx, 5]
 return(closest_genotypes)
}

for (idx in missing_idx) {
 closest_genotypes <- closest_ten(idx)
 majority_genotype <- names(which.max(table(closest_genotypes)))
 toassigndels[idx, 5] <- majority_genotype
}


closest_ten_GT <- function(idx) {
 chr <- toassigndels[idx, 1]
 subset_data <- toassigndels[toassigndels[, 1] == chr,]
 distances <- abs(subset_data[idx, 2] - subset_data[, 2])
 closest_idx <- c(order(distances[1:(idx-1)])[1:5], order(distances[(idx+1):nrow(subset_data)])[1:5] + idx)
 closest_genotypes <- subset_data[closest_idx, 4]
 return(closest_genotypes)
}


for (idx in missing_idx) {
 closest_genotypes <- closest_ten_GT(idx)
 majority_genotype <- names(which.max(table(closest_genotypes)))
 toassigndels[idx, 4] <- majority_genotype
}


#now we have the info to assign del genotype
#note the value of in dels to assign are for the majority of variants neighboring each del - not the observed GT of the del

dels_to_assign <- toassigndels[(missing_idx),]
head(dels_to_assign)


dels_hap1_mat <- dels_to_assign %>% filter((LR_GT_full == "0|1" & geno_assign == "0|1") | (LR_GT_full == "1|0" & geno_assign == "1|0"))

#all SNPs near dels have majority of variants being maternal hap1

fin <- merge(dels_to_assign, dels_for_assign, by=c("chr","snp_pos","pos_plus1"))
head(fin)

#if hap1 is always maternal for snps neighboring dels, then based on del gt , assign hap1 allele to maternal_allele
fin$geno_assign <- ifelse(fin$LR_GT_full.y == "1|0","1|0","0|1")

fin <- fin[,c(1,2,3,8)]
colnames(fin) <- c("del_chr","del_start","del_end","geno_assign")
dels_mat1 <- dels_mat1[,c(2,3,4,28)]
dels_mat2 <- dels_mat2[,c(2,3,4,28)]

dels_reassigned <- rbind(fin,dels_mat1,dels_mat2)

# write the output bed file
write.table(dels_reassigned[,c(1,2,3,4)], "delsAncestryAssigned.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


write.table(final_assigned[,c(2,3,3,34)], "SNPs_with_ancestry_assignment.txt", col.names = F, row.names = F, sep = "\t",quote=F)
