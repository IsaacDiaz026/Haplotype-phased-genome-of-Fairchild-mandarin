
library(DESeq2)
library(ggplot2)
library(stringr)
library(dplyr)

setwd("../counts_data/")


#Differential expression analysis from STAR output
meta <- read.table("FruitDev_metadata.txt",header=T)

# first read in the matrix
count_matrix <- read.delim("FruitDev_counts.txt", header=T, sep="\t", row.names=1)

count_matrix <- count_matrix[,c(-19)]

meta$Timepoint <- factor(meta$Timepoint)


#make deseq obj
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = meta,
                                         design = ~ Timepoint)

se_star_matrix <- se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]
nrow(se_star_matrix)
se_star2 <- DESeq(se_star_matrix,)


norm_counts <- log2(counts(se_star2, normalized = TRUE)+1)
vsd <- vst(se_star2)

resultsNames(se_star2)

# extract results using timepoint 1 as the baseline for comparison
de_2_1 <- data.frame(results(object = se_star2, 
              name="Timepoint_2_vs_1"))

de_3_1 <- data.frame(results(object = se_star2, 
                             name="Timepoint_3_vs_1"))

de_4_1 <- data.frame(results(object = se_star2, 
                             name="Timepoint_4_vs_1"))


de_5_1 <- data.frame(results(object = se_star2, 
                             name="Timepoint_5_vs_1"))

de_6_1 <- data.frame(results(object = se_star2, 
                             name="Timepoint_6_vs_1"))


de_list <- list(de_2_1, de_3_1,de_4_1,de_5_1,de_6_1)
names_for_columns <- list()

timepoint <- c("2","3","4","5","6")


for (i in 1:5) {
  de_list[[i]]$GeneID <- rownames(de_list[[i]])
  de_list[[i]] <- de_list[[i]] %>% select(GeneID, log2FoldChange, padj)
  de_list[[i]]$comparison <- timepoint[i]
    }

comparisons <- do.call(rbind, de_list)
write.table(comparisons, "Deseq_resultsFruit_timecourse.txt",col.names = T, row.names=T, sep="\t",quote=F )


plotCounts(se_star2, gene="FAIR_000390", intgroup="Timepoint",normalized = T)

KC1 <- as.data.frame(plotCounts(se_star2, gene="FAIR_000390", intgroup="Timepoint",returnData = T,normalized=T))


# Plot the MOV10 normalized counts, using the samplenames (rownames(d) as labels)
a <- round(de_2_1 %>% filter(rownames(de_2_1)=="FAIR_000390") %>% select(padj),digits = 6)
b <- round(de_3_1 %>% filter(rownames(de_3_1)=="FAIR_000390") %>% select(padj),digits = 6)
c <- round(de_4_1 %>% filter(rownames(de_4_1)=="FAIR_000390") %>% select(padj),digits = 6)
d <- round(de_5_1 %>% filter(rownames(de_5_1)=="FAIR_000390") %>% select(padj),digits = 6)
e <- round(de_6_1 %>% filter(rownames(de_6_1)=="FAIR_000390") %>% select(padj),digits=6)


pdf("KC1_fruit_development.pdf")
ggplot(KC1, aes(x = Timepoint, y = count)) +geom_boxplot() + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) + 
  ggtitle(paste("KC1 ",a,b,c,d,e))+ theme_classic()

dev.off()
  


## GENOTYPE COMPARISONS
meta <- read.table("species_fruit_metadata.txt",header=T)

# first read in the matrix
count_matrix <- read.delim("species_fruit_counts.txt", header=T, sep="\t",row.names=1)

#filter to include mand pummelo and swo 
tokeep <- c("Csinensis","Cmaxima","Creticulata")
meta <- meta[meta$Genotype %in% tokeep,]

count_matrix <- count_matrix[,(colnames(count_matrix) %in% meta$SampleName)]


se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = meta,
                                         design = ~Genotype)

se_star_matrix <- se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]
nrow(se_star_matrix)


se_star2 <- DESeq(se_star_matrix,)
norm_counts <- log2(counts(se_star2, normalized = TRUE)+1)

resultsNames(se_star2)

#Contast to pummelo (C.maxima)
de_ret_max <- data.frame(results(object = se_star2, 
                             contrast=c("Genotype","Creticulata","Cmaxima")))

de_sin_max <- data.frame(results(object = se_star2, 
                                       contrast=c("Genotype","Csinensis","Cmaxima")))
de_ret_max$geneID <- rownames(de_ret_max)
de_sin_max$geneID <- rownames(de_sin_max)

ret_max <- de_2_1 %>% select(geneID, log2FoldChange,padj)
sin_max <- de_3_1 %>% select(geneID, log2FoldChange,padj)

ret_max$Compare <- "reticulata_vs_maxima"
sin_max$Compare <- "sinensis_vs_maxima"

geno_comparisons <- rbind(ret_max,sin_max)

write.table(geno_comparisons, "Deseq_results_species_comparisons.txt", col.names=T, row.names = F, sep="\t",quote=F)


plotCounts(se_star2, gene="FAIR_000390", intgroup="Genotype",normalized = T)

KC1 <- as.data.frame(plotCounts(se_star2, gene="FAIR_000390", intgroup="Genotype",returnData = T,normalized=T))

head(KC1)


a <- round(de_ret_max %>% filter(rownames(de_ret_max)=="FAIR_000390") %>% select(padj),digits = 6)
b <- round(de_sin_max %>% filter(rownames(de_sin_max)=="FAIR_000390") %>% select(padj),digits = 6)

pdf("KAT3_GENOTYPE.pdf")
ggplot(KC1, aes(x = reorder(Genotype, count), y = count, color = Genotype)) + 
  theme_classic() +geom_boxplot() + geom_point() + ggtitle(paste("KC1 by Genotype",a,b)) + 
  ylab("Normalized Counts") + scale_color_manual(values=c("slategray","slateblue4","slateblue1"))

dev.off()

#load individual fruit data packline
twoYearPackline_indivFruitCLEAN <- read.delim("../PhenotypeDataProcessing/twoYearPackline_indivFruitCLEAN.txt")
RAW_packline_data_individual_fruit_2020.2021 <- read.delim("../PhenotypeDataProcessing/2021-05-18_RAW_packline_data_individual_fruit_2020-2021.txt")

RAW_packline_data_individual_fruit_2020.2021$Year <- "2020/2021"

RAW_packline_data_individual_fruit_2021.2022 <- read.delim("2022-06-29_RAW_packline_data_individual_fruit_2021-2022.txt")
RAW_packline_data_individual_fruit_2021.2022$Year <- "2021/2022"


rawpheno <- rbind(RAW_packline_data_individual_fruit_2020.2021,RAW_packline_data_individual_fruit_2021.2022)


inpheno <- unique(twoYearPackline_indivFruitCLEAN$CRC)

#select accessions that were used in RNA seq comparison
ofinterest <- c("3224","3868","4003")
toplot <- rawpheno[rawpheno$CRC %in% ofinterest,]
toplot <- toplot %>% filter(MajorDiameter > 0, Grade=="A" )

hist(toplot$Major.Diameter..mm.)


table(inpheno %in% ofinterest)

meta <- data.frame()


toplot$Major.Diameter..mm. <- as.numeric(toplot$Major.Diameter..mm.)

ggplot(toplot, aes(x = reorder(CRC, -Major.Diameter..mm.), y=Major.Diameter..mm. )) + geom_jitter() + geom_boxplot()



#Raw destructives per fruit 
setwd("../PhenotypeDataProcessing")

rawdest <- list.files(".", pattern="*FQ_data",full.names = T)
rawdest <- rawdest[-4]


rawdest_df <- lapply(rawdest, function(i){
  read.delim(i,head=T)
})

rawdest_df[[1]]$Year = "2020/2021"
rawdest_df[[2]]$Year = "2020/2021"
rawdest_df[[3]]$Year = "2020/2021"

rawdest_df[[4]]$Year = "2021/2022"
rawdest_df[[5]]$Year = "2021/2022"
rawdest_df[[6]]$Year = "2021/2022"


allrawdest <- do.call(rbind,rawdest_df)

allrawdest <- merge(allrawdest, toplot, by.x="TreeID",by.y="Location")

allrawdest <- allrawdest %>% select(TreeID, CRC,Variety, Year.x, RindFruit1, RindFruit2, RindFruit3,RindFruit4,RindFruit5,RindFruit6,RindFruit7,RindFruit8,RindFruit9,RindFruit10,RindFruit11,RindFruit12)


#concatenate rind measurements into 1 column
fruits <- list()
for (i in 5:16){
  fruits[[i]] <- allrawdest[,c(1,2,3,4,i)]
  colnames(fruits[[i]]) <- c("TreeID","CRC","Variety","Year.x","Rind")
}

dest <- do.call(rbind,fruits)

dest$Rind <- as.numeric(as.character(dest$Rind))
hist(dest$Rind)
dest <- unique(dest)

table(dest$CRC)
table(toplot$CRC)


pdf("Peel_thickness.pdf")
ggplot(dest, aes(x = CRC, y= Rind,color=CRC)) + geom_boxplot() + geom_jitter(width=0.1) + scale_color_manual(values=c("slategray","slateblue1","slateblue4")) + theme_classic() + ggtitle("n=46,45,23, p < 0.001, DunnTest")

ggplot(toplot, aes(x = reorder(CRC, -Major.Diameter..mm.), y=Major.Diameter..mm. , color=CRC)) + geom_jitter(width=0.1) + geom_boxplot() + scale_color_manual(values=c("slategray","slateblue1","slateblue4")) + theme_classic()

dev.off()

tests <- kwAllPairsDunnTest(dest$Rind,g=as.factor(dest$Variety))
tests2 <- kwAllPairsDunnTest(toplot$MajorDiameter,g=as.factor(toplot$CRC))



