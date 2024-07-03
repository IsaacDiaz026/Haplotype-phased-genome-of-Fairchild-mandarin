#Combined phenotype data from year 1 (2020-2021) and year 2 (2021-2022)
#Test CV and GP

library(ggplot2)
library(ggforce)
library(ggpubr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(lme4)
library(FactoMineR)
library(factoextra)
library(bestNormalize)

setwd("../PhenotypeDataProcessing/")

#meta data of trees in harvest
GWAS_meta <- read.table("2024-03-29_meta_data.txt",header=T)
##########
#Read in meta data
meta <- read.delim("2020-06-19_CVC_meta_data_simplified.txt",header=TRUE,sep="\t")
#remove leading zeros
meta$CRC <- str_remove(meta$CRC , "^0+")
head(meta)
dim(meta)

metadata_CVC <- read.table("GWAS_sampless_meta.txt",header=T)
###########
#Import packline data for each year

year1 <- read.table("2021-05-18_RAW_packline_data_individual_fruit_2020-2021.txt",header=TRUE,sep="\t")
year2 <- read.table("2022-06-29_RAW_packline_data_individual_fruit_2021-2022.txt",header=TRUE,sep="\t")


#Are headers the same
identical(colnames(year1),colnames(year2))
setdiff(colnames(year1),colnames(year2))
length(colnames(year1))
length(colnames(year2))
cbind(colnames(year1),colnames(year2))

#Reorder
year2 <- year2[,colnames(year1)]
identical(colnames(year1),colnames(year2))

#Add year factor before combining
year1$year <- "2020/2021"
year2$year <- "2021/2022"


#Remove things that we didn't expect to harvest - or with mismatch treeID
table(is.na(year1$Location))
table(is.na(year2$Location))

allyear <- rbind(year1,year2)
allyear <- allyear[!is.na(allyear$Location),]
dim(allyear)


#Filter out fruit that are not grade "A"
#grade D and E represent rough or touching fruit
allyear <- allyear[(allyear$Grade != "E"),]
allyear <- allyear[(allyear$Grade != "D"),]


#######
##Calculate per year mean and variance for each tree/year/traitcombo
#output for gemma - only obvious traits for now
pout <- allyear[,c("CRC","treeid","year","Harvest_date","Type","SizeName",
  "Weight","MajorDiameter","MinorDiameter","Volume","Red","Red.Orange","Dark.Orange",
    "Orange","Orange.Yellow","Yellow.Green","Green","Dark.Green","Light.Scar","Dark.Scar",
      "Yellow","Texture","Overall.Roundness","Stem.angle","Stem.Area","Smoothness","Flatness",
        "Calyx.size","Stem.Size","Rough.skin")]

#Some rows have all zeros, remove these data observations
pout$isfruit <- rowSums(pout[,c(11:27)])
pout <- pout %>% filter(isfruit > 0)

#only keep fruit where calyx is detected
pout <- pout %>% filter(Calyx.size > 0 )
pout <- pout %>% select(-isfruit)

#remove fruit that have no data for everything 
table(pout$MajorDiameter > 0)
table(pout$MinorDiameter > 0)
table(pout$Volume > 0)

#filtering based on fruit density 

pout$Density <- pout$Weight / pout$Volume

#remove fruit with volume less than 10, must be errors
pout <- pout[(pout$Volume > 10),]

density_IQR <- list()

pout$UniqueID <- paste(pout$treeid,pout$year,sep="_")
Unique_IDs <- unique(pout$UniqueID)


density_data <- pout %>% group_by(UniqueID) %>% 
  summarize(densityQ1 = quantile(Density,.25), densityQ2 = quantile(Density, 0.75), densityIQR=IQR(Density))


unique_trees <- list()
tree_density <- list()
number_Fruits <- list()


pout_clean <- pout %>%
  group_by(UniqueID) %>% mutate(IQR=IQR(Density),
                                O_upper = quantile(Density, probs=c(.75), na.rm=F) +1.5*IQR,
                                O_lower=quantile(Density, probs=c(.25), na.rm=F)-1.5*IQR) %>%
  filter(O_lower <= Density & Density <= O_upper)


BEFORE_DENSITY_FILT <- data.frame(table(pout$CRC)) 
AFTER_DENSITY_FILT <- as.data.frame(table(pout_clean$CRC)) 

density_filt <- merge(BEFORE_DENSITY_FILT,AFTER_DENSITY_FILT, by="Var1",all=T)

colnames(density_filt) <- c("CRC","BeforeDensityFilt","AfterFilt")

ggplot(density_filt, aes(x=CRC, y=(BeforeDensityFilt - AfterFilt))) + geom_bar(stat="identity") + coord_flip()

write.table(density_filt,"2024-04-25_densityFiltSummary.txt",col.names = T, row.names=T, sep="\t", quote=F)


#Output mean of reach replicate for modeling
pout_clean <- as.data.frame(pout_clean)
outrep <- vector()
for (i in 7:ncol(pout_clean)) {
  print(colnames(pout_clean)[i])
  print(table(pout_clean[,i] > 0))
  
  tempmean <- tapply(pout_clean[,i],INDEX=list(pout_clean$treeid,pout_clean$year),FUN=mean)
  
  outrep <- cbind(outrep,c(tempmean[,1],tempmean[,2]))
  colnames(outrep)[ncol(outrep)] <- paste(colnames(pout)[i],".mean",sep="")

}


outrep <- cbind(outrep,row.names(outrep))
colnames(outrep)[26] <- c("TreeID")
outrep <- cbind(outrep,c(rep("2020/2021",nrow(outrep)/2),rep("2021/2022",nrow(outrep)/2)))
colnames(outrep)[27] <- c("Year")

outrep <- as.data.frame(outrep)
head(outrep)

outrep <- merge(meta,outrep,by.x="Location",by.y="TreeID",all.y=TRUE)
head(outrep)
dim(outrep)


write.table(pout_clean, "twoYearPackline_indivFruitCLEAN.txt", row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(outrep, "two_year_packline_phenotypes_mean_per_tree.txt",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


#Load destructive trait data, average of 12 fruit per tree

dest <- read.delim("TWO_years_destructive_fruit_quality_phenotypes_per_tree.txt")
head(dest)

dest <- dest[,c(-2,-3,-4,-5,-6,-7,-8)]

head(outrep)
#merge together packline and fruit quality lab data
allpheno <- merge(outrep,dest, by=c("Location", "Year"))

#FILTER TO include only samples we have SNPs called for
allpheno <- allpheno[(allpheno$CRC %in% metadata_CVC$CRC),]

table(allpheno$Year) 
 
table(table(allpheno$CRC))
allpheno$harvest_year <-str_split_i(allpheno$Year, "/",1)

#calculate tree age at time of harvest
allpheno$Year.Planted.1 <- gsub("/26", "",allpheno$Year.Planted.1)
allpheno$Year.Planted.1 <- gsub(".*\\/", "",allpheno$Year.Planted.1)
allpheno$Year.Planted.1 <- as.numeric(allpheno$Year.Planted.1)

#some have year planted formated differently
allpheno$Year.Planted.1 <- ifelse(allpheno$Year.Planted.1 < 100 , (allpheno$Year.Planted.1 + 1900), allpheno$Year.Planted.1)

allpheno$TreeAge <- as.numeric(allpheno$harvest_year) - as.numeric(allpheno$Year.Planted.1)

length(unique(allpheno$CRC))

colors <- allpheno[,c(12:22)]
rownames(colors) <- paste(allpheno$Location,allpheno$Year)

colors <- scale(colors)

pca <- PCA(colors,ncp=10)
pcs <- get_pca_ind(pca)

pcs_for_data <- as.data.frame(pcs[["coord"]][,c(1:4)])

allpheno <- cbind(allpheno,pcs_for_data)

#remove color data
allpheno <- allpheno[,-c(12:22)]
allpheno <- allpheno[,-(15)]



colnames(allpheno)
colnames(allpheno)[1] <- "Location"
colnames(allpheno)[29] <- "Year"

hist(as.numeric(allpheno$Weight.mean))
hist(allpheno$Weight.mean)


write.table(allpheno,"GWAS_samples_phenotypes_pertree.txt", col.names = T, row.names=F,sep = "\t",quote=F)


#####
#Model the effects of scion, rootstock, age, harvest
library(lme4)


colnames(allpheno)
allvar <- vector()
allcoeff <- cbind(unique(allpheno$CRC),unique(allpheno$CRC))
colnames(allcoeff) <- c("CRC","CRC2")

allsd <- vector()

#reorder columns
allpheno <- allpheno[,c(1:7,29,30,8:28,31:34)]
for (p in 10:30) {
  print(p)
  
  ####Normalize data
  tempdata <- allpheno
  if (p%in%c(10:30)) {tempdata[,p] <- boxcox(tempdata[,p]+1)$x.t} #pseudocount
  
  ####Model each trait
  templm <- lmer(tempdata[,p]~(1|CRC) + (1|Rootstock) + (1|TreeAge) + (1|Year), data=tempdata)
  
  #extract fitted values from full model
  fit.temp <-coef(templm)$CRC
  hist(fit.temp[,1],main =colnames(allpheno[p]))
  allcoeff <- merge(allcoeff,fit.temp,by.x="CRC",by.y="row.names",all.x=TRUE)
  print(dim(allcoeff))
  
  #calculate variance components and store sd
  var.comp <- as.data.frame(VarCorr(templm))
  allvar <- cbind(allvar,as.numeric(var.comp[,4])/sum(as.numeric(var.comp[,4])))
  allsd <- c(allsd,sd(fit.temp[,1]))
}

row.names(allvar) <- c("CRC","TreeAge","Rootstock","Year","Residual")
colnames(allvar) <- colnames(allpheno[,10:30])
colnames(allvar) <- gsub(".mean","",colnames(allvar))

allcoeff <- allcoeff[,-2]
colnames(allcoeff) <- c("CRC",colnames(allpheno)[10:30])
head(allcoeff)


h2plot <- rbind(table(cut(allvar[1,],seq(0,1,0.2),include.lowest=T)),table(cut(allvar[2,],seq(0,1,0.2),include.lowest=T)),table(cut(allvar[3,],seq(0,1,0.2),include.lowest=T)),table(cut(allvar[4,],seq(0,1,0.2),include.lowest=T)))
row.names(h2plot) <- row.names(allvar[1:4,])
row.names(h2plot)[1] <- "Scion"

rowSums(h2plot)


dev.off()
pdf("Proportion_of_variance_summary.pdf",width=8,height=6)
barplot(h2plot,beside=TRUE,las=2,ylab="Number of traits",xlab="",legend=TRUE,col=c("darkslategray4","darkgray","gray83","gray97"))
mtext("Proportion of variance explained", side=1, line=5)


barplot(allvar,las=2,legend=TRUE,main = "Scion")
barplot(allvar[2,],las=2,main = "Age") #tree age
barplot(allvar[3,],las=2, main = "Rootstock") #rootstock
barplot(allvar[4,],las=2,main = "Year") #year
dev.off()


#Fit full model to extract variance components including color PCS
pheno_formodel <- tempdata

allcoeff <- cbind(unique(pheno_formodel$CRC),unique(pheno_formodel$CRC))
colnames(allcoeff) <- c("CRC","CRC2")
models <- list()

for (i in c(10:34)) {
  par(mfrow=c(1,1))
  templm <- lmer(pheno_formodel[,i]~(1|CRC) + (1|Rootstock) + (1|TreeAge) + (1|Year), data=pheno_formodel)

  #extract fitted values from full model
  fit.temp <-coef(templm)$CRC
  hist(fit.temp[,1],main =colnames(pheno_formodel[i]),)
  allcoeff <- merge(allcoeff,fit.temp,by.x="CRC",by.y="row.names",all.x=TRUE)
  print(dim(allcoeff))
  }


colnames(allcoeff) <- c("CRC","CRC2",colnames(pheno_formodel[,c(10:34)]))



#reorder phenotype data / filter out samples that aren't in final vcf
samples_in_finalVCF <- read.table("../samples_in_finalVCF.txt", quote="\"", comment.char="")


allcoeffsub <- allcoeff[allcoeff$CRC %in% samples_in_finalVCF$V1,]
notin <- allcoeff$CRC %in% samples_in_finalVCF$V1
NOT <- allcoeff[-(which(notin)),]$CRC
notinREV <- samples_in_finalVCF$V1 %in% allcoeff$CRC
NOTrev <- samples_in_finalVCF[-(which(notinREV)),]


rownames(allcoeffsub) <- allcoeffsub$CRC

x <- as.vector(samples_in_finalVCF$V1)


allcoeff_ordered <- allcoeffsub %>%
  mutate(CRC =  factor(CRC, levels = x)) %>%
  arrange(CRC) 

#remove brix post thaw bc it has 1 NA value and we already have a measurement for brix
allcoeff_ordered <- allcoeff_ordered[,c(-(19))]


write.table(allcoeff_ordered, "pheno_fitted_wHeader.txt" , col.names=T, row.names = F, sep="\t",quote=F)

write.table(allcoeff_ordered, "pheno_fitted_traits.txt" , col.names=F, row.names = F, sep="\t",quote=F)

#export inputs for gwas
for (i in 3:length(allcoeff_ordered)){
  write.table(allcoeff_ordered[,i], paste(names(allcoeff_ordered[i]),"fitted.txt",sep=""),row.names = F, col.names = F, sep="\t",quote=F )
}