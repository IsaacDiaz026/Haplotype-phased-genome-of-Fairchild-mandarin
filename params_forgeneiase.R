#!/usr/bin/env Rscript


Genes_with_counts <- read.delim("Genes_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes missing estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "overdispersion_param.txt" , quote=F, row.names=F, col.names=F)



#####
print("ACRs")
Genes_with_counts <- read.delim("ACRs_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes missing estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "ACRs_mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "ACRs_overdispersion_param.txt" , quote=F, row.names=F, col.names=F)


#####
print("H3K4me3")
Genes_with_counts <- read.delim("H3K4me3_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes weird estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "H3K4me3_mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "H3K4me3_overdispersion_param.txt" , quote=F, row.names=F, col.names=F)


####
#####
print("H3K36me3")
Genes_with_counts <- read.delim("H3K36me3_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes weird estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "H3K36me3_mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "H3K36me3_overdispersion_param.txt" , quote=F, row.names=F, col.names=F)


####
#####
print("H3K27me3")
Genes_with_counts <- read.delim("H3K27me3_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes weird estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "H3K27me3_mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "H3K27me3_overdispersion_param.txt" , quote=F, row.names=F, col.names=F)

####
#####
print("H3K56ac")
Genes_with_counts <- read.delim("H3K56ac_with_counts.txt", header=FALSE)

#remove any sites with NAs, this fixes weird estimates of fitted and rho
#rho is called "variation" in output,
#there are 3685 sites with NAs

no_NA_genes <- Genes_with_counts[complete.cases(Genes_with_counts),]
colnames(no_NA_genes) <- c("Chr","Start","End","ID","alpha","beta","fitted.values","dispersion","upper.ci","lower.ci","pvalues","error","chr","start","end","y","n")

a <- summary(no_NA_genes$fitted.values)
b <- summary(no_NA_genes$dispersion)

write.table(as.numeric(a[[4]]), "H3K56ac_mean_param.txt" , quote=F, row.names=F, col.names=F)

write.table(as.numeric(b[[4]]), "H3K56ac_overdispersion_param.txt" , quote=F, row.names=F, col.names=F)
