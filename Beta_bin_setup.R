#!/usr/bin/env Rscript

#Fairchild genomic allele counts

#### Beta-Binomial fitting
# In order to obtain the estimates for alpha and beta, the log-likelihood function needs to be defined.  Consult Bolker (In Press) for further information regarding the determination of the likelihood.

lbetabin = function (data, inits) {
 n <- data[,1] # This corresponds to the group_size
 y <-data[,2] # This corresponds to the data on incidence
 alpha <- inits[1] # This corresponds to the initial starting parameter for alpha
 beta <- inits[2] # This corresponds to the initial starting parameter for beta

 # Because optim minimizes a function, the negative log-likelihood is used. Also, the factorial is not necessary, as it does not depend on the parameters of interest
 sum(-lgamma(alpha+y)-lgamma(beta+n-y)+lgamma(alpha+beta+n)+lgamma(alpha)+lgamma(beta)-lgamma(alpha+beta))
}

#### Define initial starting parameters:
# For alpha, use the estimated p from the binomial and for beta, use 1-p
#bglm <- glm(formula = input$y/input$n ~ 1, family = binomial, weights = input$n)

#inits <- c(bglm$fitted.values[1],1-bglm$fitted.values[1])
# ?optim for further information
# The method used was BFGS = Broyden-Fletcher-Goldfarb-Shanno algorithm

#optim.bb <- optim(inits, lbetabin, method='BFGS',hessian=T, data=cbind(input$n,input$y))
#optim.bb
# optim.bb$par provided the estimated alpha and beta,which will be used to obtain the new p
# optim.bb$value is the estimated log-likelihood value, which matches results from Pethybridge using the BBD program of Hughes and Madden



#### An estimate of the mean incidence
####    may be obtained as:
#new.p <- optim.bb$par[1]/(optim.bb$par[1]+optim.bb$par[2])
#new.p

#Estimate of variation
#var.pi = 1/(optim.bb$par[1]+optim.bb$par[2])
#var.pi



#Create function for beta-binomial model of mapping distortion
#MODEL


#This is the correct function
###################
#################
#####################
Fit_BB<-function(input,coord) {

  #Get data
  pos <- split(input$pos,input$chr)

  n <- split(input$n,input$chr)
  y <- split(input$y,input$chr)

  #Results storage
  err <- list()

  gene <- list()

  alpha <- list()
  beta <- list()

  fitted.v <- list()
  variation <- list()
  upper.ci <- list()
  lower.ci <- list()

  pval <- list()

  #Loop through coordinates and subset for gene of interest
  for (g in 1:nrow(coord)) {
    if (g%%1000==0) {print(paste("Working on gene: ",g,sep=""))}

    geneid = paste(coord[g,1],coord[g,2],coord[g,3],sep="_")

    chr <- coord[g,1]

    p <- pos[[chr]]

    dd <- cbind(n[[chr]],y[[chr]])
    #print(dd)
    w <- seq(coord[g,2],coord[g,3],1)

    #Fit model for each gene
    ddt <- matrix(dd[which(p%in%w),],ncol=2)
    inits <- c(0.5,0.5)

    tt <- try(optim(inits, lbetabin, method='L-BFGS-B',data=ddt,control=list(maxit=50)))

    if (inherits(tt,"try-error")) {
      gene[[geneid]] <- geneid
      err[[geneid]] <- 1
      alpha[[geneid]] <- NA
      beta[[geneid]] <- NA
      fitted.v[[geneid]] <- NA
      variation[[geneid]] <- NA
      upper.ci[[geneid]] <- NA
      lower.ci[[geneid]] <- NA
      pval[[geneid]] <- NA
    }

    else if (tt$convergence==0) {
      gene[[geneid]] <- geneid
      err[[geneid]] <- tt$convergence
      alpha[[geneid]] <- tt$par[1]
      #print(alpha[[geneid]])
      beta[[geneid]] <- tt$par[2]
      #print(beta[[geneid]])
      fitted.v[[geneid]] <- tt$par[1]/(tt$par[1]+tt$par[2])
      variation[geneid] <- 1/(tt$par[1]+tt$par[2])
      upper.ci[[geneid]] <- qbeta(1-0.025,shape1=alpha[[geneid]],shape2=beta[[geneid]])
      lower.ci[[geneid]] <- qbeta(0.025,shape1=alpha[[geneid]],shape2=beta[[geneid]])
      pval[[geneid]] <- ks.test(0.5,function(x) pbeta(x,alpha[[geneid]],beta[[geneid]]))$p.value
    }

    else {
      gene[[geneid]] <- geneid
      err[[geneid]] <- tt$convergence
      alpha[[geneid]] <- NA
      beta[[geneid]] <- NA
      fitted.v[[geneid]] <- NA
      variation[geneid] <- NA
      upper.ci[[geneid]] <- NA
      lower.ci[[geneid]] <- NA
      pval[[geneid]] <- NA
    }
  }

  #Output results
  results <- list()
  results[["gene"]] <- unlist(gene,use.names=FALSE)
  results[["alpha"]] <- unlist(alpha,use.names=FALSE)
  results[["beta"]] <- unlist(beta,use.names=FALSE)
  results[["fitted.values"]] <- unlist(fitted.v,use.names=FALSE)
  results[["variation"]] <- unlist(variation,use.names = FALSE)
  results[["upper.ci"]] <- unlist(upper.ci,use.names=FALSE)
  results[["lower.ci"]] <- unlist(lower.ci,use.names=FALSE)
  results[["pvalues"]] <- unlist(pval,use.names=FALSE)
  results[["error"]] <- unlist(err,use.names=FALSE)
  return(results)
}


#TEST MODEL
#IMPORT DAA
Genomic_WASP_ASE_counts <- read.delim("Genomic_WASP_ASE_counts.txt")
data <- Genomic_WASP_ASE_counts
#FILTER DATA
data <- data[data[,8]>10,]

#remove non-canonical chromosomes
data <- data[grep("chr",data$contig),]


input <- data[,c(1,2,7,8)]
colnames(input) <- c("chr","pos","y","n")
input <- lapply(input,gsub, pattern = "chr",replacement = "Chr", fixed = TRUE)
input <- do.call(cbind,input)
input <- as.data.frame(input)
input$y <- as.numeric(input$y)
input$n <- as.numeric(input$n)

Fairchild_gff <- read.delim("genes.gff3", header=FALSE)

Fairchild_gff <- Fairchild_gff[,c(1,4,5)]
genes <- Fairchild_gff[grep("chr",Fairchild_gff$V1),]

genes <- lapply(genes, gsub, pattern = "chr", replacement = "Chr", fixed = TRUE)
genes <- do.call(cbind,genes)
genes <- as.data.frame(genes)
genes$V4 <- as.numeric(genes$V4)
genes$V5 <- as.numeric(genes$V5)


##ACRs
ACRs <- read.delim("ACRs_final.bed", header = FALSE)
ACRs <- ACRs[,c(1:3)]

ACRs <- ACRs[grep("chr",ACRs$V1),]
ACRs$V1 <- gsub("chr","Chr",fixed=TRUE, ACRs$V1)
ACRs$V2 <- as.numeric(ACRs$V2)
ACRs$V3 <- as.numeric(ACRs$V3)

##
H3K4me3 <- read.delim("H3K4me3_peaks.txt.bed",header= FALSE)
H3K4me3 <- H3K4me3[,c(1:3)]

H3K4me3 <- H3K4me3[grep("chr",H3K4me3$V1),]
H3K4me3$V1 <- gsub("chr","Chr",fixed=TRUE, H3K4me3$V1)
H3K4me3$V2 <- as.numeric(H3K4me3$V2)
H3K4me3$V3 <- as.numeric(H3K4me3$V3)

##
H3K36me3 <- read.delim("H3K36me3_peaks.txt.bed",header= FALSE)
H3K36me3 <- H3K36me3[,c(1:3)]

H3K36me3 <- H3K36me3[grep("chr",H3K36me3$V1),]
H3K36me3$V1 <- gsub("chr","Chr",fixed=TRUE, H3K36me3$V1)
H3K36me3$V2 <- as.numeric(H3K36me3$V2)
H3K36me3$V3 <- as.numeric(H3K36me3$V3)

###
H3K27me3 <- read.delim("H3K27me3_peaks.txt.bed",header= FALSE)
H3K27me3 <- H3K27me3[,c(1:3)]

H3K27me3 <- H3K27me3[grep("chr",H3K27me3$V1),]
H3K27me3$V1 <- gsub("chr","Chr",fixed=TRUE, H3K27me3$V1)
H3K27me3$V2 <- as.numeric(H3K27me3$V2)
H3K27me3$V3 <- as.numeric(H3K27me3$V3)

##
H3K56ac <- read.delim("H3K56ac_peaks.txt.bed",header= FALSE)
H3K56ac <- H3K56ac[,c(1:3)]

H3K56ac <- H3K56ac[grep("chr",H3K56ac$V1),]
H3K56ac$V1 <- gsub("chr","Chr",fixed=TRUE, H3K56ac$V1)
H3K56ac$V2 <- as.numeric(H3K56ac$V2)
H3K56ac$V3 <- as.numeric(H3K56ac$V3)



#run model
out <- Fit_BB(input,genes)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, "Genes_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)


#make bed file of input to model
counts_bed <- input[,c(1,2,2,3,4)]

write.table(counts_bed, "Genomic_allele_counts.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


print("Run on ACRs")
#run model
out <- Fit_BB(input,ACRs)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, ACRs_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)

##
print("run on H3K4me3")
#run model
out <- Fit_BB(input,H3K4me3)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, "H3K4me3_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)

###
print("run on H3K6me3")
#run model
out <- Fit_BB(input,H3K36me3)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, "H3K36me3_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)

##
print("run on H3K27me3")
out <- Fit_BB(input,H3K27me3)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, "H3K27me3_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)

print("run on H3K56ac")
out <- Fit_BB(input,H3K56ac)

#maybe I should filter to sites that overlap my genes
out_df <- do.call(cbind,out)
out_df <- as.data.frame(out_df)

sites <- t(as.data.frame(strsplit(out_df$gene,"_")))


out_df$chr <- sites[,c(1)]
out_df$start <- sites[,c(2)]
out_df$end <- sites[,c(3)]

out_df <- out_df[,c(10,11,12,1:9)]
write.table(out_df, "H3K56ac_BB_fitted.bed", col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)
