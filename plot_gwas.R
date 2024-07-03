#!/usr/bin/Rscript


library(qqman)
library(pacman)
p_load(plyr, dplyr, ggplot2, cowplot)
library(gtools)
library(stringr)
library(ggrastr)
library(ggplot2)
library(patchwork)
library(data.table)

setwd("../Fitted_values/")

GWAS_results <- list.files("/output", pattern = ".assoc.txt", full.names = TRUE)

All_results <- lapply(GWAS_results, function(i){
  read.delim(i,header=TRUE)
})


for (i in seq_along(All_results)) {
  All_results[[i]]$file <- gsub("_gwas.relax.assoc.txt","",str_split_i(GWAS_results[i],"output/",2))
}


sig_snp_pos <- list()
sig_snp_info <- list()

colour1="slateblue1"
colour2="slategrey"

for (i in All_results) {
  filenames <- i$file[1]
  print(filenames)
  i <- i[complete.cases(i),]
  thresh= -log10(0.05/nrow(i))
  i$chr <- as.numeric(gsub("chr","",i$chr))
  i$pos2=i$ps

  for(chr in 2:9){
    tmp=chr-1
    tmp2=max(as.numeric(i[i$chr==tmp,"pos2"]),na.rm=T)
    i[i$chr==chr,"pos2"]=i[i$chr==chr,"ps"]+tmp2
  } 
  df <- i

  df2=data.frame(chr=1:9,center=NA,stringsAsFactors=F)
  for(chr in 1:9){
    start=min(df[df$chr==chr,"pos2"],na.rm=T)
    end=max(df[df$chr==chr,"pos2"],na.rm=T)
    center=mean(c(start,end))
    df2[df2$chr==chr,"center"]=center
  }
  

  g1=ggplot(df,aes(x=pos2,y=-log10(p_wald),colour=as.character(chr)))+
  geom_hline(yintercept=thresh,linetype="dashed")+
  geom_point_rast(size=0.5)+
  scale_colour_manual(values=c("1"=colour1,"2"=colour2,"3"=colour1,"4"=colour2,"5"=colour1,"6"=colour2,"7"=colour1,"8"=colour2,"9"=colour1))+ # not cool though
  theme_classic()+theme(legend.position="NONE")+
  scale_x_continuous(breaks=df2$center,labels=c(1:9))+
  xlab("Chromosome")+
  ylab(expression(paste(-log[10],"(",italic(P), " value)")))+
  geom_hline(yintercept=thresh,linetype="dashed") + ggtitle(paste(filenames))
  ggsave(paste("Manhattan_plots/Raster/",filenames,".pdf",sep=""),g1,width=9,height=3)


   
}

#Export significant SNPs
top_snp_pos <- list() 

for (i in (1:length(All_results))) {
  filenames <- All_results[[i]]$file[1]
  print(filenames)
  All_results[[i]] <- All_results[[i]][complete.cases(All_results[[i]]),]
  thresh = 0.05 / nrow(All_results[[i]])
  top_snp_pos[[i]] <- All_results[[i]] %>% top_n(n=-1000,wt=p_wald) %>% select(chr, ps,beta, se,p_wald,file)
  sig_snp_info[[i]] <- All_results[[i]] %>% filter(p_wald <= thresh ) %>% group_by(chr) %>% select(chr, ps,beta, se,p_wald,file) 
  sig_snp_pos[[i]] <- sig_snp_info[[i]] %>% select(chr, ps)
  
}

top_snps <- do.call(rbind,top_snp_pos)
sig_snp_full <- do.call(rbind, sig_snp_info)
sig_snps <- do.call(rbind,sig_snp_pos)

print("export")
write.table(top_snps, "TOP1000_snps_PERTRAIT_metadata.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(sig_snp_full, "significant_snps_metadata.txt",sep="\t",quote = F, col.names = T, row.names = F)

write.table(sig_snps,"significant_snps_pos.txt",sep="\t",quote = F, col.names = F, row.names = F)

