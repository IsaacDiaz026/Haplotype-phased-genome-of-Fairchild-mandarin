#!/usr/bin/env Rscript

library(caret)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(dplyr)

###2022-11-25 TESTING
#add column to average coverage


All_GENES <- read.delim("ALL_genes_foroverall_expression.txt")

All_GENES[is.na(All_GENES)] <- 0


#load genes determined to be expressed (log(TPM > -3))

expressed <- read.delim("/bigdata/seymourlab/idiaz026/Results/Fairchild/RNA_TPMcounts/ExpressedGenesTPMboth.txt")[,1]



ALL_EXP <- All_GENES[(All_GENES$geneid %in% expressed),]


All_ASE <- read.delim("ASE_genes_for_overall_expression.txt")
All_ASE[is.na(All_ASE)] <- 0


ALL_GENES_forASE <- read.delim("ALL_genes_for_ASE.txt")

ALL_GENES_forASE[is.na(ALL_GENES_forASE)] <- 0

ASE_GENES_forASE <- read.delim("ASE_genes_for_ASE.txt")
ASE_GENES_forASE[is.na(ASE_GENES_forASE)] <- 0


#Use SV haplotypes before ancestry assignment to correct misassigned parentage
SS <-read.delim("ASE_GENES_without_ancestry_assignmentOfSVs.txt")

b <- merge(ASE_GENES_forASE,SS, by=c("geneid"))
b[is.na(b)] <- 0
ggplot(b, aes(x=RNA_BOTH.x,RNA_BOTH.y, color = as.factor(HAP1_promoter_del.x))) + geom_point(size = .5, shape = as.factor(b$HAP1_promoter_del.y))

G <- b %>% filter(HAP1_promoter_del.x == 1 | HAP2_promoter_del == 1)

G$parental_del <- ifelse(G$HAP1_promoter_del.x == 1, "Yes", "No")
G$hap_del <- ifelse(G$HAP1_promoter_del.y == 1 , "Yes", "No")

G$pat_rna <- ifelse(G$parental_del == "Yes", G$RNA1.x, G$RNA2.x)
G$hap_rna <- ifelse(G$ha)


ggplot(G, aes(x=RNA_BOTH.x,RNA_BOTH.y, color = as.factor(HAP1_promoter_del.x))) + geom_point(size = 5, shape = as.factor(G$HAP2_promoter_del)) + ylab("HAP1 /HAP2") + xlab("Maternal / Paternal") + labs(color='SV MATERNAL ALLELE') + theme_bw()


missassigned_dels <- G %>% filter(RNA_BOTH.x != RNA_BOTH.y & HAP1_promoter_del.x == 1)
missassigned_geneids <- missassigned_dels$geneid

toswitch <- ASE_GENES_forASE$geneid %in% missassigned_geneids

rm_toswitch <- ASE_GENES_forASE[-which(toswitch),]
genes_toswitch <- ASE_GENES_forASE[which(toswitch),]

write.table(genes_toswitch$geneid, "/bigdata/seymourlab/idiaz026/Results/Fairchild/TE_SV_ACR_analysis/2023-06-14_genestoswitchDelGT.txt",col.names = F, row.names = F)


#swap response for misassigned genes
genes_toswitch$RNA_BOTH <- genes_toswitch$RNA_BOTH * -1


ASE_GENES_forASE <- rbind(genes_toswitch,rm_toswitch)

ggplot(missassigned_dels, aes(x=RNA_BOTH.x,RNA_BOTH.y, color = as.factor(HAP1_promoter_del.x))) + geom_point(size = 5, shape = as.factor(missassigned_dels$HAP2_promoter_del)) + ylab("HAP1 /HAP2") + xlab("Maternal / Paternal") + labs(color='SV MATERNAL ALLELE') + theme_bw()


table(G$HAP1_promoter_del.x,G$HAP2_promoter_del)
table(G$HAP2_promoter_del)


#outlier filtering for All expressed genes

set.seed(1)

train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)



#predictors
x_df_ase <- ASE_GENES_forASE %>% dplyr::select(-geneid,-RNA1,-RNA2,-RNA_BOTH,-RNA1_cov,-RNA2_cov) %>% scale(center=TRUE, scale=F)

x_df_all <- ALL_GENES_forASE %>% dplyr::select(-geneid,-RNA1,-RNA2,-RNA_BOTH,-RNA1_cov,-RNA2_cov) %>% scale(center=TRUE, scale=F)

x_df_ase_COV <- All_ASE %>% dplyr::select(-geneid,-RNA1,-RNA2,-RNA_BOTH,-RNA_cov_both,-RNA1_cov,-RNA2_cov) %>% scale(center=TRUE, scale=F)

x_df_all_COV <- All_GENES %>% dplyr::select(-geneid,-RNA1,-RNA2,-RNA_BOTH,-RNA_cov_both,-RNA1_cov,-RNA2_cov) %>% scale(center=TRUE, scale =F)

x_df_all_exp_COV <- ALL_EXP %>% dplyr::select(-geneid,-RNA1,-RNA2,-RNA_BOTH,-RNA_cov_both,-RNA1_cov,-RNA2_cov) %>% scale(center=T, scale = F)


#responses
y_df_ase <- ASE_GENES_forASE %>% dplyr::select(RNA_BOTH) %>% scale(center=TRUE, scale=F) %>% as.matrix()

y_df_all_ase <- ALL_GENES_forASE %>% dplyr::select(RNA_BOTH) %>% scale(center=TRUE, scale=F)

y_df_COV_ase <- All_ASE %>% dplyr::select(RNA_cov_both) %>% scale(center=TRUE, scale=F) %>% as.matrix()

y_df_COV_all <- All_GENES %>% dplyr::select(RNA_cov_both) %>% scale(center=TRUE, scale=F) %>% as.matrix()

y_df_COV_all_exp <- ALL_EXP %>% dplyr::select(RNA_cov_both) %>% scale(center=T, scale=F) %>% as.matrix()


#model ase of ase genes
set.seed(1)

elastic_net_model <- train(RNA_BOTH ~ .,
                           data = cbind(x_df_ase,y_df_ase),
                           method = "glmnet",
                           preProcess = c("center"),
                           tuneLength = 25,
                           trControl = train_control)



elastic_net_model


print("ASE ~ ase-genes")
ALPHA <- elastic_net_model[["bestTune"]][["alpha"]]
print("Alpha determined by 5 fold cross-validation")
print(ALPHA)
LAMBDA <- elastic_net_model[["bestTune"]][["lambda"]]
print("Lamda determined by 5 fold cross-validation")
print(LAMBDA)


# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, x_df_ase)
rsq_enet <- cor(y_df_ase, y_hat_enet)^2

new_fit <- glmnet(x_df_ase,y_df_ase ,alpha = ALPHA,lambda = LAMBDA)

res <-glmnet(x_df_ase,y_df_ase ,alpha = 1)
plot(res,label = T)


print(new_fit)

ALPHA <- round(ALPHA,4)
LAMBDA <- round(LAMBDA,4)
R <-round((new_fit[["dev.ratio"]]*100),2)

all_coefs <- coef(new_fit)
only_coefs <- all_coefs[all_coefs != 0]


#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
 geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
 xlab(label="") + ggtitle(paste0("ASE of ASE genes with Lambda = ", LAMBDA, " and Alpha = ", ALPHA)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))

coef_sub <- coef_table[1:12,]


 ggplot(data=coef_sub) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12,L= ", LAMBDA, "A= ", ALPHA, " R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=13))



set.seed(1)
##ALL GENE ASE
elastic_net_model <- train(RNA_BOTH ~ .,
                           data = cbind(x_df_all,y_df_all_ase),
                           method = "glmnet",
                           preProcess = c("center"),
                           tuneLength = 25,
                           trControl = train_control)

elastic_net_model


print("ASE ~ all-genes")
ALPHA <- elastic_net_model[["bestTune"]][["alpha"]]
print("Alpha determined by 5 fold cross-validation")
print(ALPHA)
LAMBDA <- elastic_net_model[["bestTune"]][["lambda"]]
print("Lamda determined by 5 fold cross-validation")
print(LAMBDA)

round(LAMBDA,5)
# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, x_df_all)
rsq_enet <- cor(y_df_all_ase, y_hat_enet)^2

new_fit <- glmnet(x_df_all,y_df_all_ase ,alpha = ALPHA,lambda = LAMBDA)

print(new_fit)

ALPHA <- round(ALPHA,5)
LAMBDA <- round(LAMBDA,5)
R <-round((new_fit[["dev.ratio"]]*100),2)

all_coefs <- coef(new_fit)
only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
 geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
 xlab(label="") + ggtitle(paste0("ASE of all genes with Lambda = ", LAMBDA, " and Alpha = ", ALPHA)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))

coef_sub <- coef_table[1:12,]


ggplot(data=coef_sub) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12,L = ", LAMBDA, " A = ", ALPHA,"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))






set.seed(1)
##ASE_gene_overall expression
elastic_net_model <- train(RNA_cov_both ~ .,
                           data = cbind(x_df_ase_COV,y_df_COV_ase),
                           method = "glmnet",
                           preProcess = c("center"),
                           tuneLength = 25,
                           trControl = train_control)

elastic_net_model


print("Overall expression ~ ASE-genes")
ALPHA <- elastic_net_model[["bestTune"]][["alpha"]]
print("Alpha determined by 5 fold cross-validation")
print(ALPHA)
LAMBDA <- elastic_net_model[["bestTune"]][["lambda"]]
print("Lamda determined by 5 fold cross-validation")
print(LAMBDA)

# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, x_df_ase_COV)
rsq_enet <- cor(y_df_COV_ase, y_hat_enet)^2

new_fit <- glmnet(x_df_ase_COV,y_df_COV_ase ,alpha = ALPHA,lambda = LAMBDA)

print(new_fit)

ALPHA <- round(ALPHA,4)
LAMBDA <- round(LAMBDA,4)
R <-round(new_fit[["dev.ratio"]],2)

all_coefs <- coef(new_fit)
only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
 geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
 xlab(label="") + ggtitle(paste0("Overall Expression ASE genes with Lambda = ", LAMBDA, " and Alpha = ", ALPHA)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))

coef_sub <- coef_table[1:12,]


ggplot(data=coef_sub) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12, L = ", LAMBDA, "A =", ALPHA,"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=13))




set.seed(1)
##All_gene_overall expression
elastic_net_model <- train(RNA_cov_both ~ .,
                           data = cbind(x_df_all_COV,y_df_COV_all),
                           method = "glmnet",
                           preProcess = c("center"),
                           tuneLength = 25,
                           trControl = train_control,verbose=F)

elastic_net_model



print("Overall expression ~ All-genes")
ALPHA <- elastic_net_model[["bestTune"]][["alpha"]]
print("Alpha determined by 5 fold cross-validation")
print(ALPHA)
LAMBDA <- elastic_net_model[["bestTune"]][["lambda"]]
print("Lamda determined by 5 fold cross-validation")
print(LAMBDA)


# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model, x_df_all_COV)
rsq_enet <- cor(y_df_COV_all, y_hat_enet)^2

new_fit <- glmnet(x_df_all_COV,y_df_COV_all ,alpha = ALPHA,lambda = LAMBDA)

print(new_fit)

ALPHA <- round(ALPHA,5)
LAMBDA <- round(LAMBDA,5)
R <-round(new_fit[["dev.ratio"]]*100,2)

all_coefs <- coef(new_fit)
only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
 geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
 xlab(label="") + ggtitle(paste0("Overall Expression All genes with Lambda = ", LAMBDA, " and Alpha = ", ALPHA)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))

coef_sub <- coef_table[1:12,]


ggplot(data=coef_sub) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12,L = ", LAMBDA, "A = ", ALPHA,"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))



set.seed(1)


### PICK LAMBDA
ALPHA <- elastic_net_model[["bestTune"]][["alpha"]]
cvdd <- cv.glmnet(x_df_all_COV,y_df_COV_all ,alpha = ALPHA)

plot(cvdd)
LAMBDA=round(cvdd$lambda.1se,5)



plot(coef(cvdd, s = "lambda.1se"))
index=cvdd[["index"]][[2]]
R <-round(cvdd$glmnet.fit$dev.ratio[index]*100,2)

all_coefs <- coef(cvdd, s = "lambda.1se")

only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("1SE L = ", LAMBDA,", ", "A = ", round(ALPHA,5),"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))


set.seed(1)
##All_EXP gene_overall expression
elastic_net_model1 <- train(RNA_cov_both ~ .,
                           data = cbind(x_df_all_exp_COV,y_df_COV_all_exp),
                           method = "glmnet",
                           preProcess = c("center"),
                           tuneLength = 25,
                           trControl = train_control)

elastic_net_model1


print("Overall expression ~ All-expressed genes")
ALPHA1 <- elastic_net_model1[["bestTune"]][["alpha"]]
print("Alpha determined by 5 fold cross-validation")
print(ALPHA)
LAMBDA1 <- elastic_net_model1[["bestTune"]][["lambda"]]
print("Lamda determined by 5 fold cross-validation")


plot(elastic_net_model1[["results"]]$alpha,elastic_net_model1[["results"]]$RMSE)

plot(elastic_net_model1[["results"]]$lambda,elastic_net_model1[["results"]]$RMSE)

plot(elastic_net_model1[["results"]]$lambda,elastic_net_model1[["results"]]$RMSE,xlim=c(0,0.1))

plot(elastic_net_model1[["results"]]$alpha,elastic_net_model1[["results"]]$RMSE)



# Check multiple R-squared
y_hat_enet <- predict(elastic_net_model1, x_df_all_exp_COV)
rsq_enet <- cor(y_df_COV_all_exp, y_hat_enet)^2

new_fit <- glmnet(x_df_all_exp_COV,y_df_COV_all_exp ,alpha = ALPHA1,lambda=LAMBDA1)

print(new_fit)

ALPHA1 <- round(ALPHA1,5)
LAMBDA1 <- round(LAMBDA1,5)
R <-round((rsq_enet[1]*100),2)

all_coefs <- coef(new_fit)
only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
 geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
 xlab(label="") + ggtitle(paste0("Overall All genes Lambda = ", LAMBDA1, " and Alpha = ", ALPHA1)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))

coef_sub <- coef_table[1:12,]


ggplot(data=coef_sub) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12,L = ", LAMBDA1, "A = ", ALPHA1,"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))


write.table(coef_table, "/bigdata/seymourlab/idiaz026/Data/Fairchild/ASE_modeling_reports/coefs/2023-01-04_all_EXP_gene_expression_coefs.txt", sep = "\t", row.names = T, col.names =T, quote=F)



set.seed(1)

cvdd <- cv.glmnet(x_df_all_exp_COV,y_df_COV_all_exp ,alpha = ALPHA1)

plot(cvdd)


## [1] 0.06284188
plot(coef(cvdd, s = "lambda.1se"))

index=cvdd[["index"]][[2]]
R <-round(cvdd$glmnet.fit$dev.ratio[index]*100,2)

LAMBDA <- cvdd$lambda.1se


all_coefs <- coef(cvdd, s = "lambda.1se")

only_coefs <- all_coefs[all_coefs != 0]

#https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
coef_table<-data.frame(Factor= all_coefs@Dimnames[[1]][all_coefs@i +1], Coefficient = all_coefs@x)

coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = T),]

#lock in order
coef_table$Factor <- factor(coef_table$Factor, levels = coef_table$Factor)

ggplot(data=coef_table) +
  geom_col(aes(x=Factor, y =Coefficient,fill = {Coefficient >0})) + scale_fill_manual(name = 'Coefficient > 0', values = setNames(c('slateblue1','slategrey'),c(T,F)))+
  xlab(label="") + ggtitle(paste0("Top 12,L = ", LAMBDA, "A = ", ALPHA1,"R=",R)) +theme_classic() + theme(axis.text.x = element_text(angle=70, hjust = 1), legend.position = "none") + theme(text=element_text(size=16))
