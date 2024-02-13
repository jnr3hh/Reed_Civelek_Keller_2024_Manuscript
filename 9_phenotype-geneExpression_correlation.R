# Correlations between phenotypes exp 2 and RNA seq


setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)

#load data
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")
load("nutrigenetics_phenotypeData.RData")


VAT_aov_p = read.table("pvals_VAT_B_75.txt", header= T)
VAT_aov_pv = VAT_aov_p[!is.na(match(rownames(VAT_aov_p), annotsV$ensembl_gene_id)),]
rownames(VAT_aov_pv) = annotsV$external_gene_name
VAT_sigGenes = VAT_aov_pv %>% filter(strain.diet.p.adj <= 0.05)
VAT_expr = VATs %>% dplyr::select(c(rownames(VAT_sigGenes)))



# Remove Duplicate Column Names
VAT_expr = VAT_expr[!duplicated(colnames(VAT_expr))]

#select samples
exp2_half = exp2[!is.na(match(exp2$MouseNum, gsub("VAT_","",rownames(VAT_expr)))),]
VAT_expr$MouseNum = as.numeric(gsub("VAT_","",rownames(VAT_expr)))
rownames(VAT_expr) = VAT_expr$MouseNum
VAT_expr = VAT_expr %>% arrange(MouseNum)
VAT_expr = VAT_expr[,1:417]






cor = data.frame(cor = 0,pval = 0,trait = 0,gene = 0)
for (i in 4:ncol(exp2_half)){
  for (j in 1:ncol(VAT_expr)){
    cor1 = data.frame(cor = 0,pval = 0,trait = 0,gene = 0)
    print(i)
    print(j)
    l = cor.test(as.numeric(VAT_expr[,j]), as.numeric(exp2_half[,i]))
    cor1$cor = l$estimate 
    cor1$pval = l$p.value
    cor1$gene = colnames(VAT_expr)[j]
    cor1$trait = colnames(exp2_half)[i]
    cor = rbind(cor,cor1)
  }
}
cor$padj = p.adjust(cor$pval, method = "fdr", n = length(cor$pval))

cor_sig = cor %>% filter(padj <= 0.05)


write.table(cor_sig, "sig_correlations_gene_phen_outlier in.txt", quote = F, sep = "\t")
