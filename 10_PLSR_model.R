
setwd("~/")

#BiocManager::install("mixOmics", force = T)
library(mixOmics)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

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
rownames(exp2_half) = exp2_half$MouseNum
VAT_expr$MouseNum = as.numeric(gsub("VAT_","",rownames(VAT_expr)))
rownames(VAT_expr) = VAT_expr$MouseNum
VAT_expr = VAT_expr %>% arrange(MouseNum)
VAT_expr = VAT_expr[,1:417]



X = VAT_expr
Y = exp2_half[,4:28]
rownames(Y)= paste0("m",rownames(Y))
rownames(X)= paste0("m",rownames(X))
Y = Y[!is.na(match(rownames(Y), rownames(X))),]
X = X[!is.na(match(rownames(X), rownames(Y))),]

X = X %>% arrange(rownames(X))
Y = Y %>% arrange(rownames(Y))
X = data.matrix(X)
Y = data.matrix(Y)
Y1 = data.matrix(Y[,5:24])

#Y2 = as.data.frame(Y) %>% dplyr::select(!Week_5_MRI) %>% na.omit()
#X1 = X[!is.na(match(rownames(X), rownames(Y2))),]

pls.nutri <- mixOmics::pls(X, Y, ncomp = 4) # undergo PLS regression

#model validation
#valdcr <- perf(pls.nutri, validation = 'Mfold',
               folds = 10, nrepeat = 10) 
#plot(valdcr, criterion = 'Q2.total')
#list.keepX <- c(seq(10, 40, 5))
# set range of test values for number of variables to use from Y dataframe
#list.keepY <- c(3:9) 
#tune.spls.liver <- tune.spls(X, Y1, ncomp = 2,
#                             test.keepX = list.keepX,
#                             test.keepY = list.keepY,
#                             nrepeat = 1, folds = 10, # use 10 folds
#                             mode = 'regression', measure = 'cor') 
#plot(tune.spls.liver)         # use the correlation measure 

#optimal.keepX <- tune.spls.liver$choice.keepX 
# extract optimal number of variables for Y datafram
#optimal.keepY <- tune.spls.liver$choice.keepY
#optimal.ncomp <-  length(optimal.keepX)



Xvar = pls.nutri[["variates"]][["X"]]
Yvar = pls.nutri[["variates"]][["Y"]]
XYspace = (Xvar+Yvar)/2
normXYspace = as.data.frame(XYspace)
normXYspace$comp1 = (rescale(XYspace[,1], c(-1,1)))
normXYspace$comp2 = (rescale(XYspace[,2], c(-1,1)))
normXYspace$sample = metaVAT$Sample
normXYspace$Strain = metaVAT$Strain
normXYspace$Diet = metaVAT$Diet
normXYspace$X1 = rescale(Xvar[,1], c(-1,1))
normXYspace$Y1 = rescale(Yvar[,1], c(-1,1))
normXYspace$X2 = rescale(Xvar[,2], c(-1,1))
normXYspace$Y2 = rescale(Yvar[,2], c(-1,1))

Xload_ = pls.nutri[["loadings.star"]][[1]][,1:2]
Yload_ = pls.nutri[["loadings.star"]][[2]][,1:2]

Xload = as.data.frame(rescale(Xload_[,1], c(-1,1)))
Yload = as.data.frame(rescale(Yload_[,1], c(-1,1)))
Xload$comp2 = rescale(Xload_[,2], c(-1,1))
Yload$comp2 = rescale(Yload_[,2], c(-1,1))
colnames(Xload) = c("comp1","comp2")
colnames(Yload) = c("comp1","comp2")
Xload$module = pls.nutri[["names"]][["colnames"]][["X"]]
Yload$module = pls.nutri[["names"]][["colnames"]][["Y"]]
Xload = Xload %>% filter(!grepl('Gm|Rik', module)) %>% filter(comp1 <= -.35|comp1 >=.35|comp2>=.35|comp2<=-.35)
Xload$space = "X"
Yload$space = "Y"

YX = rbind(Xload, Yload)

colo1 = c("grey60","slateblue","#fa76ef","red3","black","blue")

colo2 = c(colo2[-5], "black","blue")
normXYspace$sma = paste0(normXYspace$Strain, normXYspace$Diet)
normXYspace$Diet = factor(normXYspace$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan","X","Y"))

bipl = normXYspace %>% ggplot(aes(x = comp1, y = comp2))+
  scale_color_manual(values = colo1 )+
  geom_point(size = 5, aes(color = Strain))+
  geom_text_repel(data = YX, size = 3, aes(x = comp1, y = comp2,label = module, color = space))+
  theme_bw()+
  xlab("Component 1")+ ylab("Component 2")+
  theme(text = element_text(size = 10))

bipl


c(colo2[1], colo2[2],colo2[4],colo2[3],colo2[5],colo2[6])



bipl = normXYspace %>% ggplot(aes(x = comp1, y = comp2))+
  scale_color_manual(values = c(colo2[1], colo2[2],colo2[4],colo2[3],colo2[5],colo2[6]))+
  geom_point(size = 3, aes(color = Diet))+
  geom_text_repel(data = YX, size = 2.8, aes(x = comp1, y = comp2,label = module, color = space), max.overlaps = 20)+
  theme_bw()+
  xlab("Scaled Component 1")+ ylab("Scaled Component 2")+
  theme(text = element_text(size = 10), legend.position = "none")

bipl

pdf("biplot_dietcomppredict genes.pdf", width = 8, height= 5) 
bipl
dev.off()




Xload$external_gene_name = Xload$module


x2 = left_join(Xload, annots2, by = "external_gene_name")
