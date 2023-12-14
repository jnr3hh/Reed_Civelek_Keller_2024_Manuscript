
### PCA on phenotype data from study 1 and study 2



setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)


load("nutrigenetics_phenotypeData_exp1_outlier_removed.RData")
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")

#load("nutrigenetics_phenotypeData.RData")

#PCA
#remove traits with >10% missing data
#remove samples with any missing data
sums = rep(0, ncol(exp1))
for(i in 4:(ncol(exp1))){
  sums[i] = sum(is.na(exp1[,i]))/nrow(exp1)
}
all1 = na.omit(exp1[,(sums < 0.1)] %>% dplyr::select(!contains("Week_0"))) 

sums = rep(0, ncol(exp2))
for(i in 4:(ncol(exp2))){
  sums[i] = sum(is.na(exp2[,i]))/nrow(exp2)
}
all2 = na.omit(exp2[,(sums < 0.1)] %>% dplyr::select(!contains("Week_0")))

PCA2 <- prcomp((all2[,4:ncol(all2)]), center = TRUE,scale. = TRUE)
PCA1 <- prcomp((all1[,4:ncol(all1)]), center = TRUE,scale. = TRUE)

colo1 = c("grey","slateblue","lightpink","red3")
colo2 = c("yellow","orange","forestgreen","blue")

#pull out the data
PCi1<-data.frame(PCA1$x,Strain=all1$Strain, Diet = all1$Diet)
PCi2<-data.frame(PCA2$x,Strain=all2$Strain, Diet = all2$Diet)
eigs1 <- PCA1$sdev^2
eigs2 <- PCA2$sdev^2

PCA_1_S = ggplot(PCi1,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2, alpha = .75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+theme(legend.position = "none", text = element_text(size = 9))+
  xlab(paste("PC1 (",round(100*eigs1[1]/sum(eigs1),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs1[2]/sum(eigs1),2) ,"% Variance Explained)"))




PCA_2_S =ggplot(PCi2,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2, alpha = .75)+ #Size and alpha just for fun
  scale_color_manual(values =  colo1)+ #your colors here
  theme_classic()+theme(legend.position = "none", text = element_text(size = 9))+
  xlab(paste("PC1 (",round(100*eigs2[1]/sum(eigs2),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs2[2]/sum(eigs2),2) ,"% Variance Explained)"))




PCA_1_D =ggplot(PCi1,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2, alpha = .75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+theme(legend.position = "none", text = element_text(size = 9))+
  xlab(paste("PC1 (",round(100*eigs1[1]/sum(eigs1),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs1[2]/sum(eigs1),2) ,"% Variance Explained)"))



PCA_2_D =ggplot(PCi2,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2, alpha = .75)+ #Size and alpha just for fun
  scale_color_manual(values =  colo2)+ #your colors here
  theme_classic()+theme(legend.position = "none", text = element_text(size = 9))+
  xlab(paste("PC1 (",round(100*eigs2[1]/sum(eigs2),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs2[2]/sum(eigs2),2) ,"% Variance Explained)"))


PCA_1_S
PCA_1_D
PCA_2_S
PCA_2_D

setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest/nutri_figures")
pdf("PCA_1_S_outlier_removed.pdf", height = 2, width = 2)
PCA_1_S
dev.off()
pdf("PCA_1_D_outlier_removed.pdf", height = 2, width = 2)
PCA_1_D
dev.off()
pdf("PCA_2_S_outlier_removed.pdf", height = 2, width = 2)
PCA_2_S
dev.off()
pdf("PCA_2_D_outlier_removed.pdf", height = 2, width = 2)
PCA_2_D
dev.off()




