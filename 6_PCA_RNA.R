# PCA RNAseq

setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)

#samples are the rows now
meta = read.table("meta_Nutri_batch.txt", header =T)
rownames(meta) = meta[,1]
meta[,1] = as.factor(meta[,1])
meta[,2] = as.factor(meta[,2])
meta[,3] = as.factor(meta[,3])
meta$Tissue = sub("_.*", "", meta$Sample)
meta$mouseNum = sub(".*_", "", meta$Sample)

metaLiver = meta %>% filter(Tissue == "Liver")
metaQuad = meta %>% filter(Tissue == "Muscle")
metaSAT = meta %>% filter(Tissue == "SAT")
metaVAT = meta %>% filter(Tissue == "VAT") %>% filter(Sample != "VAT_50") #outlier
metaBAT = meta %>% filter(Tissue == "BAT")


TPMs = read.table("Nutrigenetics_TPMs.txt", header = T)
#log scale
TPMs_log = log(TPMs+1)
TPMs_log = TPMs_log[rowSums(TPMs_log != 0) > 0,]


tpms = as.data.frame(t(TPMs_log))
tpms$Strain = meta$Strain
tpms$Diet = meta$Diet
tpms$Sample = meta$Sample
tpms$Tissue = meta$Tissue
tpms$Diet = factor(tpms$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
tpms = na.omit(tpms)
#PCA
PCA <- prcomp((tpms[,1:44842]), center = TRUE,scale. = T)

#pull out the data
PCi<-data.frame(PCA$x,Strain=tpms$Strain, Diet = tpms$Diet, Tissue= tpms$Tissue)
summary(PCA)
eigs <- PCA$sdev^2


colo3 = c("#000430","#d5674b","#74683b","#44c4b8","#cb6683")
colo1 = c("grey","slateblue","lightpink","red3")
colo2 = c("yellow","orange","forestgreen","navy", "purple")



PCA12_tissue = ggplot(PCi,aes(x=PC1,y=PC2,col=Tissue))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo3)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigs[1]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs[2]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

PCA12_diet =ggplot(PCi,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigs[1]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs[2]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

PCA12_strain =ggplot(PCi,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigs[1]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigs[2]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

PCA34_strain =ggplot(PCi,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigs[3]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigs[4]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

PCA34_tissue =ggplot(PCi,aes(x=PC3,y=PC4,col=Tissue))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo3)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigs[3]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigs[4]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")
PCA34_diet =ggplot(PCi,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigs[3]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigs[4]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

PCA56_strain =ggplot(PCi,aes(x=PC5,y=PC6,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC5 (",round(100*eigs[5]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC6 (",round(100*eigs[6]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")
PCA56_diet =ggplot(PCi,aes(x=PC5,y=PC6,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC5 (",round(100*eigs[5]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC6 (",round(100*eigs[6]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")
PCA56_tissue =ggplot(PCi,aes(x=PC5,y=PC6,col=Tissue))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo3)+ #your colors here
  theme_classic()+
  xlab(paste("PC5 (",round(100*eigs[5]/sum(eigs),2) ,"% Variance Explained)"))+
  ylab(paste("PC6 (",round(100*eigs[6]/sum(eigs),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10), legend.position = "none")

pdf(height = 2, width = 2, "rna_pca_strain12.pdf")
PCA12_strain
dev.off()

pdf(height = 2, width = 2, "rna_pca_strain56.pdf")
PCA56_strain
dev.off()

pdf(height = 2, width = 2,"rna_pca_tissue56.pdf")
PCA56_tissue
dev.off()

pdf(height = 2, width = 2, "rna_pca_diet56.pdf")
PCA56_diet
dev.off()

pdf(height = 2, width = 2, "rna_pca_strain34.pdf")
PCA34_strain
dev.off()

pdf(height = 2, width = 2, "rna_pca_tissue34.pdf")
PCA34_tissue
dev.off()

pdf(height = 2, width = 2, "rna_pca_diet34.pdf")
PCA34_diet
dev.off()

pdf(height = 2, width = 2, "rna_pca_tissue12.pdf")
PCA12_tissue
dev.off()

pdf(height = 2, width = 2, "rna_pca_diet12.pdf")
PCA12_diet
dev.off()





tpmsLiver = tpms %>% filter(Tissue == "Liver")
tpmsLiver = tpmsLiver[,colSums(tpmsLiver[,1:44842], na.rm = T) != 0]
PCALiver <- prcomp((tpmsLiver[,1:34930]), center = TRUE,scale. = T)
PCLiver<-data.frame(PCALiver$x,Strain=tpmsLiver$Strain, Diet = tpmsLiver$Diet)
eigsL <- PCALiver$sdev^2

summary(PCALiver)
liverpca_strain = ggplot(PCLiver,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsL[1]/sum(eigsL),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsL[2]/sum(eigsL),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
PC34liverstrain = ggplot(PCLiver,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsL[3]/sum(eigsL),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsL[4]/sum(eigsL),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

liverpca_diet = ggplot(PCLiver,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsL[1]/sum(eigsL),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsL[2]/sum(eigsL),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
liverpca34_diet = ggplot(PCLiver,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsL[3]/sum(eigsL),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsL[4]/sum(eigsL),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

pdf(height = 2, width = 2, "liver pca strain 12.pdf")
liverpca_strain
dev.off()
pdf(height = 2, width = 2, "liver pca diet 12.pdf")
liverpca_diet
dev.off()
pdf(height = 3, width = 3, "liver pca strain 34.pdf")
PC34liverstrain
dev.off()
pdf(height = 3, width = 3, "liver pca diet 34.pdf")
liverpca34_diet
dev.off()





tpmsQuad = tpms %>% filter(Tissue == "Muscle")
tpmsQuad = tpmsQuad[,colSums(tpmsQuad[,1:44842]) != 0]
PCAQuad <- prcomp((tpmsQuad[,1:37274]), center = TRUE,scale. = T)
PCQuad<-data.frame(PCAQuad$x,Strain=tpmsQuad$Strain, Diet = tpmsQuad$Diet)
eigsQ <- PCAQuad$sdev^2

quadpca_strain34 = ggplot(PCQuad,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()++
  xlab(paste("PC3 (",round(100*eigsQ[3]/sum(eigsQ),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsQ[4]/sum(eigsQ),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = 'none')

quadpca_strain = ggplot(PCQuad,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsQ[1]/sum(eigsQ),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsQ[2]/sum(eigsQ),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = 'none')

quadpca_diet = ggplot(PCQuad,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsQ[1]/sum(eigsQ),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsQ[2]/sum(eigsQ),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = 'none')

quadpca_diet34 = ggplot(PCQuad,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsQ[3]/sum(eigsQ),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsQ[4]/sum(eigsQ),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = 'none')

pdf(height = 2, width = 2, "quad pca strain 12.pdf")
quadpca_strain
dev.off()
pdf(height = 2, width = 2, "quad pca diet 12.pdf")
quadpca_diet
dev.off()
pdf(height = 3, width = 3, "quad pca strain 34.pdf")
quadpca_strain34
dev.off()
pdf(height = 3, width = 3, "quad pca diet 34.pdf")
quadpca_diet34
dev.off()


tpmsBAT = tpms %>% filter(Tissue == "BAT")
tpmsBAT = tpmsBAT[,colSums(tpmsBAT[,1:44842]) != 0]
PCABAT <- prcomp((tpmsBAT[,1:37110]), center = TRUE,scale. = T)
PCBAT<-data.frame(PCABAT$x,Strain=tpmsBAT$Strain, Diet = tpmsBAT$Diet)
eigsS <- PCABAT$sdev^2

BATpca_strain =ggplot(PCBAT,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values =colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsB[1]/sum(eigsB),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsB[2]/sum(eigsB),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
BATpca_diet = ggplot(PCBAT,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsB[1]/sum(eigsB),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsB[2]/sum(eigsB),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
BATpca_strain34 =ggplot(PCBAT,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values =colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsB[3]/sum(eigsB),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsB[4]/sum(eigsB),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
BATpca_diet34 = ggplot(PCBAT,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsB[3]/sum(eigsB),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsB[4]/sum(eigsB),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")





pdf(height = 2, width = 2, "BAT pca strain 12.pdf")
BATpca_strain
dev.off()
pdf(height = 2, width = 2, "BAT pca diet 12.pdf")
BATpca_diet
dev.off()
pdf(height = 3, width = 3, "BAT pca strain 34.pdf")
BATpca_strain34
dev.off()
pdf(height = 3, width = 3, "BAT pca diet 34.pdf")
BATpca_diet34
dev.off()


tpmsVAT = tpms %>% filter(Tissue == "VAT")  %>% filter(Sample != "VAT_50")
tpmsVAT = tpmsVAT[,colSums(tpmsVAT[,1:44842]) != 0]
PCAVAT <- prcomp((tpmsVAT[,1:40161]), center = TRUE,scale. = T)
PCVAT<-data.frame(PCAVAT$x,Strain=tpmsVAT$Strain, Diet = tpmsVAT$Diet)
eigsV <- PCAVAT$sdev^2


VATpca_strain = ggplot(PCVAT,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsV[1]/sum(eigsV),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsV[2]/sum(eigsV),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

VATpca_diet = ggplot(PCVAT,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values =colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsV[1]/sum(eigsV),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsV[2]/sum(eigsV),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

VATpca_strain34 = ggplot(PCVAT,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsV[3]/sum(eigsV),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsV[4]/sum(eigsV),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

VATpca_diet34 = ggplot(PCVAT,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values =colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsV[3]/sum(eigsV),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsV[4]/sum(eigsV),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
pdf(height = 2, width = 2, "VAT pca strain 12.pdf")
VATpca_strain
dev.off()
pdf(height = 2, width =2, "VAT pca diet 12.pdf")
VATpca_diet
dev.off()

pdf(height = 3, width = 3, "VAT pca strain 34.pdf")
VATpca_strain34
dev.off()
pdf(height = 3, width = 3, "VAT pca diet 34.pdf")
VATpca_diet34
dev.off()


tpmsSAT = tpms %>% filter(Tissue == "SAT")
tpmsSAT = tpmsSAT[,colSums(tpmsSAT[,1:44842]) != 0]
PCASAT <- prcomp((tpmsSAT[,1:41056]), center = TRUE,scale. = T)
PCSAT<-data.frame(PCASAT$x,Strain=tpmsSAT$Strain, Diet = tpmsSAT$Diet)
eigsS <- PCABAT$sdev^2


SATpca_strain =ggplot(PCSAT,aes(x=PC1,y=PC2,col=Strain))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsS[1]/sum(eigsS),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsS[2]/sum(eigsS),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

SATpca_diet = ggplot(PCSAT,aes(x=PC1,y=PC2,col=Diet))+
  geom_point(size=2,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC1 (",round(100*eigsS[1]/sum(eigsS),2) ,"% Variance Explained)"))+
  ylab(paste("PC2 (",round(100*eigsS[2]/sum(eigsS),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

SATpca_strain34 =ggplot(PCSAT,aes(x=PC3,y=PC4,col=Strain))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo1)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsS[3]/sum(eigsS),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsS[4]/sum(eigsS),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")

SATpca_diet34 = ggplot(PCSAT,aes(x=PC3,y=PC4,col=Diet))+
  geom_point(size=3,alpha=0.75)+ #Size and alpha just for fun
  scale_color_manual(values = colo2)+ #your colors here
  theme_classic()+
  xlab(paste("PC3 (",round(100*eigsS[3]/sum(eigsS),2) ,"% Variance Explained)"))+
  ylab(paste("PC4 (",round(100*eigsS[4]/sum(eigsS),2) ,"% Variance Explained)"))+
  theme(text = element_text(size = 10),legend.position = "none")
pdf(height = 2, width =2, "SAT pca strain 12.pdf")
SATpca_strain
dev.off()
pdf(height = 2, width = 2, "SAT pca diet 12.pdf")
SATpca_diet
dev.off()
pdf(height = 3, width = 3, "SAT pca strain 34.pdf")
SATpca_strain34
dev.off()
pdf(height = 3, width = 3, "SAT pca diet 34.pdf")
SATpca_diet34
dev.off()









