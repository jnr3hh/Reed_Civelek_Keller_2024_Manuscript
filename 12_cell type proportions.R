library(tidyr)
library(dplyr)
library(ggplot2)



meta = read.table("meta_Nutri_batch.txt", header =T)
rownames(meta) = meta[,1]
meta[,1] = as.factor(meta[,1])
meta[,2] = as.factor(meta[,2])
meta[,3] = as.factor(meta[,3])
meta$Tissue = sub("_.*", "", meta$Sample)
meta$mouseNum = sub(".*_", "", meta$Sample)

meta$Group = factor(paste0(meta$Strain, meta$Diet))

metaLiver = meta %>% filter(Tissue == "Liver")
metaQuad = meta %>% filter(Tissue == "Muscle")
metaSAT = meta %>% filter(Tissue == "SAT")
metaVAT = meta %>% filter(Tissue == "VAT") %>% filter(Sample != "VAT_50")
metaBAT = meta %>% filter(Tissue == "BAT")

SATcell = read.table("CIBERSORTx_satperc_Results.txt", header = T, sep = "\t")
SATcell$Sample = metaSAT$Sample
SATcell$Diet = metaSAT$Diet
SATcell$Strain = metaSAT$Strain

SATcell$Strain = factor(SATcell$Strain)
SATcell$Diet = factor(SATcell$Diet)

VATcell = read.table("CIBERSORTx_vatperc_Results.txt", header = T, sep = "\t")
VATcell$Sample = metaVAT$Sample
VATcell$Diet = metaVAT$Diet
VATcell$Strain = metaVAT$Strain

VATcell$Strain = factor(VATcell$Strain)
VATcell$Diet = factor(VATcell$Diet)

BATcell = read.table("CIBERSORTx_batperc_Results.txt", header = T, sep = "\t")
BATcell$Sample = metaBAT$Sample
BATcell$Diet = metaBAT$Diet
BATcell$Strain = metaBAT$Strain

BATcell$Strain = factor(BATcell$Strain)
BATcell$Diet = factor(BATcell$Diet)


strain.p <- NA
diet.p <- NA
`strain:diet` <- NA

#####

exp1_AOV <- rbind(BATcell, strain.p, diet.p,`strain:diet`)
rownames(exp1_AOV)[(nrow(BATcell)+1):(nrow(BATcell)+3)] <- c("strain.p", "diet.p", "strain_diet")

for(i in 2:15){
  print(i)
  fm <- aov(unlist(exp1_AOV[1:nrow(BATcell),i]) ~ Strain * Diet, data =exp1_AOV[1:nrow(BATcell),])
  exp1_AOV[(nrow(BATcell)+1):(nrow(BATcell)+3), i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}



exp1_AOVres <- as.data.frame(t(exp1_AOV))
exp1_AOVres$strain.p  = as.numeric(exp1_AOVres$strain.p)
exp1_AOVres$diet.p  = as.numeric(exp1_AOVres$diet.p)
exp1_AOVres$strain_diet  = as.numeric(exp1_AOVres$strain_diet)



exp1_AOVres$strain.p.adj <- NA
exp1_AOVres$diet.p.adj <- NA
exp1_AOVres$strain.diet.p.adj <- NA

exp1_AOVres$strain.p.adj[2:(nrow(exp1_AOVres)-7)] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain.p[2:(nrow(exp1_AOVres)-7)])), method = 'fdr')
exp1_AOVres$diet.p.adj[2:(nrow(exp1_AOVres)-7)] <- p.adjust(as.numeric(as.character(exp1_AOVres$diet.p[2:(nrow(exp1_AOVres)-7)])), method = 'fdr')
exp1_AOVres$strain.diet.p.adj[2:(nrow(exp1_AOVres)-7)] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain_diet[2:(nrow(exp1_AOVres)-7)])), method = 'fdr')

#####

twecol = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "gold1",  "palegreen2", "#FDBF6F", "gray70", "maroon", "darkturquoise", "brown" ,"orchid1", "grey30","orange")


SATcelll = SATcell %>% pivot_longer(cols = 2:15, names_to = "celltype", values_to = "absProp")

SATcellll = SATcelll %>% group_by(Diet, Strain, celltype) %>% summarise(meann = mean(absProp))
SATcellll$Diet = factor(SATcellll$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
SATcellll$Strain = factor(SATcellll$Strain, levels = c("AJ","C57","DBA","SJL"))

pdf("SATcelltypes.pdf",width=3, height = 3)
SATcellll %>% ggplot(aes(fill = celltype, y = meann, x = Strain))+
  geom_bar(position = "stack", stat = "identity")+facet_wrap(~Diet)+ theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ylab("Absolute Cell Type Proportion in SAT")+scale_fill_manual(values = twecol)
dev.off()

VATcelll = VATcell %>% pivot_longer(cols = 2:15, names_to = "celltype", values_to = "absProp")

VATcellll = VATcelll %>% group_by(Diet, Strain, celltype) %>% summarise(meann = mean(absProp))
VATcellll$Diet = factor(VATcellll$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
VATcellll$Strain = factor(VATcellll$Strain, levels = c("AJ","C57","DBA","SJL"))

pdf("VATcelltypes.pdf",width=3, height = 3)
VATcellll %>% ggplot(aes(fill = celltype, y = meann, x = Strain))+
  geom_bar(position = "stack", stat = "identity")+facet_wrap(~Diet)+ theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ylab("Absolute Cell Type Proportion in VAT")+scale_fill_manual(values = twecol)
dev.off()

BATcelll = BATcell %>% pivot_longer(cols = 2:15, names_to = "celltype", values_to = "absProp")
BATcellll = BATcelll %>% group_by(Diet, Strain, celltype) %>% summarise(meann = mean(absProp))
BATcellll$Diet = factor(BATcellll$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
BATcellll$Strain = factor(BATcellll$Strain, levels = c("AJ","C57","DBA","SJL"))


pdf("BATcelltypes.pdf",width=3, height = 3)
BATcellll %>% ggplot(aes(fill = celltype, y = meann, x = Strain))+
  geom_bar(position = "stack", stat = "identity")+facet_wrap(~Diet)+ theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ylab("Absolute Cell Type Proportion in BAT")+scale_fill_manual(values = twecol)
dev.off()

quadcelll = quadcell %>% pivot_longer(cols = 2:13, names_to = "celltype", values_to = "absProp")

pdf("quadcelltypes.pdf",width=3, height = 3)
quadcelll %>% ggplot(aes(fill = celltype, y = absProp, x = Strain))+
  geom_bar(position = "stack", stat = "identity")+facet_wrap(~Diet)+ theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ylab("Absolute Cell Type Proportion in Muscle")+scale_fill_manual(values = twecol)
dev.off()

livercelll = livercell %>% pivot_longer(cols = 2:26, names_to = "celltype", values_to = "absProp")

pdf("livrecelltypes.pdf",width=3, height = 3)
livercelll %>% ggplot(aes(fill = celltype, y = absProp, x = Strain))+
  geom_bar(position = "stack", stat = "identity")+facet_wrap(~Diet)+ theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ylab("Absolute Cell Type Proportion in Liver")+scale_fill_manual(values = twecol)
dev.off()





VATcell = VATcell %>% dplyr::select(Diet,Strain,Monocyte,M0.Macrophage,M1.Macrophage,M2.Macrophage)
VATcell$Sample = metaVAT$Sample
VATcell$Diet = metaVAT$Diet
VATcell$Strain = metaVAT$Strain

VATcell$Strain = factor(VATcell$Strain)
VATcell$Diet = factor(VATcell$Diet)

VATcelll = VATcell %>% pivot_longer(cols = 3:6, names_to = "celltype", values_to = "absProp")
VATcelll$Diet = factor(VATcelll$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

VATcelll$celltype = factor(VATcelll$celltype, levels = c("Monocyte","M0.Macrophage","M1.Macrophage","M2.Macrophage"))
pdf("VATcelltypesM.pdf",width=3, height = 3)
VATcelll %>% group_by(Strain, Diet, celltype)  %>% summarise(meanPC = mean((as.numeric(absProp))),sePC = sd((as.numeric(absProp)))/sqrt(length((as.numeric(absProp))))) %>% ggplot(aes(x= Strain, y = meanPC, fill = celltype))+
  geom_bar(position = position_dodge(width = 0.9), stat = "identity")+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), position = position_dodge(width = 0.8),width = .5)+
  facet_wrap(~Diet)+ 
  theme_classic()+ 
  theme(legend.position = "none",text = element_text(size = 10))+
  ylab("Absolute Cell Type Proportion in VAT")+  
  scale_fill_manual(values = c( twecol[8],twecol[4:6]))
dev.off()




VATcell = VATcell %>% dplyr::select(Diet,Strain,T_CD4_average, T_CD8_average, Th_average)
VATcell$Sample = metaVAT$Sample
VATcell$Diet = metaVAT$Diet
VATcell$Strain = metaVAT$Strain

VATcell$Strain = factor(VATcell$Strain)
VATcell$Diet = factor(VATcell$Diet)

VATcelll = VATcell %>% pivot_longer(cols = 3:5, names_to = "celltype", values_to = "absProp")
VATcelll$Diet = factor(VATcelll$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

#VATcelll$celltype = factor(VATcelll$celltype, levels = c("Monocyte","M0.Macrophage","M1.Macrophage","M2.Macrophage"))
pdf("VATcelltypest.pdf",width=3, height = 3)
VATcelll %>% group_by(Strain, Diet, celltype)  %>% summarise(meanPC = mean((as.numeric(absProp))),sePC = sd((as.numeric(absProp)))/sqrt(length((as.numeric(absProp))))) %>% ggplot(aes(x= Strain, y = meanPC, fill = celltype))+
  geom_bar(position = position_dodge(width = 0.9), stat = "identity")+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), position = position_dodge(width = 0.8),width = .5)+
  facet_wrap(~Diet)+ 
  theme_classic()+ 
  theme(legend.position = "none",text = element_text(size = 10))+
  ylab("Absolute Cell Type Proportion in VAT")+
  scale_fill_manual(values = twecol[12:14])
dev.off()
