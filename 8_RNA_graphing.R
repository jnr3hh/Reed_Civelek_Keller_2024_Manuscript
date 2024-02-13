#
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

Liverc=read.table("liver_deseq_norm_counts.txt",header = T, sep = "\t")
Quadc=read.table("Quad_deseq_norm_counts.txt", header = T, sep = "\t")
BATc=read.table("BAT_deseq_norm_counts.txt", header = T, sep = "\t")
SATc=read.table("SAT_deseq_norm_counts.txt",header = T, sep = "\t")
VATc=read.table("VAT_deseq_norm_counts.txt", header = T, sep = "\t")



BAT <- as.data.frame(t(BATc))
SAT <- as.data.frame(t(SATc))
VAT <- as.data.frame(t(VATc))
Liver <- as.data.frame(t(Liverc))
Quad <- as.data.frame(t(Quadc))

meta = read.table("meta_Nutri_batch.txt", header =T)
rownames(meta) = meta[,1]
meta[,1] = as.factor(meta[,1])
meta[,2] = as.factor(meta[,2])
meta[,3] = as.factor(meta[,3])
meta$Tissue = sub("_.*", "", meta$Sample)

metaVAT = meta %>% filter(Tissue == "VAT") %>% filter(Sample != "VAT_50")
metaQuad = meta %>% filter(Tissue == "Muscle")
metaLiver = meta %>% filter(Tissue == "Liver")
metaBAT = meta %>% filter(Tissue == "BAT")
metaSAT = meta %>% filter(Tissue == "SAT") 

BAT$strain = metaBAT$Strain
BAT$diet = metaBAT$Diet
BAT$strain = factor(BAT$strain)
BAT$diet = factor(BAT$diet)
BAT$batch = metaBAT$batch
BAT$batch = factor(BAT$batch)

Liver$strain = metaLiver$Strain
Liver$diet = metaLiver$Diet
Liver$strain = factor(Liver$strain)
Liver$diet = factor(Liver$diet)
Liver$batch = metaLiver$batch
Liver$batch = factor(Liver$batch)

Quad$strain = metaQuad$Strain
Quad$diet = metaQuad$Diet
Quad$strain = factor(Quad$strain)
Quad$diet = factor(Quad$diet)
Quad$batch = metaQuad$batch
Quad$batch = factor(Quad$batch)

VAT$strain = metaVAT$Strain
VAT$diet = metaVAT$Diet
VAT$strain = factor(VAT$strain)
VAT$diet = factor(VAT$diet)
VAT$batch = metaVAT$batch
VAT$batch = factor(VAT$batch)

SAT$strain = metaSAT$Strain
SAT$diet = metaSAT$Diet
SAT$strain = factor(SAT$strain)
SAT$diet = factor(SAT$diet)
SAT$batch = metaSAT$batch
SAT$batch = factor(SAT$batch)

library(biomaRt)
mart =
  useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl')

annotsB = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(BAT), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsB = as.data.frame(annotsB) %>% filter(external_gene_name != "") %>% na.omit()
annotsB$external_gene_name = make.unique(annotsB$external_gene_name)


annotsQ = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(Quad), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsQ = as.data.frame(annotsQ) %>% filter(external_gene_name != "") %>% na.omit()
annotsQ$external_gene_name = make.unique(annotsQ$external_gene_name)


annotsL = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(Liver), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsL = as.data.frame(annotsL) %>% filter(external_gene_name != "") %>% na.omit()
annotsL$external_gene_name = make.unique(annotsL$external_gene_name)


annotsS = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(SAT), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsS = as.data.frame(annotsS) %>% filter(external_gene_name != "") %>% na.omit()
annotsS$external_gene_name = make.unique(annotsS$external_gene_name)


annotsV = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(VAT), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsV = as.data.frame(annotsV) %>% filter(external_gene_name != "") %>% na.omit()
annotsV$external_gene_name = make.unique(annotsV$external_gene_name)




BAT_EXP = BAT[,!is.na(match(colnames(BAT), annotsB$ensembl_gene_id))]
colnames(BAT_EXP) = annotsB$external_gene_name
BATE = BAT_EXP[(colSums(BAT_EXP[,1:(ncol(BAT_EXP)-3)], na.rm=T) > 0)]
BATs = cbind(BAT[,(ncol(BAT)-2):ncol(BAT)],BATE)

SAT_EXP = SAT[,!is.na(match(colnames(SAT), annotsS$ensembl_gene_id))]
colnames(SAT_EXP) = annotsS$external_gene_name
SATE = SAT_EXP[(colSums(SAT_EXP[,1:(ncol(SAT_EXP)-3)], na.rm=T) > 0)]
SATs = cbind(SAT[,(ncol(SAT)-2):ncol(SAT)],SATE)

VAT_EXP = VAT[,!is.na(match(colnames(VAT), annotsV$ensembl_gene_id))]
colnames(VAT_EXP) = annotsV$external_gene_name
VATE = VAT_EXP[(colSums(VAT_EXP[,1:(ncol(VAT_EXP)-3)], na.rm=T) > 0)]
VATs = cbind(VAT[,(ncol(VAT)-2):ncol(VAT)],VATE)

Liver_EXP = Liver[,!is.na(match(colnames(Liver), annotsL$ensembl_gene_id))]
colnames(Liver_EXP) = annotsL$external_gene_name
LiverE = Liver_EXP[(colSums(Liver_EXP[,1:(ncol(Liver_EXP)-3)], na.rm=T) > 0)]
Livers = cbind(Liver[,(ncol(Liver)-2):ncol(Liver)],LiverE)

Quad_EXP = Quad[,!is.na(match(colnames(Quad), annotsQ$ensembl_gene_id))]
colnames(Quad_EXP) = annotsQ$external_gene_name
QuadE = Quad_EXP[(colSums(Quad_EXP[,1:(ncol(Quad_EXP)-3)], na.rm=T) > 0)]
Quads = cbind(Quad[,(ncol(Quad)-2):ncol(Quad)],QuadE)

BATs$diet = factor(BATs$diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
SATs$diet = factor(SATs$diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
VATs$diet = factor(VATs$diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
Livers$diet = factor(Livers$diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
Quads$diet = factor(Quads$diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))




#single graph RNA expression

fbf_graph = BATs %>% dplyr::select(strain, diet, Fbf1) %>% group_by(strain,diet) %>% summarise(meanPC = mean(Fbf1),sePC = sd(Fbf1)/sqrt(length(Fbf1))) %>% 
  ggplot(aes(x= diet, y = meanPC, group = strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.2, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =strain,color = strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =BATs %>% dplyr::select(strain, diet, Fbf1) %>% group_by(strain,diet), position = position_dodge(width = 0.8),aes(x = diet, y = Fbf1, group= strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("Fbf1 Norm. Expression in BAT")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
fbf_graph
pdf("fbfsat_strain.pdf", width = 3, height=2)
fbf_graph
dev.off()




#facet graph RNA expression (supplemental)


neg_cor_genes_VAT = VATs %>% dplyr::select(strain, diet, Bivm, Klk10, Pcsk4, Unc93a2) %>% pivot_longer(cols=3:6, names_to = "Gene", values_to="Exp") %>%  group_by(strain, diet, Gene)  %>% summarise(meanPC = mean((as.numeric(Exp))),sePC = sd((as.numeric(Exp)))/sqrt(length((as.numeric(Exp))))) %>% ggplot(aes(x= diet, y = meanPC, group = strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.2, fill=colo2[4])+
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = 0.5, aes(color = strain, fill = strain))+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), position = position_dodge(width = 0.8),width = .5)+
  theme_classic()+
  theme(text = element_text(size = 10),axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+
  facet_wrap(~Gene,  scales = "free_y", ncol = 2)+
  xlab("")+
  ylab("Gene Expression")+
  scale_color_manual(values = colo1)+
  scale_fill_manual(values = colo1)

neg_cor_genes_VAT
pdf("neg_cor_genes_VAT.pdf",height = 4, width=6.5)
neg_cor_genes_VAT
dev.off()



# heatmap VAT gene expression


VAT_aov_p = read.table("pvals_VAT_B_75.txt", header= T)
VAT_aov_pv = VAT_aov_p[!is.na(match(rownames(VAT_aov_p), annotsV$ensembl_gene_id)),]
rownames(VAT_aov_pv) = annotsV$external_gene_name
VAT_sigGenes = VAT_aov_pv %>% filter(strain.diet.p.adj <= 0.05)
VAT_expr = VATs %>% dplyr::select(c(rownames(VAT_sigGenes)))

ann1 = metaVAT %>% filter(Sample != "VAT_50")
annn = as.data.frame(ann1$Strain)
colnames(annn) = "Strain"
annn$Diet = ann1$Diet 

annoCol<-list(Diet=c(American=colo2[1], Mediterranean=colo2[2], Vegetarian=colo2[3], Vegan=colo2[4]),Strain =c(AJ=colo1[1], C57=colo1[2], DBA=colo1[3], SJL=colo1[4]))

rownames(annn) = rownames(VAT_expr)
vatheat = pheatmap((VAT_expr), scale = "column", annotation_row = annn, treeheight_row = 5, treeheight_col = 5, fontsize = 10, annotation_colors = annoCol,show_rownames = F)

pdf("VAT heatmap gxe.pdf", height = 3, width = 6.75)
vatheat = pheatmap((VAT_expr), scale = "column", annotation_row = annn, treeheight_row = 5, treeheight_col = 5, fontsize = 10, annotation_colors = annoCol,show_rownames = F)
dev.off()











