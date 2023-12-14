# ANOVA RNA seq


library(ggplot2)
library(dplyr)
library(tidyr)

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
annotsV$external_gene_name = make.unique(annotsB$external_gene_name)




BAT_EXP = BAT[,!is.na(match(colnames(BAT), annotsB$ensembl_gene_id))]
colnames(BAT_EXP) = annotsB$external_gene_name
BATE = BAT_EXP[(colSums(BAT_EXP[,1:(ncol(BAT_EXP)-3)], na.rm=T) >= 320)]
BATs = cbind(BAT[,(ncol(BAT)-2):ncol(BAT)],BATE)

SAT_EXP = SAT[,!is.na(match(colnames(SAT), annotsS$ensembl_gene_id))]
colnames(SAT_EXP) = annotsS$external_gene_name
SATE = SAT_EXP[(colSums(SAT_EXP[,1:(ncol(SAT_EXP)-3)], na.rm=T) >= 320)]
SATs = cbind(SAT[,(ncol(SAT)-2):ncol(SAT)],SATE)

VAT_EXP = VAT[,!is.na(match(colnames(VAT), annotsV$ensembl_gene_id))]
colnames(VAT_EXP) = annotsV$external_gene_name
VATE = VAT_EXP[(colSums(VAT_EXP[,1:(ncol(VAT_EXP)-3)], na.rm=T) >= 310)]
VATs = cbind(VAT[,(ncol(VAT)-2):ncol(VAT)],VATE)

Liver_EXP = Liver[,!is.na(match(colnames(Liver), annotsL$ensembl_gene_id))]
colnames(Liver_EXP) = annotsL$external_gene_name
LiverE = Liver_EXP[(colSums(Liver_EXP[,1:(ncol(Liver_EXP)-3)], na.rm=T) >= 320)]
Livers = cbind(Liver[,(ncol(Liver)-2):ncol(Liver)],LiverE)

Quad_EXP = Quad[,!is.na(match(colnames(Quad), annotsQ$ensembl_gene_id))]
colnames(Quad_EXP) = annotsQ$external_gene_name
QuadE = Quad_EXP[(colSums(Quad_EXP[,1:(ncol(Quad_EXP)-3)], na.rm=T) >= 320)]
Quads = cbind(Quad[,(ncol(Quad)-2):ncol(Quad)],QuadE)

strain.p <- NA
diet.p <- NA

batch.p <- NA
`strain:diet` <- NA





BAT_1 <- rbind(BATs, strain.p, diet.p,batch.p,`strain:diet`)
rownames(BAT_1)[33:36] <- c("strain.p", "diet.p","batch.p", "strain_diet")
colnames(BAT_1) = colnames(BATs)

SAT_1 <- rbind(SATs, strain.p, diet.p,batch.p,`strain:diet`)
rownames(SAT_1)[33:36] <- c("strain.p", "diet.p","batch.p", "strain_diet")
colnames(SAT_1) = colnames(SATs)

VAT_1 <- rbind(VATs, strain.p, diet.p,batch.p,`strain:diet`)
rownames(VAT_1)[32:35] <- c("strain.p", "diet.p","batch.p", "strain_diet")
colnames(VAT_1) = colnames(VATs)

Liver_1 <- rbind(Livers, strain.p, diet.p,batch.p,`strain:diet`)
rownames(Liver_1)[33:36] <- c("strain.p", "diet.p","batch.p", "strain_diet")
colnames(Liver_1) = colnames(Livers)

Quad_1 <- rbind(Quads, strain.p, diet.p,batch.p,`strain:diet`)
rownames(Quad_1)[33:36] <- c("strain.p", "diet.p","batch.p", "strain_diet")
colnames(Quad_1) = colnames(Quads)




for(i in 5:(ncol(BAT_1))){
  print(i)
  fm <- aov(unlist(BAT_1[1:32,i]) ~ strain * diet +batch, data =BAT_1[1:32,])
  BAT_1[33:36, i] <- summary(fm)[[1]]$`Pr(>F)`[1:4]
}

for(i in 5:(ncol(SAT_1))){
  print(i)
  fm <- aov(unlist(SAT_1[1:32,i]) ~ strain * diet +batch, data =SAT_1[1:32,])
  SAT_1[33:36, i] <- summary(fm)[[1]]$`Pr(>F)`[1:4]
}

for(i in 5:(ncol(VAT_1))){
  print(i)
  fm <- aov(unlist(VAT_1[1:31,i]) ~ strain * diet +batch, data =VAT_1[1:31,])
  VAT_1[32:35, i] <- summary(fm)[[1]]$`Pr(>F)`[1:4]
}

for(i in 5:(ncol(Liver_1))){
  print(i)
  fm <- aov(unlist(Liver_1[1:32,i]) ~ strain * diet +batch, data =Liver_1[1:32,])
  Liver_1[33:36, i] <- summary(fm)[[1]]$`Pr(>F)`[1:4]
}

for(i in 5:(ncol(Quad_1))){
  print(i)
  fm <- aov(unlist(Quad_1[1:32,i]) ~ strain * diet +batch, data =Quad_1[1:32,])
  Quad_1[33:36, i] <- summary(fm)[[1]]$`Pr(>F)`[1:4]
}
#save.image("batb_plantvsan.RData")
#write.table(BAT_1, "1_bat_plant.txt", quote = F, sep = "\t")

BAT_2 <- as.data.frame(t(BAT_1))
VAT_2 <- as.data.frame(t(VAT_1))
SAT_2 <- as.data.frame(t(SAT_1))
Liver_2 <- as.data.frame(t(Liver_1))
Quad_2 <- as.data.frame(t(Quad_1))

BAT_2$strain.p.adj <- NA
BAT_2$diet.p.adj <- NA
BAT_2$strain.diet.p.adj <- NA

VAT_2$strain.p.adj <- NA
VAT_2$diet.p.adj <- NA
VAT_2$strain.diet.p.adj <- NA

SAT_2$strain.p.adj <- NA
SAT_2$diet.p.adj <- NA
SAT_2$strain.diet.p.adj <- NA

Liver_2$strain.p.adj <- NA
Liver_2$diet.p.adj <- NA
Liver_2$strain.diet.p.adj <- NA

Quad_2$strain.p.adj <- NA
Quad_2$diet.p.adj <- NA
Quad_2$strain.diet.p.adj <- NA

BAT_2$strain.p.adj[5:(nrow(BAT_2))] <- p.adjust(as.numeric(as.character(BAT_2$strain.p[5:(nrow(BAT_2))])), method = 'fdr')
BAT_2$diet.p.adj[5:(nrow(BAT_2))] <- p.adjust(as.numeric(as.character(BAT_2$diet.p[5:(nrow(BAT_2))])), method = 'fdr')
BAT_2$strain.diet.p.adj[5:(nrow(BAT_2))] <- p.adjust(as.numeric(as.character(BAT_2$strain_diet[5:(nrow(BAT_2))])), method = 'fdr')

VAT_2$strain.p.adj[5:(nrow(VAT_2))] <- p.adjust(as.numeric(as.character(VAT_2$strain.p[5:(nrow(VAT_2))])), method = 'fdr')
VAT_2$diet.p.adj[5:(nrow(VAT_2))] <- p.adjust(as.numeric(as.character(VAT_2$diet.p[5:(nrow(VAT_2))])), method = 'fdr')
VAT_2$strain.diet.p.adj[5:(nrow(VAT_2))] <- p.adjust(as.numeric(as.character(VAT_2$strain_diet[5:(nrow(VAT_2))])), method = 'fdr')

SAT_2$strain.p.adj[5:(nrow(SAT_2))] <- p.adjust(as.numeric(as.character(SAT_2$strain.p[5:(nrow(SAT_2))])), method = 'fdr')
SAT_2$diet.p.adj[5:(nrow(SAT_2))] <- p.adjust(as.numeric(as.character(SAT_2$diet.p[5:(nrow(SAT_2))])), method = 'fdr')
SAT_2$strain.diet.p.adj[5:(nrow(SAT_2))] <- p.adjust(as.numeric(as.character(SAT_2$strain_diet[5:(nrow(SAT_2))])), method = 'fdr')

Liver_2$strain.p.adj[5:(nrow(Liver_2))] <- p.adjust(as.numeric(as.character(Liver_2$strain.p[5:(nrow(Liver_2))])), method = 'fdr')
Liver_2$diet.p.adj[5:(nrow(Liver_2))] <- p.adjust(as.numeric(as.character(Liver_2$diet.p[5:(nrow(Liver_2))])), method = 'fdr')
Liver_2$strain.diet.p.adj[5:(nrow(Liver_2))] <- p.adjust(as.numeric(as.character(Liver_2$strain_diet[5:(nrow(Liver_2))])), method = 'fdr')

Quad_2$strain.p.adj[5:(nrow(Quad_2))] <- p.adjust(as.numeric(as.character(Quad_2$strain.p[5:(nrow(Quad_2))])), method = 'fdr')
Quad_2$diet.p.adj[5:(nrow(Quad_2))] <- p.adjust(as.numeric(as.character(Quad_2$diet.p[5:(nrow(Quad_2))])), method = 'fdr')
Quad_2$strain.diet.p.adj[5:(nrow(Quad_2))] <- p.adjust(as.numeric(as.character(Quad_2$strain_diet[5:(nrow(Quad_2))])), method = 'fdr')

rownames(BAT_2)[1:(nrow(BAT_2)-4)] = colnames(BATs)
Bpvals = BAT_2[5:(nrow(BAT_2)),37:39]
rownames(VAT_2)[1:(nrow(VAT_2)-4)] = colnames(VATs)
Vpvals = VAT_2[5:(nrow(VAT_2)),36:38]
rownames(SAT_2)[1:(nrow(SAT_2)-4)] = colnames(SATs)
Spvals = SAT_2[5:(nrow(SAT_2)),37:39]
rownames(Liver_2)[1:(nrow(Liver_2)-4)] = colnames(Livers)
Lpvals = Liver_2[5:(nrow(Liver_2)),37:39]
rownames(Quad_2)[1:(nrow(Quad_2)-4)] = colnames(Quads)
Qpvals = Quad_2[5:(nrow(Quad_2)),37:39]

write.table(Bpvals, "pvals_BAT_ANOVA_RNA.txt", quote = F, sep = "\t")
#write.table(BAT_2, "exp_BAT.txt", quote = F, sep = "\t")
write.table(Spvals, "pvals_SAT_ANOVA_RNA.txt", quote = F, sep = "\t")
write.table(Vpvals, "pvals_VAT_ANOVA_RNA.txt", quote = F, sep = "\t")
write.table(Lpvals, "pvals_Liver_ANOVA_RNA.txt", quote = F, sep = "\t")
write.table(Qpvals, "pvals_Quad_ANOVA_RNA.txt", quote = F, sep = "\t")




























