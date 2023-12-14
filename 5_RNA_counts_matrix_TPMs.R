# process RNA seq counts

setwd("/scratch/jnr3hh/htseq2")


### unite the 160 counts vectors into one matrix
#install.packages("BiocManager")
#BiocManager::install("edgeR")
#install.packages("stringr")
library(edgeR)
library(stringr)
library(DESeq2)

files <- list.files("/scratch/jnr3hh/htseq", pattern = ".*geneCounts.bam")
samples = str_remove(files, "_out_geneCounts.bam")

dge_counts <- readDGE(files, labels=samples, header=FALSE)

write.csv(dge_counts$counts, 'Nutrigenetics_counts_fraction.csv')


#TPMs for PCA

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

l = listAttributes(mouse)

annots = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
               filters = 'ensembl_gene_id', 
               values = rownames(counts), 
               mart = mouse)

library(GenomicFeatures)
mouseTxDB = makeTxDbFromGFF("Mus_musculus.GRCm39.108.gtf",format=("gtf"))
exons.list.per.gene <- exonsBy(mouseTxDB,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
exonic.gene.sizes$ensembl_gene_id = rownames(exonic.gene.sizes)

library(dplyr)
sizes = inner_join(exonic.gene.sizes, annots, by = "ensembl_gene_id")
colnames(sizes) = c("gene_length","ensembl_gene_id", "external_gene_name")

counts_RPK = counts[,2:161]
for (i in 1:nrow(counts)){
  counts_RPK[i,]= counts[i,2:161]/sizes$gene_length[i]
  print(i)
}
counts_RPK_scaled = counts_RPK*1000
scaling_factor = colSums(counts_RPK_scaled)/ 1000000

TPMs = counts_RPK_scaled
for (i in 1:ncol(counts_RPK_scaled)){
  TPMs[,i]= counts_RPK_scaled[,i]/scaling_factor[i]
  print(i)
}

write.table(TPMs, "Nutrigenetics_TPMs.txt", quote = F, row.names = T, sep = "\t")


#log scale
TPMs_log = log(TPMs+1)
TPMs_log = TPMs_log[rowSums(TPMs_log != 0) > 0,]


tpms = as.data.frame(t(TPMs_log))
tpms$Strain = meta$Strain
tpms$Diet = meta$Diet
meta$Tissue = sub("_.*", "", meta$Sample)
tpms$Tissue = meta$Tissue


#################################
#DESeq normalization for ANOVA
#read in csv
counts = read.table("fraction_Nutri.txt", header = T)
row.names(counts) = counts[,1]
counts[,1] = as.factor(counts[,1])
counts[,2:161] = round(counts[,2:161])


#samples are the rows now
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
metaVAT = meta %>% filter(Tissue == "VAT")%>% filter(Sample !="VAT_50")
metaBAT = meta %>% filter(Tissue == "BAT")

countsLiver = cbind(counts[,1],counts[,34:65])
countsQuad = cbind(counts[,1],counts[,66:97])
countsBAT = cbind(counts[,1],counts[,2:33])
countsSAT = cbind(counts[,1],counts[,98:129])
countsVAT = cbind(counts[,1],counts[,130:161])
countsVAT = countsVAT %>% select(!VAT_50)




library(edgeR)
library(stringr)
library(DESeq2)
library(dplyr)
library(tidyr)



DE_obj_Liver <-DESeqDataSetFromMatrix(countData=countsLiver, colData=metaLiver, design=~Strain*Diet+batch, tidy = T)
DE_obj_Quad <-DESeqDataSetFromMatrix(countData=countsQuad, colData=metaQuad, design=~Strain*Diet+batch, tidy = T)
DE_obj_SAT <-DESeqDataSetFromMatrix(countData=countsSAT, colData=metaSAT, design=~Strain*Diet+batch, tidy = T)
DE_obj_VAT <-DESeqDataSetFromMatrix(countData=countsVAT, colData=metaVAT, design=~Strain*Diet+batch, tidy = T)
DE_obj_BAT <-DESeqDataSetFromMatrix(countData=countsBAT, colData=metaBAT, design=~Strain*Diet+batch, tidy = T)


#######

DE_res_Liver <-DESeq(DE_obj_Liver, fitType = "local", full = ~Strain*Diet+batch)
DE_res_Quad <-DESeq(DE_obj_Quad, fitType = "local", full = ~Strain*Diet+batch)
DE_res_BAT <-DESeq(DE_obj_BAT, fitType = "local", full = ~Strain*Diet+batch)
DE_res_SAT <-DESeq(DE_obj_SAT, fitType = "local", full = ~Strain*Diet+batch)
DE_res_VAT <-DESeq(DE_obj_VAT, fitType = "local", full = ~Strain*Diet+batch)



liver_norm_counts = as.data.frame(counts(DE_res_Liver, normalized=TRUE))
Quad_norm_counts = as.data.frame(counts(DE_res_Quad, normalized=TRUE))
SAT_norm_counts = as.data.frame(counts(DE_res_SAT, normalized=TRUE))
VAT_norm_counts = as.data.frame(counts(DE_res_VAT, normalized=TRUE))
BAT_norm_counts = as.data.frame(counts(DE_res_BAT, normalized=TRUE))

BAT_coll_zero = replace(BAT_norm_counts, BAT_norm_counts<=0.1, 0)
VAT_coll_zero = replace(VAT_norm_counts, VAT_norm_counts<=0.1, 0)
SAT_coll_zero = replace(SAT_norm_counts, SAT_norm_counts<=0.1, 0)
Liver_coll_zero = replace(liver_norm_counts, liver_norm_counts<=0.1, 0)
Quad_coll_zero = replace(Quad_norm_counts, Quad_norm_counts<=0.1, 0)

Bsums = rep(0,nrow(BAT_coll_zero))
Vsums = rep(0,nrow(VAT_coll_zero))
Ssums = rep(0,nrow(SAT_coll_zero))
Lsums = rep(0,nrow(Liver_coll_zero))
Qsums = rep(0,nrow(Quad_coll_zero))

#calculate number of zeros
for(k in 1:nrow(BAT_coll_zero)){
  print(k)
  Bsums[k] = sum(BAT_coll_zero[k,]==0)
  Vsums[k] = sum(VAT_coll_zero[k,]==0)
  Ssums[k] = sum(SAT_coll_zero[k,]==0)
  Lsums[k] = sum(Liver_coll_zero[k,]==0)
  Qsums[k] = sum(Quad_coll_zero[k,]==0)
}

BAT = BAT_coll_zero[(Bsums<= 32*0.75),]
VAT = VAT_coll_zero[(Vsums<= 31*0.75),]
SAT = SAT_coll_zero[(Ssums<= 32*0.75),]
Liver = Liver_coll_zero[(Lsums<= 32*0.75),]
Quad = Quad_coll_zero[(Qsums<= 32*0.75),]



write.table(Liver,"liver_deseq_norm_counts.txt", quote = F, sep = "\t")
write.table(Quad,"Quad_deseq_norm_counts.txt", quote = F, sep = "\t")
write.table(BAT,"BAT_deseq_norm_counts.txt", quote = F, sep = "\t")
write.table(SAT,"SAT_deseq_norm_counts.txt", quote = F, sep = "\t")
write.table(VAT,"VAT_deseq_norm_counts.txt", quote = F, sep = "\t")






