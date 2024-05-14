#DESeq normalization for ANOVA


library(edgeR)
library(stringr)
library(DESeq2)
library(dplyr)
library(tidyr)

#read in csv
counts = read.table("fraction_Nutri.txt", header = T)
row.names(counts) = counts[,1]
counts[,1] = as.factor(counts[,1])
counts[,2:161] = round(counts[,2:161])

annotsV = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = rownames(counts), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsV = as.data.frame(annotsV) %>% filter(external_gene_name != "") %>% na.omit()
annotsV$external_gene_name = make.unique(annotsV$external_gene_name)
counnts = counts[!is.na(match(rownames(counts), annotsV$ensembl_gene_id)),]
rownames(counnts) = annotsV$external_gene_name
counnts[,1] = row.names(counnts) 
counnts[,1] = as.factor(counnts[,1])
counts = counnts

#samples are the rows now
meta = read.table("meta_Nutri_batch.txt", header =T)
rownames(meta) = meta[,1]
meta[,1] = as.factor(meta[,1])
meta[,2] = as.factor(meta[,2])
meta[,3] = as.factor(meta[,3])
meta$Tissue = sub("_.*", "", meta$Sample)
meta$mouseNum = sub(".*_", "", meta$Sample)

meta$Group = factor(paste0(meta$Strain, meta$Diet))

metaVAT = meta %>% filter(Tissue == "VAT")%>% filter(Sample !="VAT_50")
countsVAT = cbind(counts[,1],counts[,130:161])
countsVAT = countsVAT %>% dplyr::select(!VAT_50)
metaVAT$Diet = as.factor(metaVAT$Diet)
metaVAT$Strain = factor(metaVAT$Strain, levels = c("C57","AJ","DBA","SJL"))
metaVAT$Strain = as.factor(metaVAT$Strain)

#DESeq object
DE_obj_VAT <-DESeqDataSetFromMatrix(countData=countsVAT, colData=metaVAT, design=~Strain*Diet+batch, tidy = T)

# run DESeq
DE_res_VAT <-DESeq(DE_obj_VAT, fitType = "local", full = ~Strain*Diet+batch)

resultsNames(DE_res_VAT)
resAM_V <- results(DE_res_VAT,
                   contrast = c('Diet','American','Mediterranean'), alpha = 0.05)
resAVt_V <- results(DE_res_VAT,
                    contrast = c('Diet','American','Vegetarian'), alpha = 0.05)
resAVn_V <- results(DE_res_VAT,
                    contrast = c('Diet','American','Vegan'), alpha = 0.05)
summary(resAM_V)
summary(resAVt_V)
summary(resAVn_V)

pdf("DE_AVt_V.pdf", width = 3, height= 3.5)
EnhancedVolcano(resAVt_V,
                lab = rownames(resAVt_V),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 5e-2,
                FCcutoff = 1, labSize = 3.0)+theme_classic()+  theme(legend.position = "none",text = element_text(size = 10))+ 
  geom_text(data = ann_text, aes( x=15, y=11, label="Vegetarian"), color="black", size=3)+  
  geom_text(data = ann_text, aes( x=-15, y=11, label="American"), color="black", size=3)+
  xlim(-30,30)+
  ylim(0,11)
dev.off()

pdf("DE_AVn_V.pdf", width = 3, height= 3.5)
EnhancedVolcano(resAVn_V,
                lab = rownames(resAVn_V),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 5e-2,
                FCcutoff = 1, labSize = 3.0)+theme_classic()+ theme(legend.position = "none",text = element_text(size = 10))+ 
  geom_text(data = ann_text, aes( x=15, y=35, label="Vegan"), color="black", size=3)+  
  geom_text(data = ann_text, aes( x=-15, y=35, label="American"), color="black", size=3)+
  xlim(-30,30)+
  ylim(0,35)
dev.off()

pdf("DE_AM_V.pdf", width = 3, height= 3.5)
EnhancedVolcano(resAM_V,
                lab = rownames(resAM_V),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 5e-2,
                FCcutoff = 1, labSize = 3.0)+theme_classic()+  theme(legend.position = "none",text = element_text(size = 10))+ 
  geom_text(data = ann_text, aes( x=15, y=16, label="Mediterranean"), color="black", size=3)+  
  geom_text(data = ann_text, aes( x=-15, y=16, label="American"), color="black", size=3)+
  xlim(-30,30)+
  ylim(0,16)
dev.off()







# KEGG pathway enrichment of differentially expressed genes
library(KEGGREST)


pathways.list <- keggList("pathway", "mmu")
head(pathways.list)
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

resAVn_V =na.omit(resAVn_V)
geneList = (resAVn_V$pvalue)
ggg = (resAVn_V$padj)
fold = resAVn_V$log2FoldChange
names(geneList) = (rownames(resAVn_V))


nearestBFD = as.data.frame(names(geneList)[(abs(fold)>=1 & ggg<=0.05)])

colnames(nearestBFD) = 'gene_symbol'
head(geneList)
background = as.data.frame(names(geneList))
colnames(background) = 'gene_symbol'


kegg = as.data.frame(pathway.codes)
kegg$pathscore = 1
colnames(kegg) = c("pathway","pval")




annotsQ = getBM(attributes=c('external_gene_name', 'entrezgene_id'), 
                filters = 'external_gene_name', 
                values = as.character(background$gene_symbol), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsQ = as.data.frame(annotsQ) %>% filter(entrezgene_id != "") %>% na.omit()

library(bc3net)

s=1
for(kd in kegg$pathway){
  pathwaygenes = genes.by.pathway[kd]
  pathwaygenes = pathwaygenes[[1]]
  pathwaygeness = annotsQ$external_gene_name[!is.na(match(as.character(annotsQ$entrezgene_id),as.numeric(pathwaygenes)))]
  
  GWAS <- list(GWASg = nearestBFD$gene_symbol)
  enrich = enrichment(pathwaygeness, as.character(background$gene_symbol), GWAS, adj = "fdr", verbose = FALSE)
  kegg$pval[kegg$pathway == kd] = enrich$padj
  s=s+1
  print(s)
}


enrich = kegg[kegg$pval<0.05,]
paths = as.data.frame(pathways.list)
paths$pathway = rownames(paths)

enriched = inner_join(enrich, paths, by = 'pathway')
enriched$namess = gsub("- Mus musculus \\(house mouse\\)", "", enriched$pathways.list)


enriched %>% ggplot(aes(x = reorder(namess, -log10(pval)), y = -log10(pval)))+
  geom_bar(stat = "identity")+coord_flip()+xlab("Adjusted P-value")+ylab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 10))
enriched = enriched[enriched$pval<0.005,]

pdf("kegg_nutri_VATveganDEres.pdf",height = 3.5, width = 3)

enriched %>% ggplot(aes(x = reorder(namess, -log10(pval)), y = -log10(pval)))+
  geom_bar(stat = "identity")+coord_flip()+ylab("Adjusted P-value")+xlab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 8))
dev.off()


