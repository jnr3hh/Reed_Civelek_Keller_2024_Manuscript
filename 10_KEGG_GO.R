#### KEGG + GO pathway enrichment

library(ggplot2)
library(dplyr)
library(tidyr)
library(KEGGREST)
library(org.Mm.eg.db)
library(bc3net)

#KEGG #adapted from https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
#load kegg pathways
pathways.list <- keggList("pathway", "mmu")
head(pathways.list)
pathway.codes <- sub("path:", "", names(pathways.list)) 
#get genes in each pathway
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

#load a list of genes in visceral adipose
DE.table1 = read.table("VAT_AOV_table.txt", header = T, sep = "\t")

geneList <- DE.table1$strain.diet.p.adj
names(geneList) <- DE.table1$external_gene_name
head(geneList)
background = as.data.frame(names(geneList))
colnames(background) = 'gene_symbol'

geneList <- DE.table1$Gene
names(geneList) <- DE.table1$external_gene_name
head(geneList)
background = as.data.frame(names(geneList))
colnames(background) = 'gene_symbol'



kegg = as.data.frame(pathway.codes)
kegg$pathscore = 1
colnames(kegg) = c("pathway","pval")

VAT_aov_p = read.table("pvals_VAT_B_75.txt", header= T)
VAT_aov_pv = VAT_aov_p[!is.na(match(rownames(VAT_aov_p), annotsV$ensembl_gene_id)),]
rownames(VAT_aov_pv) = annotsV$external_gene_name
VAT_sigGenes = VAT_aov_pv %>% filter(strain.diet.p.adj <= 0.05)

nearestBFD = as.data.frame(VAT_sigGenes)
colnames(nearestBFD) = 'gene_symbol'

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
  geom_bar(stat = "identity")+coord_flip()+ylab("Adjusted P-value")+xlab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 10))

pdf("kegg_nutri_p5.pdf",height = 1, width = 4)

enriched %>% ggplot(aes(x = reorder(namess, -log10(pval)), y = -log10(pval)))+
  geom_bar(stat = "identity")+coord_flip()+ylab("Adjusted P-value")+xlab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 10))
dev.off()

x





########### GO enrichment
VAT_aov_p = read.table("pvals_VAT_B_75.txt", header= T)
VAT_aov_pv = VAT_aov_p[!is.na(match(rownames(VAT_aov_p), annotsV$ensembl_gene_id)),]
rownames(VAT_aov_pv) = annotsV$external_gene_name
VAT_sigGenes = VAT_aov_pv %>% filter(strain.diet.p.adj <= 0.05)



library(topGO)
library(ALL)
library(org.Mm.eg.db)

DSGO_module <- data.frame(GO.ID=character(),
                             Term=character(), 
                             Annotated=numeric(),
                             Significant=numeric(),
                             Expected=numeric(),
                             classicFisher=numeric(),
                             classicFisherA = numeric(),
                             mod=character(),
                             stringsAsFactors=FALSE) 



geneNames = colnames(VAT_EXP)

#downstreamGenes = Moccur3

downstreamGenes = VAT_sigGenes
if(length(downstreamGenes != 0)){
  myInterestingGenes <- (t(downstreamGenes)) # int genes ( all the genes in each modules)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames # gene names and 1 or 0 in their cell s
  
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultFisher
  allGO = usedGO(object = GOdata) 
  allRes <- GenTable(GOdata, classicFisher = resultFisher,ranksOf = "classicFisher", topNodes = length(allGO))
  allRes$classicFisher = sub("< ","", allRes$classicFisher)
  a <- p.adjust(allRes$classicFisher, method = "fdr", n = length(allRes$classicFisher))
  allRes$classicFisherA <- a 
  # GO term result file 
  
  sig = allRes[allRes$classicFisherA <= 0.05,]
  if(nrow(sig) != 0){
    
    sig$mod = mods
    DSGO_module = rbind(DSGO_module, sig)
  }
}    


























