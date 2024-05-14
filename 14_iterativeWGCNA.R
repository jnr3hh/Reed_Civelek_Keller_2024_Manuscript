##################################values = 
## WGCNA module investigation- nutrigenetics

setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")

library(dplyr)
library(ggplot2)
library(tidyr)
library(WGCNA)

#####
# 1. Preprocessing

#1.1 remove zeros
TPMs = read.table("Nutrigenetics_TPMs.txt", header = T)
Nutri_PC_coll_zero = replace(TPMs, TPMs<=0.1, 0)

###repeat for other tissues
VAT_PC_coll_zero = Nutri_PC_coll_zero[,129:160]
Vsums = rep(0,nrow(VAT_PC_coll_zero))

#calculate number of zeros
for(k in 1:nrow(VAT_PC_coll_zero)){
  Vsums[k] = sum(VAT_PC_coll_zero[k,]==0)
  print(k)
}

#remove rows with more than 80% zero
VAT_PC_coll_nozero = VAT_PC_coll_zero[(Vsums< 32*0.8),]
#log 
VAT_gene_exp_processed = log2(VAT_PC_coll_nozero+1)
VAT_gene_exp_processed = VAT_gene_exp_processed %>% dplyr::select(!VAT_50)
write.table(VAT_gene_exp_processed, "VAT_processed_geneExp.txt", quote = F, sep = "\t")



#1.2 Choose iterativeWGCNA power

library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)

library(WGCNA)
library(sva)
library(colorspace)
library(qvalue)
library(cluster)
library(flashClust)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
allowWCGNAThreads()


input = t(Quad_gene_exp_processed) #repeat for others

# check for genes and samples with too many missing values in the null data
# if gsg$allOK returns true, all genes have passed the cuts. If not,
# we remove the offending genes and samples from the data set
gsg <- goodSamplesGenes(input, verbose=3)
dim0 = dim(input)
if(!gsg$allOK) {  
  if(sum(!gsg$goodGenes) > 0)
    printFlush(paste('Removing genes:', paste(names(input)[!gsg$goodGenes], collapse = ', ')))
  if(sum(!gsg$goodSamples) > 0)
    printFlush(paste('Removing samples:', paste(rownames(input)[!gsg$goodSamples], collapse = ', ')))
  # Remove offending genes and samples from the data
  input <- input[gsg$goodSamples, gsg$goodGenes]
}
if(dim0[1]!=nrow(input)){paste0("Null samples were removed")}

# Cluster the Null samples to ID outliers (Euclidian distance)
sampleTree = flashClust(dist(input), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# apply a cut after reviewing this plot (where?)
cut = 170 #maybe change this??????
# Plot a line to show the cut
pdf(file= "prePro_Null_sample_clustering_to_detect_outliers.pdf", width=12, height=9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
input = input[keepSamples, ]
length(which(keepSamples==F))

# range of powers for consideration
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# select power(s)
sft = pickSoftThreshold(input, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#once you get to here, look at the two graphs
# choose the power where both lines level off
# ideally the degree is above .8
# when in doubt, choose the lower number for power



###############
# 2. Run iterative
#caution, these commands must be run in python


#get to your python3 command line or anaconda window
# cd C:\Users\jnr3hh\Documents\iterative\Nutri_indiv

# load the module with:
# python3 -m pip install iterativeWGCNA

# run iterative on one sample

# python3 -m iterativeWGCNA -i /oldscratch/jnr3hh/iterative/VAT/VAT_processed_geneExp.txt --enableWGCNAThreads --wgcnaParameters maxBlockSize=50000,power=10

#swap out the power for what you found above
#make sure maxBlockSize is bigger than the number of genes you have

#if you start with 20000 or more genes, it will take 12-36 hours to run on a normal machine (can move to rivanna)

# keep the output files "merged-0.05-eigengenes.txt" and "merged-0.05-membership.txt"


################ Analysis of module genes

# read in module and eigengene files from iterative WGCNA
modules = read.table("merged-0.05-membership_VAT.txt", header = T)
eigengenes = read.table("merged-0.05-eigengenes_VAT.txt", header = T)
moduleName = eigengenes[,1]

#correlation of module eigengene to phenotypes
eigs = as.data.frame(t(eigengenes[2:32]))
colnames(eigs) = moduleName
eigs$Tissue = sub("_.*", "", rownames(eigs))
eigs$mouseNum = sub(".*_", "", rownames(eigs))

sumeigs = as.data.frame(eigs %>% group_by(mouseNum) %>% summarise(across(1:131, ~mean(.x))))
sumeigs = sumeigs %>% arrange(as.numeric(mouseNum))
all_sum = exp2[!is.na(match(exp2$MouseNum, sumeigs$mouseNum)),]

sum_eig_cor = matrix(nrow = ncol(sumeigs[,2:132]),ncol = ncol(all_sum[,4:28]))
rownames(sum_eig_cor) = colnames(sumeigs[,2:132])
colnames(sum_eig_cor) = colnames(all_sum[,4:28])
sum_eig_pval = matrix(nrow = ncol(sumeigs[,2:132]),ncol = ncol(all_sum[,4:28]))
rownames(sum_eig_pval) = colnames(sumeigs[,2:132])
colnames(sum_eig_pval) = colnames(all_sum[,4:28])


for (j in 1:ncol(sum_eig_cor)){
  for(i in 1:nrow(sum_eig_cor)){
    print(j)
    if(length(all_sum[,j+3])==length(sumeigs[,i+1])){
      l = cor.test(sumeigs[,i+1], all_sum[,j+3])
      sum_eig_cor[i,j] = l$estimate 
      sum_eig_pval[i,j] = l$p.value 
    }
  }
}
sum_eig_padj = matrix(p.adjust(as.vector(sum_eig_pval), method = "fdr", n = length(sum_eig_pval)), ncol = ncol(sum_eig_pval))
rownames(sum_eig_padj) = rownames(sum_eig_pval)
colnames(sum_eig_padj) = colnames(sum_eig_pval)
sum_eig_pstar = as.data.frame(matrix(NA_integer_, nrow = nrow(sum_eig_padj), ncol = ncol(sum_eig_padj)))
colnames(sum_eig_pstar) = colnames(sum_eig_padj)
rownames(sum_eig_pstar) = rownames(sum_eig_padj)

for (j in 1:ncol(sum_eig_cor)){
  for(i in 1:nrow(sum_eig_cor)){
    print(j)
    if(sum_eig_padj[i,j] <= 0.05){
      sum_eig_pstar[i,j] = "*"
    }
    if(sum_eig_padj[i,j] < 0.01){
      sum_eig_pstar[i,j] = "**"
    }
    if(sum_eig_padj[i,j] < 0.001){
      sum_eig_pstar[i,j] = "***"
    }
    if(sum_eig_padj[i,j] > 0.05) {
      sum_eig_pstar[i,j] = ""}
  }
}

intMods = moduleName
sum_eig_pstars = as.data.frame(sum_eig_pstar)

sum_eig_pstar = as.matrix(sum_eig_pstars)

library(pheatmap)
# plot correlations
pdf("heatmaptraitcorr_WGCNAmod.pdf", height = 4, width = 3)
sumHeat = pheatmap(t(sum_eig_cor[!is.na(match(rownames(sum_eig_cor),intMods)),]), display_numbers = t(sum_eig_pstars[!is.na(match(rownames(sum_eig_pstars),intMods)),]), fontsize = 10, treeheight_row = 5, treeheight_col = 5)
dev.off()


sum_eigsCor = as.data.frame(sum_eig_cor) %>% mutate(module = rownames(sum_eig_cor)) %>%  pivot_longer(cols=1:ncol(sum_eig_cor),names_to = "Trait", values_to = "Correlation")
sum_eigspadj = as.data.frame(sum_eig_padj) %>% mutate(module = rownames(sum_eig_cor)) %>% pivot_longer(cols=1:ncol(sum_eig_cor),names_to = "Trait", values_to = "p.adj")
sum_eigsCor$padj = sum_eigspadj$p.adj
top = sum_eigsCor %>% filter(padj < 0.05)
r = unique(top$module)




# perform 2-way ANOVA on module eigengenes
eigsann = cbind(eigs[1:131], metaVAT$Strain, metaVAT$Diet)
colnames(eigsann)[132:133] = c("Strain","Diet")


exp1 = eigsann
exp1$Strain = factor(exp1$Strain)
exp1$Diet = factor(exp1$Diet)

strain.p <- NA
diet.p <- NA
`strain:diet` <- NA

exp1_AOV <- rbind(exp1, strain.p, diet.p,`strain:diet`)
rownames(exp1_AOV)[(nrow(exp1)+1):(nrow(exp1)+3)] <- c("strain.p", "diet.p", "strain_diet")

for(i in 1:(ncol(exp1_AOV)-2)){
  print(i)
  fm <- aov(unlist(exp1_AOV[1:nrow(exp1),i]) ~ Strain * Diet, data =exp1_AOV[1:nrow(exp1),])
  exp1_AOV[(nrow(exp1)+1):(nrow(exp1)+3), i] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}



exp1_AOVres <- as.data.frame(t(exp1_AOV))
exp1_AOVres$strain.p  = as.numeric(exp1_AOVres$strain.p)
exp1_AOVres$diet.p  = as.numeric(exp1_AOVres$diet.p)
exp1_AOVres$strain_diet  = as.numeric(exp1_AOVres$strain_diet)



exp1_AOVres$strain.p.adj <- NA
exp1_AOVres$diet.p.adj <- NA
exp1_AOVres$strain.diet.p.adj <- NA

exp1_AOVres$strain.p.adj[1:(nrow(exp1_AOVres)-2)] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain.p[1:(nrow(exp1_AOVres)-2)])), method = 'fdr')
exp1_AOVres$diet.p.adj[1:(nrow(exp1_AOVres)-2)] <- p.adjust(as.numeric(as.character(exp1_AOVres$diet.p[1:(nrow(exp1_AOVres)-2)])), method = 'fdr')
exp1_AOVres$strain.diet.p.adj[1:(nrow(exp1_AOVres)-2)] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain_diet[1:(nrow(exp1_AOVres)-2)])), method = 'fdr')



sigmods = na.omit(exp1_AOVres[exp1_AOVres$strain.diet.p.adj <= 0.05,])
intMods = rownames(sigmods)



# module enrichment for KEGG pathways
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
modules$ensembl_gene_id = modules$Gene

annotsV = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = modules$Gene, 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsV = as.data.frame(annotsV) %>% filter(external_gene_name != "") %>% na.omit()
annotsV$external_gene_name = make.unique(annotsV$external_gene_name)
mods = na.omit(inner_join(modules, annotsV, by = "ensembl_gene_id"))


geneList <- mods$kME
names(geneList) <- mods$external_gene_name
head(geneList)
background = as.data.frame(names(geneList))
colnames(background) = 'gene_symbol'


kegg = as.data.frame(pathway.codes)
kegg$pathscore = 1
colnames(kegg) = c("pathway","pval")


nearestBFD = as.data.frame(mods$external_gene_name[mods$Module == intMods[3]])
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
  geom_bar(stat = "identity")+coord_flip()+ylab("Adjusted P-value")+xlab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 10))+ggtitle("P1_I9_M5")

pdf("kegg_nutri_421.pdf",height = 1, width = 4)

enriched %>% ggplot(aes(x = reorder(namess, -log10(pval)), y = -log10(pval)))+
  geom_bar(stat = "identity")+coord_flip()+ylab("Adjusted P-value")+xlab("KEGG pathway")+theme_classic()+theme(text = element_text(size = 10))+ggtitle("P2_I12_M12")
dev.off()


# plot individual genes within module 4
eiggenes = as.data.frame(eigs[,colnames(eigs)==intMods[4]])
rownames(eiggenes) = rownames(eigs)
eiggenes$Strain = metaVAT$Strain
eiggenes$Diet = metaVAT$Diet
colnames(eiggenes)[1] = "mod"

genesinmod = mods$external_gene_name[mods$Module == intMods[4]]
exp = VAT_EXP[,!is.na(match(colnames(VAT_EXP),c("Slc22a15","Tmem144","Mkx")))]
exp <- mutate_all(exp, function(x) as.numeric(as.character(x)))

exp$Strain = metaVAT$Strain
exp$Diet = metaVAT$Diet
exp3 = exp %>% pivot_longer(cols=1:3, names_to = "genes", values_to = "Exp")

exp4 = exp3 %>% group_by(Strain, Diet, genes) %>% summarise(aExp = mean(Exp))

exp3$Diet = factor(exp3$Diet, levels = c("American","Mediterranean", "Vegetarian","Vegan"))

pdf("mkx.pdf", height = 2, width = 2)
exp3 %>% filter(genes == "Mkx") %>%  group_by(Strain, Diet)  %>% summarise(meanPC = mean((as.numeric(Exp))),sePC = sd((as.numeric(Exp)))/sqrt(length((as.numeric(Exp))))) %>% ggplot(aes(x= Strain, y = meanPC, group = Strain))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),width = 0.5, aes(color = Strain, fill = Strain))+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), position = position_dodge(width = 0.8),width = .3, color = "black")+
  theme_classic()+
  theme(text = element_text(size = 8),axis.text.x = element_text(angle = 45,  hjust=1), legend.position = "none")+
  xlab("")+
  ylab("Gene Expression")+
  scale_color_manual(values = colo1)+
  scale_fill_manual(values = colo1)+ggtitle("Mkx")+
  facet_wrap(~Diet, nrow =2)
dev.off()

x













## eigengene expression per module plots

eiggenes$Diet = factor(eiggenes$Diet, levels = c("American","Mediterranean", "Vegetarian","Vegan"))

VAT_expr %>%  select(Strain, Diet, Slc22a15)


pdf("p3_eig_expression.pdf", height = 2, width = 2)
eiggenes %>% ggplot(aes(x = Strain, y = mod, color = Strain, fill = Strain))+
  geom_bar(stat = "identity")+facet_wrap(~Diet)+ylab("Relative Eigengene Expression")+theme_classic()+theme(legend.position = "none",text = element_text(size = 8))+scale_fill_manual(values = colo1)+scale_color_manual(values = colo1)+ggtitle("Module 3")
dev.off()

fda = read.table("novelwgnca.txt", header = T, sep= "\t")
pdf("p5_eig_expression.pdf", height = 2, width = 1.5)
fda %>% pivot_longer(2:3, names_to = "Type", values_to= "Number") %>% ggplot(aes(x = Module, y = Number, fill = Type))+
  geom_bar(stat = "identity",position = "stack")+theme_classic()+scale_fill_manual(values = c("grey40","grey80"))+theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), text = element_text(size = 10), axis.title.x = element_text(color = "white"))
dev.off()



















