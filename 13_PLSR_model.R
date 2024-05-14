
setwd("~/")

#BiocManager::install("mixOmics", force = T)
library(mixOmics)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

#load phenotype data
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")
#load("nutrigenetics_phenotypeData.RData")

#setup RNA data 
#samples are the rows now
meta = read.table("meta_Nutri_batch.txt", header =T)
rownames(meta) = meta[,1]
meta[,1] = as.factor(meta[,1])
meta[,2] = as.factor(meta[,2])
meta[,3] = as.factor(meta[,3])
meta$Tissue = sub("_.*", "", meta$Sample)
meta$mouseNum = sub(".*_", "", meta$Sample)
meta = meta %>% filter(Sample != "VAT_50")
metaVAT = meta %>% filter(Tissue == "VAT") %>% filter(Sample != "VAT_50") #outlier

TPMs = read.table("Nutrigenetics_TPMs.txt", header = T)
TPMs = TPMs[,-151]
#log scale
TPMs_log = log(TPMs+1)
TPMs_log = TPMs_log[rowSums(TPMs_log != 0) > 0,]

tpms = as.data.frame(t(TPMs_log))
tpms$Strain = meta$Strain
tpms$Diet = meta$Diet
tpms$Sample = meta$Sample
tpms$Tissue = meta$Tissue
tpms$Diet = factor(tpms$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

tpms = tpms %>% filter(Sample != "VAT_50")

tpms = na.omit(tpms)
tpmsVAT = tpms %>% filter(Tissue == "VAT")  %>% filter(Sample != "VAT_50")
tpmsVAT = tpmsVAT[,colSums(tpmsVAT[,1:44809]) != 0]

VAT <- as.data.frame(tpmsVAT)
VAT$strain = metaVAT$Strain
VAT$diet = metaVAT$Diet
VAT$strain = factor(VAT$strain)
VAT$diet = factor(VAT$diet)
VAT$batch = metaVAT$batch
VAT$batch = factor(VAT$batch)

annotsV = getBM(attributes=c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'), 
                filters = 'ensembl_gene_id', 
                values = colnames(VAT), 
                mart = mart) %>%
  distinct() %>%
  as_tibble()
annotsV = as.data.frame(annotsV) %>% filter(external_gene_name != "") %>% na.omit()
annotsV$external_gene_name = make.unique(annotsV$external_gene_name)

VAT_EXP = VAT[,!is.na(match(colnames(VAT), annotsV$ensembl_gene_id))]
colnames(VAT_EXP) = annotsV$external_gene_name


VAT_aov_p = read.table("pvals_VAT_B_75.txt", header= T)
VAT_aov_pv = VAT_aov_p[!is.na(match(rownames(VAT_aov_p), annotsV$ensembl_gene_id)),]
rownames(VAT_aov_pv) = annotsV$external_gene_name
VAT_sigGenes = VAT_aov_pv %>% filter(strain.diet.p.adj <= 0.05)
#VAT_sigGenes = VAT_aov_pv %>% filter(strain_diet <= 0.1)
VAT_expr = VAT_EXP %>% dplyr::select(c(rownames(VAT_sigGenes)))


# Remove Duplicate Column Names
VAT_expr = VAT_expr[!duplicated(colnames(VAT_expr))]

#select samples
exp2_half = exp2[!is.na(match(exp2$MouseNum, gsub("VAT_","",rownames(VAT_expr)))),]
rownames(exp2_half) = exp2_half$MouseNum
VAT_expr$MouseNum = as.numeric(gsub("VAT_","",rownames(VAT_expr)))
rownames(VAT_expr) = VAT_expr$MouseNum
VAT_expr = VAT_expr %>% arrange(MouseNum)
VAT_expr = VAT_expr[,1:417]
#VAT_expr = VAT_expr[,1:4437]


# load diet info
dietcomp = read.table("dietcomp.txt", header = T, sep = "\t", fill = T) %>% filter(MouseNum != 19) %>% mutate(MouseNum = as.character(MouseNum)) %>% arrange(as.numeric(MouseNum)) #supplemental table 1
rownames(dietcomp) = paste0("m", dietcomp$MouseNum)
dietsource = read.table("dietsource.txt", header = T, sep = "\t", fill = T) %>% filter(MouseNum != 19) %>% mutate(MouseNum = as.character(MouseNum)) %>% arrange(as.numeric(MouseNum)) #supplemental table 2
rownames(dietsource) = paste0("m",dietsource$MouseNum)
X1 = as.data.frame((dietsource[4:ncol(dietsource)]))
X2 = as.data.frame((dietcomp[4:ncol(dietcomp)]))



# make sure the row names are the same
Y1 = VAT_expr
Y2 = exp2_half[,4:28]
rownames(Y1)= paste0("m",rownames(Y1))
rownames(Y2)= paste0("m",rownames(Y2))
Y1 = Y1[!is.na(match(rownames(Y1), rownames(X1))),]
Y2 = Y2[!is.na(match(rownames(Y2), rownames(X1))),]
X1 = X1[!is.na(match(rownames(X1), rownames(Y1))),]
X2 = X2[!is.na(match(rownames(X2), rownames(Y1))),]

X1 = data.matrix(X1)
Y1 = data.matrix(Y1)

X2 = data.matrix(X2)
Y2 = data.matrix(Y2)


#Y2 = as.data.frame(Y) %>% dplyr::select(!Week_5_MRI) %>% na.omit()
#X1 = X[!is.na(match(rownames(X), rownames(Y2))),]

pls.nutri <- mixOmics::pls(Y1, Y2, ncomp = 3) # undergo PLS regression

#model validation
#valdcr <- perf(pls.nutri, validation = 'Mfold',
#               folds = 10, nrepeat = 10) 
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


# pull out component space info to make biplot
#samples (mice)
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

# diet component and RNA
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
Yload = Yload %>% filter(!grepl('Gm|Rik', module)) %>% filter(comp1 <= -.55|comp1 >=.55|comp2>=.51|comp2<=-.55)
Xload = Xload %>% filter(!grepl('Gm|Rik', module)) %>% filter(comp1 <= -.4|comp1 >=.4|comp2>=.4|comp2<=-.4) #filter by removing pseudogenes and lnc (optional), filter unimportant genes by component 1 and 2
Xload$space = "X"
Yload$space = "Y"

YX = rbind(Xload, Yload)

colo1 = c("grey70","mediumpurple4","#fa76ef","red3","yellow","orange","forestgreen","blue","grey20","darkblue")

#colo2 = c(colo2[-5], "black","blue")
normXYspace$sma = paste0(normXYspace$Strain, normXYspace$Diet)
normXYspace$Diet = factor(normXYspace$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan","X","Y"))
ann_text<-data.frame(
  y = c(65), x = c(7),
  label = "\U25D7")


library(extrafont) #may need to install new fonts first, half circles wont plot in basic arial or many basic fonts
loadfonts(device = "win", quiet = TRUE)

bipl = normXYspace %>% ggplot(aes(x = comp1, y = comp2))+
  scale_color_manual(values = colo1 )+
  geom_text(label = '\u25D6', family = "Lucida Sans Unicode", aes(color = Strain),  
            size=5, hjust = .75)+
  geom_text(label = '\u25D7', family = "Lucida Sans Unicode", aes(color = Diet),  
            size=5, hjust = .25)+
  geom_text_repel(data = YX, family = "Arial", size = 2.25, aes(x = comp1, y = comp2,label = module, color = space), max.overlaps = 1000)+
  
  theme_bw()+
  xlab("Component 1")+ ylab("Component 2")+
  theme(text = element_text(size = 7), legend.position = "none")

bipl
cairo_pdf("my_plot_dietran.pdf", 6,6)
bipl
dev.off()

