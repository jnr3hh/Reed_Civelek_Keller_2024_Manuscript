
setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)

#load data
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")
load("nutrigenetics_phenotypeData.RData")
exp2 = exp2 %>% filter(MouseNum != 19)

exp1$Strain = factor(exp1$Strain)
exp1$Diet = factor(exp1$Diet)

strain.p <- NA
diet.p <- NA
`strain:diet` <- NA

exp1_AOV <- rbind(exp1, strain.p, diet.p,`strain:diet`)
rownames(exp1_AOV)[(nrow(exp1)+1):(nrow(exp1)+3)] <- c("strain.p", "diet.p", "strain_diet")

for(i in 1:(ncol(exp1_AOV)-3)){
  print(i)
  fm <- aov(unlist(exp1_AOV[1:nrow(exp1),i+3]) ~ Strain * Diet, data =exp1_AOV[1:nrow(exp1),])
  exp1_AOV[(nrow(exp1)+1):(nrow(exp1)+3), i+3] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}



exp1_AOVres <- as.data.frame(t(exp1_AOV))
exp1_AOVres$strain.p  = as.numeric(exp1_AOVres$strain.p)
exp1_AOVres$diet.p  = as.numeric(exp1_AOVres$diet.p)
exp1_AOVres$strain_diet  = as.numeric(exp1_AOVres$strain_diet)



exp1_AOVres$strain.p.adj <- NA
exp1_AOVres$diet.p.adj <- NA
exp1_AOVres$strain.diet.p.adj <- NA

exp1_AOVres$strain.p.adj[4:(nrow(exp1_AOVres))] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain.p[4:(nrow(exp1_AOVres))])), method = 'fdr')
exp1_AOVres$diet.p.adj[4:(nrow(exp1_AOVres))] <- p.adjust(as.numeric(as.character(exp1_AOVres$diet.p[4:(nrow(exp1_AOVres))])), method = 'fdr')
exp1_AOVres$strain.diet.p.adj[4:(nrow(exp1_AOVres))] <- p.adjust(as.numeric(as.character(exp1_AOVres$strain_diet[4:(nrow(exp1_AOVres))])), method = 'fdr')

write.table(exp1_AOVres[,(ncol(exp1_AOVres)-5):ncol(exp1_AOVres)],"Exp1_ANOVA_trait_pvals.txt", quote = F, sep = "\t")


exp2$Strain = factor(exp2$Strain)
exp2$Diet = factor(exp2$Diet)

strain.p <- NA
diet.p <- NA
`strain:diet` <- NA


exp2_AOV <- rbind(exp2, strain.p, diet.p,`strain:diet`)
rownames(exp2_AOV)[(nrow(exp2)+1):(nrow(exp2)+3)] <- c("strain.p", "diet.p", "strain_diet")

for(i in 1:(ncol(exp2_AOV)-3)){
  print(i)
  fm <- aov(unlist(exp2_AOV[1:nrow(exp2),i+3]) ~ Strain * Diet, data =exp2_AOV[1:nrow(exp2),])
  exp2_AOV[(nrow(exp2)+1):(nrow(exp2)+3), i+3] <- summary(fm)[[1]]$`Pr(>F)`[1:3]
}



exp2_AOVres <- as.data.frame(t(exp2_AOV))
exp2_AOVres$strain.p  = as.numeric(exp2_AOVres$strain.p)
exp2_AOVres$diet.p  = as.numeric(exp2_AOVres$diet.p)
exp2_AOVres$strain_diet  = as.numeric(exp2_AOVres$strain_diet)



exp2_AOVres$strain.p.adj <- NA
exp2_AOVres$diet.p.adj <- NA
exp2_AOVres$strain.diet.p.adj <- NA

exp2_AOVres$strain.p.adj[4:(nrow(exp2_AOVres))] <- p.adjust(as.numeric(as.character(exp2_AOVres$strain.p[4:(nrow(exp2_AOVres))])), method = 'fdr')
exp2_AOVres$diet.p.adj[4:(nrow(exp2_AOVres))] <- p.adjust(as.numeric(as.character(exp2_AOVres$diet.p[4:(nrow(exp2_AOVres))])), method = 'fdr')
exp2_AOVres$strain.diet.p.adj[4:(nrow(exp2_AOVres))] <- p.adjust(as.numeric(as.character(exp2_AOVres$strain_diet[4:(nrow(exp2_AOVres))])), method = 'fdr')


write.table(exp2_AOVres[,(ncol(exp2_AOVres)-5):ncol(exp2_AOVres)],"Exp2_AOVtrait_pval.txt", quote = F, sep = "\t")




