setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)

#load data
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")
load("nutrigenetics_phenotypeData.RData")
exp2 = exp2 %>% filter(MouseNum != 19)





exp1_BW_PC = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("PC_BW"))) %>% na.omit()
exp1_BW_PC = exp1_BW_PC %>% pivot_longer(cols=4:ncol(exp1_BW_PC), names_to = "Week", values_to = "BW_PC") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PC_BW','',WeekNum))
exp2_BW_PC = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("PC_BW")))%>% na.omit()
exp2_BW_PC = exp2_BW_PC %>% pivot_longer(cols=4:ncol(exp2_BW_PC), names_to = "Week", values_to = "BW_PC") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PC_BW','',WeekNum))
exp1_Ins = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Insulin"))) %>% na.omit()
exp1_Ins = exp1_Ins %>% pivot_longer(cols=4:ncol(exp1_Ins), names_to = "Week", values_to = "Insulin") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Insulin','',WeekNum))
exp2_Ins = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Insulin"))) %>% na.omit()
exp2_Ins = exp2_Ins %>% pivot_longer(cols=4:ncol(exp2_Ins), names_to = "Week", values_to = "Insulin") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Insulin','',WeekNum))
exp1_Glu = exp1%>% dplyr::select(c(Diet, Strain, MouseNum, contains("Glucose"))) %>% na.omit()
exp1_Glu = exp1_Glu %>% pivot_longer(cols=4:ncol(exp1_Glu), names_to = "Week", values_to = "Glucose") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Glucose','',WeekNum))
exp2_Glu = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("Glucose"))) %>% na.omit()
exp2_Glu = exp2_Glu %>% pivot_longer(cols=4:ncol(exp2_Glu), names_to = "Week", values_to = "Glucose") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Glucose','',WeekNum))
exp1_TG = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("TG"))) %>% na.omit()
exp1_TG = exp1_TG %>% pivot_longer(cols=4:ncol(exp1_TG), names_to = "Week", values_to = "TG") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_TG','',WeekNum))
exp2_TG = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("TG"))) %>% na.omit()
exp2_TG = exp2_TG %>% pivot_longer(cols=4:ncol(exp2_TG), names_to = "Week", values_to = "TG") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_TG','',WeekNum))
exp1_NEFA = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("NEFA"))) %>% na.omit()
exp1_NEFA = exp1_NEFA %>% pivot_longer(cols=4:ncol(exp1_NEFA), names_to = "Week", values_to = "NEFA") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_NEFA','',WeekNum))


exp1_BW_PC$WeekNum = as.numeric(exp1_BW_PC$WeekNum)
exp2_BW_PC$WeekNum = as.numeric(exp2_BW_PC$WeekNum)
exp1_Ins$WeekNum = as.numeric(exp1_Ins$WeekNum)
exp2_Ins$WeekNum = as.numeric(exp2_Ins$WeekNum)
exp1_Glu$WeekNum = as.numeric(exp1_Glu$WeekNum)
exp2_Glu$WeekNum = as.numeric(exp2_Glu$WeekNum)
exp1_TG$WeekNum = as.numeric(exp1_TG$WeekNum)
exp2_TG$WeekNum = as.numeric(exp2_TG$WeekNum)
exp1_NEFA$WeekNum = as.numeric(exp1_NEFA$WeekNum)
Exp1_3way_Res <- data.frame(Trait=character(),
                            strain.p=numeric(),
                            diet.p=numeric(),
                            week.p=numeric(),
                            `strain:diet`=numeric(),
                            `strain:week`=numeric(),
                            `diet:week`=numeric(),
                            `strain:diet:week`=numeric(),
                                stringsAsFactors=FALSE) 

Exp2_3way_Res <- data.frame(Trait=character(),
                            strain.p=numeric(),
                            diet.p=numeric(),
                            week.p=numeric(),
                            `strain:diet`=numeric(),
                            `strain:week`=numeric(),
                            `diet:week`=numeric(),
                            `strain:diet:week`=numeric(),
                            stringsAsFactors=FALSE)

fm <- aov(BW_PC ~ Strain * Diet*WeekNum, data =exp1_BW_PC)
Exp1_3way_Res = rbind(Exp1_3way_Res, c("Body Weight", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(BW_PC ~ Strain * Diet*WeekNum, data =exp2_BW_PC)
Exp2_3way_Res = rbind(Exp2_3way_Res, c("Body Weight", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov((Insulin) ~ Strain * Diet*WeekNum, data =exp1_Ins)
Exp1_3way_Res = rbind(Exp1_3way_Res, c("Insulin", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov((Insulin) ~ Strain * Diet*WeekNum, data =exp2_Ins)
Exp2_3way_Res = rbind(Exp2_3way_Res, c("Insulin", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(Glucose ~ Strain * Diet*Week, data =exp1_Glu)
Exp1_3way_Res = rbind(Exp1_3way_Res, c("Glucose", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(Glucose ~ Strain * Diet*Week, data =exp2_Glu)
Exp2_3way_Res = rbind(Exp2_3way_Res, c("Glucose", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(TG ~ Strain * Diet*Week, data =exp1_TG)
Exp1_3way_Res = rbind(Exp1_3way_Res,c("Triglycerides", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(TG ~ Strain * Diet*Week, data =exp2_TG)
Exp2_3way_Res = rbind(Exp2_3way_Res,c("Triglycerides", summary(fm)[[1]]$`Pr(>F)`[1:7]))
fm <- aov(NEFA ~ Strain * Diet*Week, data =exp1_NEFA)
Exp1_3way_Res = rbind(Exp1_3way_Res, c("NEFA", summary(fm)[[1]]$`Pr(>F)`[1:7]))


colnames(Exp1_3way_Res) = c("Trait","strain.p","diet.p","week.p","strain:diet.p","strain:week.p","diet:week.p","strain:diet:week.p")
colnames(Exp2_3way_Res) = c("Trait","strain.p","diet.p","week.p","strain:diet.p","strain:week.p","diet:week.p","strain:diet:week.p")

Exp1_3way_Res$strain.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$strain.p)), method = 'fdr')
Exp1_3way_Res$diet.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$diet.p)), method = 'fdr')
Exp1_3way_Res$week.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$week.p)), method = 'fdr')
Exp1_3way_Res$straindiet.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$`strain:diet.p`)), method = 'fdr')
Exp1_3way_Res$dietweek.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$`diet:week.p`)), method = 'fdr')
Exp1_3way_Res$strainweek.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$`strain:week.p`)), method = 'fdr')
Exp1_3way_Res$straindietweek.padj = p.adjust(as.numeric(as.character(Exp1_3way_Res$`strain:diet:week.p`)), method = 'fdr')
Exp2_3way_Res$strain.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$strain.p)), method = 'fdr')
Exp2_3way_Res$dietweek.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$`diet:week.p`)), method = 'fdr')
Exp2_3way_Res$week.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$week.p)), method = 'fdr')
Exp2_3way_Res$strainweek.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$`strain:week.p`)), method = 'fdr')
Exp2_3way_Res$diet.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$diet.p)), method = 'fdr')
Exp2_3way_Res$straindiet.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$`strain:diet.p`)), method = 'fdr')
Exp2_3way_Res$straindietweek.padj = p.adjust(as.numeric(as.character(Exp2_3way_Res$`strain:diet:week.p`)), method = 'fdr')
write.table(Exp1_3way_Res,"Exp1_3way_ANOVAtrait_pval_outlier_removed.txt", quote = F, sep = "\t")
write.table(Exp2_3way_Res,"Exp2_3way_ANOVAtrait_pval_outlier_removed.txt", quote = F, sep = "\t")






