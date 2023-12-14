# graphing
setwd("~/")

library(ggplot2)
library(dplyr)
library(tidyr)

#load data
load("nutrigenetics_phenotypeData_both_outlier_removed.RData")
load("nutrigenetics_phenotypeData.RData")
exp2 = exp2 %>% filter(MouseNum != 19)

#time series data
############
colo1 = c("grey","slateblue","lightpink","red3")
colo2 = c("yellow","orange","forestgreen","blue")

exp1_BW_PC = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("PC_BW")))
exp2_BW_PC = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("PC_BW")))
exp1_Ins = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Insulin")))
exp2_Ins = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Insulin")))
exp1_Glu = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Glucose")))
exp2_Glu = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("Glucose")))
exp1_TG = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("TG")))
exp2_TG = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("TG")))
exp1_NEFA = exp1 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("NEFA")))

exp1_BW_PC1= exp1_BW_PC %>% pivot_longer(cols=4:ncol(exp1_BW_PC), names_to = "Week", values_to = "BW_PC") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PC_BW','',WeekNum)) %>% na.omit()
exp1_BW_PC1$WeekNum = factor(exp1_BW_PC1$WeekNum, levels = c("1","2","3","4","5","6","7","8","9","10","16"))
exp1_BW_PC1$Diet = factor(exp1_BW_PC1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp2_BW_PC1= exp2_BW_PC %>% pivot_longer(cols=4:ncol(exp2_BW_PC), names_to = "Week", values_to = "BW_PC") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PC_BW','',WeekNum)) %>% na.omit()
exp2_BW_PC1$WeekNum = factor(exp2_BW_PC1$WeekNum, levels = c("1","2","3","4","5","6"))
exp2_BW_PC1$Diet = factor(exp2_BW_PC1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp1_Ins1= exp1_Ins %>% pivot_longer(cols=4:ncol(exp1_Ins), names_to = "Week", values_to = "Insulin") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Insulin','',WeekNum)) %>% na.omit()
exp1_Ins1$WeekNum = factor(exp1_Ins1$WeekNum, levels = c("0","4","8","16"))
exp1_Ins1$Diet = factor(exp1_Ins1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp2_Ins1= exp2_Ins %>% pivot_longer(cols=4:ncol(exp2_Ins), names_to = "Week", values_to = "Insulin") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Insulin','',WeekNum)) %>% na.omit()
exp2_Ins1$WeekNum = factor(exp2_Ins1$WeekNum, levels = c("0","6"))
exp2_Ins1$Diet = factor(exp2_Ins1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp1_TG1= exp1_TG %>% pivot_longer(cols=4:ncol(exp1_TG), names_to = "Week", values_to = "TG") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_TG','',WeekNum)) %>% na.omit()
exp1_TG1$WeekNum = factor(exp1_TG1$WeekNum, levels = c("0","4","8","16"))
exp1_TG1$Diet = factor(exp1_TG1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp2_TG1= exp2_TG %>% pivot_longer(cols=4:ncol(exp2_TG), names_to = "Week", values_to = "TG") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_TG','',WeekNum)) %>% na.omit()
exp2_TG1$WeekNum = factor(exp2_TG1$WeekNum, levels = c("0","6"))
exp2_TG1$Diet = factor(exp2_TG1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp1_Glu1= exp1_Glu %>% pivot_longer(cols=4:ncol(exp1_Glu), names_to = "Week", values_to = "Glucose") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Glucose','',WeekNum)) %>% na.omit()
exp1_Glu1$WeekNum = factor(exp1_Glu1$WeekNum, levels = c("0","4","8","16"))
exp1_Glu1$Diet = factor(exp1_Glu1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp2_Glu1= exp2_Glu %>% pivot_longer(cols=4:ncol(exp2_Glu), names_to = "Week", values_to = "Glucose") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_Glucose','',WeekNum)) %>% na.omit()
exp2_Glu1$WeekNum = factor(exp2_Glu1$WeekNum, levels = c("0","6"))
exp2_Glu1$Diet = factor(exp2_Glu1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

exp1_NEFA1= exp1_NEFA %>% pivot_longer(cols=4:ncol(exp1_NEFA), names_to = "Week", values_to = "NEFA") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_NEFA','',WeekNum)) %>% na.omit()
exp1_NEFA1$WeekNum = factor(exp1_NEFA1$WeekNum, levels = c("0","4","8","16"))
exp1_NEFA1$Diet = factor(exp1_NEFA1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))


BW_PC_graph1 = exp1_BW_PC1  %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(BW_PC),sePC = sd(BW_PC)/sqrt(length(BW_PC))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_line(aes(color = Strain))+
  geom_point(data =exp1_BW_PC1, aes(x = WeekNum, y = BW_PC, color = Strain),size = .4, alpha = .5)+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5,color = "grey20", alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5),legend.position = "none")+
  xlab("Week")+
  ylab("% Change Body Weight")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+facet_wrap(~Diet)

BW_PC_graph1
pdf("Exp1_BW_PC_graph.pdf",width = 5.5, height = 4)
BW_PC_graph1
dev.off()

BW_PC_graph2 = exp2_BW_PC1  %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(BW_PC),sePC = sd(BW_PC)/sqrt(length(BW_PC))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_line(aes(color = Strain))+
  geom_point(data =exp2_BW_PC1, aes(x = WeekNum, y = BW_PC, color = Strain),size = .4, alpha = .5)+
  geom_errorbar(aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5,color = "grey20", alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5),legend.position = "none")+
  xlab("Week")+
  ylab("% Change Body Weight")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+facet_wrap(~Diet)

BW_PC_graph2
pdf("Exp2_BW_PC_graph.pdf",width = 5.5, height = 4)
BW_PC_graph2
dev.off()



TG_graph = exp1_TG1 %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(TG),sePC = sd(TG)/sqrt(length(TG))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width = .4, aes(fill = Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .25)+
  geom_point(data =exp1_TG1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = TG, group= Strain),size = .3,alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Triglycerides (mg/dL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, nrow = 1)
TG_graph

pdf("TG_strain.pdf", width = 6, height=2)
TG_graph
dev.off()




TG_graph2 = exp2_TG1 %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(TG),sePC = sd(TG)/sqrt(length(TG))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width = .4, aes(fill = Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .25)+
  geom_point(data =exp2_TG1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = TG, group= Strain),size = .3,alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Triglycerides (mg/dL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, nrow = 1)
TG_graph2

pdf("Exp2_TGs.pdf", width = 6, height=2)
TG_graph2
dev.off()



Ins_graph = exp1_Ins1 %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean((Insulin)),sePC = sd((Insulin))/sqrt(length((Insulin)))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = .5,aes(fill = Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp1_Ins1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = (Insulin), group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Insulin (ng/mL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, ncol = 4)+ggbreak::scale_y_cut(10)+scale_y_continuous(breaks = c(0,4,8,50,100))

Ins_graph

cairo_pdf("Ins_strain.pdf", width = 6.15, height=2)
Ins_graph
dev.off()

Ins_graph2 = exp2_Ins1 %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean((Insulin)),sePC = sd((Insulin))/sqrt(length((Insulin)))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = .5,aes(fill = Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_Ins1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = (Insulin), group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Insulin (ng/mL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, ncol = 4)+ggbreak::scale_y_cut(10)+scale_y_continuous(breaks = c(0,4,8,50,100))

Ins_graph2

cairo_pdf("Ins_exp2.pdf", width = 6.15, height=2)
Ins_graph2
dev.off()



NEFA_graph = exp1_NEFA1 %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(NEFA),sePC = sd(NEFA)/sqrt(length(NEFA))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp1_NEFA1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = NEFA, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("NEFA (mEq/L)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, ncol = 4)
NEFA_graph


Glu_graph = exp1_Glu1  %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(Glucose),sePC = sd(Glucose)/sqrt(length(Glucose))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp1_Glu1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = Glucose, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Glucose (mg/dL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, ncol = 4)
Glu_graph

pdf("GLu_strain.pdf", width = 6, height=2)
Glu_graph
dev.off()

Glu_graph2 = exp2_Glu1  %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(Glucose),sePC = sd(Glucose)/sqrt(length(Glucose))) %>% 
  ggplot(aes(x= WeekNum, y = meanPC, group = Strain))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_Glu1, position = position_dodge(width = 0.8),aes(x = WeekNum, y = Glucose, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none")+
  xlab("Week")+
  ylab("Glucose (mg/dL)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)+ 
  facet_wrap(~Diet, ncol = 4)
Glu_graph2

pdf("GLu_exp2.pdf", width = 6, height=2)
Glu_graph2
dev.off()



##########
#single timepoint data
#######
#single timepoint exp2
exp2_echo = exp2 %>% dplyr::select(c(Diet, Strain, MouseNum, contains("MRI"))) 

exp2_echo = exp2_echo %>% pivot_longer(cols=4:ncol(exp2_echo), names_to = "Week", values_to = "MRI") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_MRI','',WeekNum))%>% na.omit()
exp2_echo$Diet = factor(exp2_echo$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

MRI_graph = exp2_echo  %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(MRI),sePC = sd(MRI)/sqrt(length(MRI))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = .5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_echo, position = position_dodge(width = 0.8),aes(x = Diet, y = MRI, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("Total Fat MRI (volume/BW)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
MRI_graph
pdf("MRI_exp2.pdf", width = 3, height=2.5)
MRI_graph
dev.off()




exp2_PVAT = exp2_added %>% dplyr::select(c(Diet, Strain, MouseNum, contains("PET_VAT"))) 

exp2_PVAT = exp2_PVAT %>% pivot_longer(cols=4:ncol(exp2_PVAT), names_to = "Week", values_to = "PET") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PET_VAT','',WeekNum))%>% na.omit()
exp2_PVAT$Diet = factor(exp2_PVAT$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

PETVAT_graph = exp2_PVAT %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(PET),sePC = sd(PET)/sqrt(length(PET))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = .5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_PVAT, position = position_dodge(width = 0.8),aes(x = Diet, y = PET, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("FDG-Uptake VAT (SUV %ID)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
PETVAT_graph

pdf("PEtVAT_strain.pdf", width = 3, height=2.5)
PETVAT_graph
dev.off()

exp2_PBAT = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("PET_BAT"))) 
exp2_PBAT = exp2_PBAT %>% pivot_longer(cols=4:ncol(exp2_PBAT), names_to = "Week", values_to = "PET") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PET_BAT','',WeekNum)) %>% na.omit()
exp2_PBAT$Diet = factor(exp2_PBAT$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

PETBAT_graph = exp2_PBAT %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(PET),sePC = sd(PET)/sqrt(length(PET))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_PBAT, position = position_dodge(width = 0.8),aes(x = Diet, y = PET, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("FDG-Uptake BAT (SUV %ID)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
PETBAT_graph
pdf("PEtBAT_strain.pdf", width = 3, height=2.5)
PETBAT_graph
dev.off()


exp2_PQuad = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("PET_Quad")))

exp2_PQuad = exp2_PQuad %>% pivot_longer(cols=4:ncol(exp2_PQuad), names_to = "Week", values_to = "PET") %>% mutate(WeekNum = gsub('Week_','',Week)) %>% mutate(WeekNum = gsub('_PET_Quad','',WeekNum))%>% na.omit()
exp2_PQuad$Diet = factor(exp2_PQuad$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

PETQuad_graph = exp2_PQuad %>% group_by(Strain,Diet,WeekNum) %>% summarise(meanPC = mean(PET),sePC = sd(PET)/sqrt(length(PET))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = .5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_PQuad, position = position_dodge(width = 0.8),aes(x = Diet, y = PET, group= Strain),size = .3, alpha = .8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("FDG-Uptake Quad (SUV %ID)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
PETQuad_graph
pdf("PEtQuad_strain.pdf", width = 3, height=2.5)
PETQuad_graph
dev.off()




exp2_VATw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("VATw"))) %>% na.omit()
exp2_VATw$Diet = factor(exp2_VATw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_BATw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("BATw"))) %>% na.omit()
exp2_BATw$Diet = factor(exp2_BATw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_SATw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("SATw"))) %>% na.omit()
exp2_SATw$Diet = factor(exp2_SATw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_Liverw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("Liverw"))) %>% na.omit()
exp2_Liverw$Diet = factor(exp2_Liverw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_Musclew = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("Musclew"))) %>% na.omit()
exp2_Musclew$Diet = factor(exp2_Musclew$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_pancreasw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("pancreasw"))) %>% na.omit()
exp2_pancreasw$Diet = factor(exp2_pancreasw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_brainw = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("brainw"))) %>% na.omit()
exp2_brainw$Diet = factor(exp2_brainw$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2_gastrow = exp2%>% dplyr::select(c(Diet, Strain, MouseNum, contains("gatrow"))) %>% na.omit()
colnames(exp2_gastrow)[4] = "gastrow_per_BW"
exp2_gastrow$Diet = factor(exp2_gastrow$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))

BATw_graph = exp2_BATw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(BATw_per_BW),sePC = sd(BATw_per_BW)/sqrt(length(BATw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_BATw, position = position_dodge(width = 0.8),aes(x = Diet, y = BATw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("BAT Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
BATw_graph

pdf("BATwperBW_strain.pdf", width = 2.6, height=2.2)
BATw_graph
dev.off()

VATw_graph = exp2_VATw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(VATw_per_BW),sePC = sd(VATw_per_BW)/sqrt(length(VATw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_VATw, position = position_dodge(width = 0.8),aes(x = Diet, y = VATw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("VAT Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
VATw_graph

pdf("VATwperBW_strain.pdf", width = 2.6, height=2.2)
VATw_graph
dev.off()

SATw_graph = exp2_SATw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(SATw_per_BW),sePC = sd(SATw_per_BW)/sqrt(length(SATw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_SATw, position = position_dodge(width = 0.8),aes(x = Diet, y = SATw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("SAT Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
SATw_graph

pdf("SATwperBW_strain.pdf", width = 2.6, height=2.2)
SATw_graph
dev.off()

Liverw_graph = exp2_Liverw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(Liverw_per_BW),sePC = sd(Liverw_per_BW)/sqrt(length(Liverw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_Liverw, position = position_dodge(width = 0.8),aes(x = Diet, y = Liverw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("Liver Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
Liverw_graph

pdf("LiverwperBW_strain.pdf", width = 2.6, height=2.2)
Liverw_graph
dev.off()

Musclew_graph = exp2_Musclew %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(Musclew_per_BW),sePC = sd(Musclew_per_BW)/sqrt(length(Musclew_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_Musclew, position = position_dodge(width = 0.8),aes(x = Diet, y = Musclew_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("Muscle Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
Musclew_graph

pdf("MusclewperBW_strain.pdf", width = 2.6, height=2.2)
Musclew_graph
dev.off()

pancreasw_graph = exp2_pancreasw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(pancreasw_per_BW),sePC = sd(pancreasw_per_BW)/sqrt(length(pancreasw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_pancreasw, position = position_dodge(width = 0.8),aes(x = Diet, y = pancreasw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("pancreas Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
pancreasw_graph

pdf("pancreaswperBW_strain.pdf", width = 2.6, height=2.2)
pancreasw_graph
dev.off()

brainw_graph = exp2_brainw %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(brainw_per_BW),sePC = sd(brainw_per_BW)/sqrt(length(brainw_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_brainw, position = position_dodge(width = 0.8),aes(x = Diet, y = brainw_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("brain Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
brainw_graph

pdf("brainwperBW_strain.pdf", width = 2.6, height=2.2)
brainw_graph
dev.off()

gastrow_graph = exp2_gastrow %>% group_by(Strain,Diet) %>% summarise(meanPC = mean(gastrow_per_BW),sePC = sd(gastrow_per_BW)/sqrt(length(gastrow_per_BW))) %>% 
  ggplot(aes(x= Diet, y = meanPC, group = Strain))+
  geom_rect(xmin = 3.6, xmax = 4.4, ymin = 0, ymax = 0, fill = colo2[4], alpha = 0.2)+
  annotate("rect", xmin=0.6, xmax=1.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[1])+
  
  annotate("rect", xmin=1.6, xmax=2.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[2])+
  annotate("rect", xmin=2.6, xmax=3.4, ymin=0, ymax=Inf, alpha=0.6, fill=colo2[3])+
  
  annotate("rect", xmin=3.6, xmax=4.4, ymin=0, ymax=Inf, alpha=0.3, fill=colo2[4])+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),width =.5, aes(fill =Strain,color = Strain))+
  geom_errorbar(position = position_dodge(width = 0.8),aes(ymin=meanPC-sePC, ymax=meanPC+sePC), width=.15, size = .5)+
  geom_point(data =exp2_gastrow, position = position_dodge(width = 0.8),aes(x = Diet, y = gastrow_per_BW, group= Strain),size = .3, alpha =.8)+
  theme_classic()+
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5), legend.position = "none",axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1) )+
  xlab("Diet")+
  ylab("gastro Weight (g) / BW (g)")+ 
  scale_color_manual(values = colo1)+ 
  scale_fill_manual(values = colo1)
gastrow_graph

pdf("gastrowperBW_strain.pdf", width = 2.6, height=2.2)
gastrow_graph
dev.off()



