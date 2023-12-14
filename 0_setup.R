setwd("~/Civelek Lab/RIMBANET/Final Net Data/REDO MET SEQ/STARandGTExNetsFinal/Newest")


library(dplyr)
library(ggplot2)
library(tidyr)


exp1_phen = read.csv("nutri_exp1_nofood.csv", fill = T, header = T) %>% filter(Diet != "Regular Chow")
exp2_phen = read.csv("nutri_exp2_no_food2.csv", fill = T, header = T)
exp2_organ = read.table("organweights1024.txt", fill = T, header = T) %>% pivot_wider(names_from = "Organ", values_from = "OrganWeight") %>% dplyr::select(!mouseWeight) %>% dplyr::select(!intestine)
exp2_full = full_join(exp2_phen, exp2_organ, by = c("Strain", "Diet", "MouseNum"))





exp2_added = exp2_full %>% mutate(VATw_per_BW = VAT/Week_6_BW)
exp2_added = exp2_added %>% mutate(BATw_per_BW = BAT/Week_6_BW)
exp2_added = exp2_added %>% mutate(SATw_per_BW = SAT/Week_6_BW)
exp2_added = exp2_added %>% mutate(Liverw_per_BW = Liver/Week_6_BW)
exp2_added = exp2_added %>% mutate(Musclew_per_BW = Muscle/Week_6_BW)
exp2_added = exp2_added %>% mutate(brainw_per_BW = brain/Week_6_BW)
exp2_added = exp2_added %>% mutate(gatrow_per_BW = gastro/Week_6_BW)
exp2_added = exp2_added %>% mutate(pancreasw_per_BW = pancreas/Week_6_BW)

exp2_added = exp2_added %>% mutate(Week_1_PC_BW = 100*((Week_1_BW - Week_0_BW)/Week_0_BW))
exp2_added = exp2_added %>% mutate(Week_2_PC_BW = 100*((Week_2_BW - Week_0_BW)/Week_0_BW))
exp2_added = exp2_added %>% mutate(Week_3_PC_BW = 100*((Week_3_BW - Week_0_BW)/Week_0_BW))
exp2_added = exp2_added %>% mutate(Week_4_PC_BW = 100*((Week_4_BW - Week_0_BW)/Week_0_BW))
exp2_added = exp2_added %>% mutate(Week_5_PC_BW = 100*((Week_5_BW - Week_0_BW)/Week_0_BW))
exp2_added = exp2_added %>% mutate(Week_6_PC_BW = 100*((Week_6_BW - Week_0_BW)/Week_0_BW))


exp1_added = exp1_phen %>% mutate(Week_1_PC_BW = 100*((Week_1_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_2_PC_BW = 100*((Week_2_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_3_PC_BW = 100*((Week_3_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_4_PC_BW = 100*((Week_4_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_5_PC_BW = 100*((Week_5_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_6_PC_BW = 100*((Week_6_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_7_PC_BW = 100*((Week_7_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_8_PC_BW = 100*((Week_8_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_9_PC_BW = 100*((Week_9_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_10_PC_BW = 100*((Week_10_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))
exp1_added = exp1_added %>% mutate(Week_16_PC_BW = 100*((Week_16_BW_Mass - Week_0_BW_Mass)/Week_0_BW_Mass))

exp1 = exp1_added %>% select(!Week_1_BW_Mass) %>% select(!Week_2_BW_Mass) %>% select(!Week_3_BW_Mass) %>% select(!Week_4_BW_Mass) %>% select(!Week_5_BW_Mass) %>% select(!Week_6_BW_Mass) %>% select(!Week_7_BW_Mass) %>% select(!Week_8_BW_Mass) %>% select(!Week_9_BW_Mass) %>% select(!Week_10_BW_Mass) %>% select(!Week_16_BW_Mass) %>% filter (MouseNum != 94)

exp2 = exp2_added %>% select(!Week_1_BW) %>% select(!Week_2_BW) %>% select(!Week_3_BW) %>% select(!Week_4_BW) %>% select(!Week_5_BW) %>% select(!Week_6_BW) %>% select(!VAT)%>% select(!SAT)%>% select(!BAT)%>% select(!pancreas)%>% select(!gastro)%>% select(!brain)%>% select(!Liver)%>% select(!Muscle) %>% filter (MouseNum != 14)

exp1$Diet = factor(exp1$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2$Diet = factor(exp2$Diet, levels = c("American","Mediterranean","Vegetarian","Vegan"))
exp2$Strain = factor(exp2$Strain)
exp1$Strain = factor(exp1$Strain)

rm(exp2_added, exp2_full, exp2_organ, exp2_phen, exp1_added, exp1_phen)
save.image("nutrigenetics_phenotypeData_both_outlier_removed.RData")





