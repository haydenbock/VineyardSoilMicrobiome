library(tidyverse)
library(ggpubr)

MSE_Imp_DF <- readxl::read_xlsx("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/RandomForest/MSE_Output.xlsx") #raw %MSE

#rename weird column
MSE_Imp_DF <- MSE_Imp_DF %>% rename(`Experimental Factor` = `...1`)



#plot Bacterial Alpha
BA <- MSE_Imp_DF %>% 
  ggplot(aes(y = reorder(`Experimental Factor`, Bacteria_Alpha), x= Bacteria_Alpha)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     legend.position = "none") + 
  ggtitle("Factors Driving Bacterial Alpha Diversity") + xlab("%MSE Increase")



#plot Fungal Alpha
FA <- MSE_Imp_DF %>% 
  ggplot(aes(y = reorder(`Experimental Factor`, Fungi_Alpha), x= Fungi_Alpha)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     legend.position = "none") + 
  ggtitle("Factors Driving Fungal Alpha Diversity") + xlab("%MSE Increase")



#plot Bacterial Beta
BB <- MSE_Imp_DF %>% 
  ggplot(aes(y = reorder(`Experimental Factor`, Bacteria_Beta), x= Bacteria_Beta)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     legend.position = "none") + 
  ggtitle("Factors Driving Bacterial Beta Diversity") + xlab("%MSE Increase")



#plot Fungal Beta
FB <- MSE_Imp_DF %>% 
  ggplot(aes(y = reorder(`Experimental Factor`, Fungi_Beta), x= Fungi_Beta)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     legend.position = "none") + 
  ggtitle("Factors Driving Fungal Beta Diversity") + xlab("%MSE Increase")




#Figure
RF_Figure <- ggarrange(BA, FA, BB, FB, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2) 
