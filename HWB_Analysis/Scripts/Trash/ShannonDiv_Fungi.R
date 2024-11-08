#Fungal Shannon diversity calc for bulk soil ONLY
#Created by HWB for ccxrs microbiome analysis

#### Setup and Data Import ####

setwd("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs")

### Clear workspace ###
rm(list=ls())

## load packages we may need
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("gridExtra","nlme","lme4", "effects","corncob","magrittr","car", "interactions","RVAideMemoire","tibble","doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq",
              "dplyr","Biostrings","ShortRead","ape","ade4","adegraphics","adespatial","tmap","dada2","jtools","rcompanion","ggnewscale","RColorBrewer","emmeans", "MuMIn")
ipak(packages)



##load augmented PS and spatial data ##
load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")



####Select PS objects for this iteration ####
ps.f.19<-ps.its.2019.phy #fungi 2019
ps.f.20<-ps.its.2020.phy #fungi 2020

#Subset soil only compartment
ps.soil.f.19 <-subset_samples(ps.f.19, type=="soil")
ps.soil.f.20 <-subset_samples(ps.f.20, type=="soil")

#RAREFY
sample_sums(ps.soil.f.19)
sample_sums(ps.soil.f.20)

ps.soil.f.19<-prune_samples(sample_sums(ps.soil.f.19)>0, ps.soil.f.19)
ps.soil.f.20<-prune_samples(sample_sums(ps.soil.f.20)>0, ps.soil.f.20)

set.seed(500)
ps.soil.rare.f.19<-rarefy_even_depth(ps.soil.f.19)
sample_sums(ps.soil.rare.f.19)

ps.soil.rare.f.20<-rarefy_even_depth(ps.soil.f.20)
sample_sums(ps.soil.rare.f.20) 


####breakdown ps rare into components####
#2019
counts.f.19<-as.data.frame(ps.soil.rare.f.19@otu_table)
meta.f.19<-as.data.frame(ps.soil.rare.f.19@sam_data)
meta.f.19$depth<-as.factor(meta.f.19$depth)

#2020
counts.f.20<-as.data.frame(ps.soil.rare.f.20@otu_table)
meta.f.20<-as.data.frame(ps.soil.rare.f.20@sam_data)
meta.f.20$depth<-as.factor(meta.f.20$depth)

#### Calculate alpha diversity ####
#2019
rich.f.19 <- estimate_richness(ps.soil.rare.f.19,measures="Shannon")
richtable.f.19<-cbind(meta.f.19,rich.f.19)
richtable.f.19$interact<-paste(richtable.f.19$gc,richtable.f.19$rs)

#2020
rich.f.20 <- estimate_richness(ps.soil.rare.f.20,measures="Shannon")
richtable.f.20<-cbind(meta.f.20,rich.f.20)
richtable.f.20$interact<-paste(richtable.f.20$gc,richtable.f.20$rs)


#### Examine for normality ####
#2019
qqnorm(richtable.f.19$Shannon, pch = 1, frame = FALSE)
qqline(richtable.f.19$Shannon, col = "steelblue", lwd = 2) #Not that great looking...

#2020
qqnorm(richtable.f.20$Shannon, pch = 1, frame = FALSE)
qqline(richtable.f.20$Shannon, col = "steelblue", lwd = 2) #Decent!



#### Model Building ####
#Data Prep
#bind 2019 and 2020 dataframes 
richtable.f.19 <- richtable.f.19 %>% select(-Row, -Column, -sample2)
richtable.f.20 <- richtable.f.20 %>% select(-newroot)
SoilFungi <- rbind(richtable.f.19, richtable.f.20)

#Found out that lme freaks out if your response variable (shannon div) column has any "0" values. So, lets remove...
SoilFungi <- SoilFungi[SoilFungi$Shannon != 0, ]
#SAVE CSV
write.csv(SoilFungi, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilFungiALPHA.csv")


#Full model
Full_Model <- lme(Shannon ~ depth * gc * rs * year, random = ~ 1 | blk/rs, data = SoilFungi)
summary(Full_Model)


# Generate all possible combinations of fixed effects and compare models using AIC
all_models <- dredge(Full_Model)

# View the best models ranked by AIC
print(all_models)

# Extract the best models within 2 AIC units of the top model
best_models_list <- get.models(all_models, subset = delta < 2)

# View the number of best models found
length(best_models_list)

# Loop through the best models and print the summary of each
for (i in seq_along(best_models_list)) {
  cat("\nSummary of Best Model", i, ":\n")
  print(summary(best_models_list[[i]]))
}

# Once you have looked at the model summaries, select the model you want from the best model list
best_model_1 <- best_models_list[[1]]  # Extract first best model
summary(best_model_1)  # summary
anova(best_model_1)    #ANOVA table for the model


#Fungi best model is Shannon ~ depth + year + depth*year + ~ 1 | blk/rs