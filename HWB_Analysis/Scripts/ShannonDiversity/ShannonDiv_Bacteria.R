#Bacterial Shannon diversity calc for bulk soil ONLY
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
              "dplyr","MuMIn","Biostrings","ShortRead","ape","ade4","adegraphics","adespatial","tmap","dada2","jtools","rcompanion","ggnewscale","RColorBrewer","emmeans")
ipak(packages)


##load augmented PS and spatial data ##
load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")





#ASV LEVEL ----
####Select PS objects for this iteration ####
ps.b.19<-ps.16s.2019 #bacteria 2019
ps.b.20<-ps.16s.2020 #bacteria 2020


#Subset soil only compartment
ps.soil.b.19 <-subset_samples(ps.b.19, type=="soil")
ps.soil.b.20 <-subset_samples(ps.b.20, type=="soil")

#RAREFY
sample_sums(ps.soil.b.19)
sample_sums(ps.soil.b.20)


ps.soil.b.19<-prune_samples(sample_sums(ps.soil.b.19)>0, ps.soil.b.19)
ps.soil.b.20<-prune_samples(sample_sums(ps.soil.b.20)>0, ps.soil.b.20)


set.seed(500)
ps.soil.rare.b.19<-rarefy_even_depth(ps.soil.b.19)
sample_sums(ps.soil.rare.b.19)

ps.soil.rare.b.20<-rarefy_even_depth(ps.soil.b.20)
sample_sums(ps.soil.rare.b.20)
####breakdown ps rare into components####
#2019
counts.b.19<-as.data.frame(ps.soil.rare.b.19@otu_table)
meta.b.19<-as.data.frame(ps.soil.rare.b.19@sam_data)
meta.b.19$depth<-as.factor(meta.b.19$depth)

#2020
counts.b.20<-as.data.frame(ps.soil.rare.b.20@otu_table)
meta.b.20<-as.data.frame(ps.soil.rare.b.20@sam_data)
meta.b.20$depth<-as.factor(meta.b.20$depth)


#### Calculate alpha diversity ####
#2019
rich.b.19 <- estimate_richness(ps.soil.rare.b.19,measures="Shannon")
richtable.b.19<-cbind(meta.b.19,rich.b.19)
richtable.b.19$interact<-paste(richtable.b.19$gc,richtable.b.19$rs)

#2020
rich.b.20 <- estimate_richness(ps.soil.rare.b.20,measures="Shannon")
richtable.b.20<-cbind(meta.b.20,rich.b.20)
richtable.b.20$interact<-paste(richtable.b.20$gc,richtable.b.20$rs)


#### Examine for normality ####
#2019
qqnorm(richtable.b.19$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.19$Shannon, col = "steelblue", lwd = 2) #Not that great looking...

#2020
qqnorm(richtable.b.20$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.20$Shannon, col = "steelblue", lwd = 2) #Decent!


#### Model Building ####
#bind 2019 and 2020 dataframes 
richtable.b.19 <- richtable.b.19 %>% select(-Row, -Column)
richtable.b.20 <- richtable.b.20 %>% select(-newroot)
SoilBac <- rbind(richtable.b.19, richtable.b.20)
#SAVE CSV
write.csv(SoilBac, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

SoilBac <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

#Full model
Full_Model <- lme(Shannon ~ depth * gc * rs * year, random = ~ 1 | blk/rs, data = SoilBac)
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


#Bacteria best model is Shannon ~ depth + year + depth*year + ~ 1 | blk/rs


#GENERA LEVEL
####Select PS objects for this iteration ####
ps.b.19<-ps.16s.2019.gen #bacteria 2019
ps.b.20<-ps.16s.2020.gen #bacteria 2020


#Subset soil only compartment
ps.soil.b.19 <-subset_samples(ps.b.19, type=="soil")
ps.soil.b.20 <-subset_samples(ps.b.20, type=="soil")

#RAREFY
sample_sums(ps.soil.b.19)
sample_sums(ps.soil.b.20)


ps.soil.b.19<-prune_samples(sample_sums(ps.soil.b.19)>0, ps.soil.b.19)
ps.soil.b.20<-prune_samples(sample_sums(ps.soil.b.20)>0, ps.soil.b.20)


set.seed(500)
ps.soil.rare.b.19<-rarefy_even_depth(ps.soil.b.19)
sample_sums(ps.soil.rare.b.19)

ps.soil.rare.b.20<-rarefy_even_depth(ps.soil.b.20)
sample_sums(ps.soil.rare.b.20)

####breakdown ps rare into components####
#2019
counts.b.19<-as.data.frame(ps.soil.rare.b.19@otu_table)
meta.b.19<-as.data.frame(ps.soil.rare.b.19@sam_data)
meta.b.19$depth<-as.factor(meta.b.19$depth)

#2020
counts.b.20<-as.data.frame(ps.soil.rare.b.20@otu_table)
meta.b.20<-as.data.frame(ps.soil.rare.b.20@sam_data)
meta.b.20$depth<-as.factor(meta.b.20$depth)


#### Calculate alpha diversity ####
#2019
rich.b.19 <- estimate_richness(ps.soil.rare.b.19,measures="Shannon")
richtable.b.19<-cbind(meta.b.19,rich.b.19)
richtable.b.19$interact<-paste(richtable.b.19$gc,richtable.b.19$rs)

#2020
rich.b.20 <- estimate_richness(ps.soil.rare.b.20,measures="Shannon")
richtable.b.20<-cbind(meta.b.20,rich.b.20)
richtable.b.20$interact<-paste(richtable.b.20$gc,richtable.b.20$rs)


#### Examine for normality ####
#2019
qqnorm(richtable.b.19$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.19$Shannon, col = "steelblue", lwd = 2) #Not that great looking...

#2020
qqnorm(richtable.b.20$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.20$Shannon, col = "steelblue", lwd = 2) #Decent!


#### Model Building ####
#bind 2019 and 2020 dataframes 
richtable.b.19 <- richtable.b.19 %>% select(-Row, -Column)
richtable.b.20 <- richtable.b.20 %>% select(-newroot)
SoilBac <- rbind(richtable.b.19, richtable.b.20)
#SAVE CSV
write.csv(SoilBac, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

SoilBac <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

#Full model
Full_Model <- lme(Shannon ~ depth * gc * rs * year, random = ~ 1 | blk/rs, data = SoilBac)
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


#Bacteria best model is Shannon ~ depth + year + depth*year + ~ 1 | blk/rs



#PHYLUM LEVEL ----
#RAREFY
sample_sums(ps.soil.b.19)
sample_sums(ps.soil.b.20)


ps.soil.b.19<-prune_samples(sample_sums(ps.soil.b.19)>0, ps.soil.b.19)
ps.soil.b.20<-prune_samples(sample_sums(ps.soil.b.20)>0, ps.soil.b.20)


set.seed(500)
ps.soil.rare.b.19<-rarefy_even_depth(ps.soil.b.19)
sample_sums(ps.soil.rare.b.19)

ps.soil.rare.b.20<-rarefy_even_depth(ps.soil.b.20)
sample_sums(ps.soil.rare.b.20)

####Select PS objects for this iteration ####
ps.b.19<-ps.16s.2019.gen #bacteria 2019
ps.b.20<-ps.16s.2020.gen #bacteria 2020


#Subset soil only compartment
ps.soil.b.19 <-subset_samples(ps.b.19, type=="soil")
ps.soil.b.20 <-subset_samples(ps.b.20, type=="soil")
#RAREFY
sample_sums(ps.soil.b.19)
sample_sums(ps.soil.b.20)


ps.soil.b.19<-prune_samples(sample_sums(ps.soil.b.19)>0, ps.soil.b.19)
ps.soil.b.20<-prune_samples(sample_sums(ps.soil.b.20)>0, ps.soil.b.20)


set.seed(500)
ps.soil.rare.b.19<-rarefy_even_depth(ps.soil.b.19)
sample_sums(ps.soil.rare.b.19)

ps.soil.rare.b.20<-rarefy_even_depth(ps.soil.b.20)
sample_sums(ps.soil.rare.b.20)
####breakdown ps rare into components####
#2019
counts.b.19<-as.data.frame(ps.soil.rare.b.19@otu_table)
meta.b.19<-as.data.frame(ps.soil.rare.b.19@sam_data)
meta.b.19$depth<-as.factor(meta.b.19$depth)

#2020
counts.b.20<-as.data.frame(ps.soil.rare.b.20@otu_table)
meta.b.20<-as.data.frame(ps.soil.rare.b.20@sam_data)
meta.b.20$depth<-as.factor(meta.b.20$depth)



#### Calculate alpha diversity ####
#2019
rich.b.19 <- estimate_richness(ps.soil.rare.b.19,measures="Shannon")
richtable.b.19<-cbind(meta.b.19,rich.b.19)
richtable.b.19$interact<-paste(richtable.b.19$gc,richtable.b.19$rs)

#2020
rich.b.20 <- estimate_richness(ps.soil.rare.b.20,measures="Shannon")
richtable.b.20<-cbind(meta.b.20,rich.b.20)
richtable.b.20$interact<-paste(richtable.b.20$gc,richtable.b.20$rs)


#### Examine for normality ####
#2019
qqnorm(richtable.b.19$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.19$Shannon, col = "steelblue", lwd = 2) #Not that great looking...

#2020
qqnorm(richtable.b.20$Shannon, pch = 1, frame = FALSE)
qqline(richtable.b.20$Shannon, col = "steelblue", lwd = 2) #Decent!


#### Model Building ####
#bind 2019 and 2020 dataframes 
richtable.b.19 <- richtable.b.19 %>% select(-Row, -Column)
richtable.b.20 <- richtable.b.20 %>% select(-newroot)
SoilBac <- rbind(richtable.b.19, richtable.b.20)
#SAVE CSV
write.csv(SoilBac, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

SoilBac <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilBacALPHA.csv")

#Full model
Full_Model <- lme(Shannon ~ depth * gc * rs * year, random = ~ 1 | blk/rs, data = SoilBac)
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


#Bacteria best model is Shannon ~ depth + year + depth*year + ~ 1 | blk/rs
