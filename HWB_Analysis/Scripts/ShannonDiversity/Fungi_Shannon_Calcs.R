
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


#ASV ----
  ####Select PS objects for this iteration ####
ps.f.19<-ps.its.2019 #fungi 2019
ps.f.20<-ps.its.2020 #fungi 2020

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

#bind 2019 and 2020 dataframes 
richtable.f.19 <- richtable.f.19 %>% select(-sample2)
#richtable.f.20 <- richtable.f.20 %>% select(-newroot)
SoilFungi_ASV <- rbind(richtable.f.19, richtable.f.20)

#SAVE CSV
write.csv(SoilFungi_ASV, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilFungi_ASV.csv")




#Genera ----
####Select PS objects for this iteration ####
ps.f.19<-ps.its.2019.gen #fungi 2019
ps.f.20<-ps.its.2020.gen #fungi 2020

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

#bind 2019 and 2020 dataframes 
richtable.f.19 <- richtable.f.19 %>% select(-Row, -Column)
#richtable.f.20 <- richtable.f.20 %>% select()
SoilFungi_Gen<- rbind(richtable.f.19, richtable.f.20)

#SAVE CSV
write.csv(SoilFungi_Gen, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilFungi_Gen.csv")



#Phylum ----

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

#bind 2019 and 2020 dataframes 
richtable.f.19 <- richtable.f.19 %>% select(-sample2)
richtable.f.20 <- richtable.f.20 %>% select(-newroot)
SoilFungi_Phy<- rbind(richtable.f.19, richtable.f.20)

#SAVE CSV
write.csv(SoilFungi_Phy, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/SoilFungi_Phy.csv")

