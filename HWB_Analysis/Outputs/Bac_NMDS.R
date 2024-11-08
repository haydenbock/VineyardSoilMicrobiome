setwd(
  "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs"
)

### Clear workspace ###
rm(list = ls())

## load packages we may need
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c(
  "gridExtra",
  "nlme",
  "lme4",
  "effects",
  "corncob",
  "magrittr",
  "car",
  "interactions",
  "RVAideMemoire",
  "tibble",
  "doParallel",
  "phyloseq",
  "ggplot2",
  "vegan",
  "microbiome",
  "microbiomeSeq",
  "dplyr",
  "MuMIn",
  "Biostrings",
  "ShortRead",
  "ape",
  "ade4",
  "adegraphics",
  "adespatial",
  "tmap",
  "dada2",
  "jtools",
  "rcompanion",
  "ggnewscale",
  "RColorBrewer",
  "emmeans",
  "readxl",
  "pairwiseAdonis"
)

ipak(packages)

remotes::install_github('schuyler-smith/phylosmith')
library(phylosmith)


##load augmented PS and spatial data ##
load(
  "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata"
)


#Data Loading ----
#Bac ASV
    #load phyloseq objects
    ps.b.19 <- ps.16s.2019 #bacteria 2019
    ps.b.20 <- ps.16s.2020 #bacteria 2020
    #Subset soil only compartment
    ps.soil.b.19 <- subset_samples(ps.b.19, type == "soil")
    ps.soil.b.20 <- subset_samples(ps.b.20, type == "soil")
    # Step 1: Extract the OTU table from the phyloseq object
    bac.ASV.2019  <- as.data.frame(otu_table(ps.soil.b.19))
    bac.ASV.2020 <- as.data.frame(otu_table(ps.soil.b.20))
    
#Bac Gen
    #load phyloseq objects
    ps.b.19 <- ps.16s.2019.gen #bacteria 2019
    ps.b.20 <- ps.16s.2020.gen #bacteria 2020
    #Subset soil only compartment
    ps.soil.b.19 <- subset_samples(ps.b.19, type == "soil")
    ps.soil.b.20 <- subset_samples(ps.b.20, type == "soil")
    # Step 1: Extract the OTU table from the phyloseq object
    bac.Gen.2019  <- as.data.frame(otu_table(ps.soil.b.19))
    bac.Gen.2020 <- as.data.frame(otu_table(ps.soil.b.20))
    
#Bac Phy
    #load phyloseq objects
    ps.b.19 <- ps.16s.2019.phy #bacteria 2019
    ps.b.20 <- ps.16s.2020.phy #bacteria 2020
    #Subset soil only compartment
    ps.soil.b.19 <- subset_samples(ps.b.19, type == "soil")
    ps.soil.b.20 <- subset_samples(ps.b.20, type == "soil")
    # Step 1: Extract the OTU table from the phyloseq object
    bac.Phy.2019  <- as.data.frame(otu_table(ps.soil.b.19))
    bac.Phy.2020 <- as.data.frame(otu_table(ps.soil.b.20))

#NMDS ----
    treatments <- ps.16s.2019@sam_data
    treatments <- data.frame(treatments)
    
    nmds_phyloseq(ps.16s.2019, treatment = sample_data(ps.16s.2019)[, c("year", "Row", "Column", "blk", "gc", "rs", "depth", "mass")], 
                  method = 'bray', dimensions = 2, trymax = 100, 
                  circle = 0.95, labels = NULL, colors = 'default', verbose = TRUE)
  
    
    nmds_phyloseq(ps.16s.2019, method = "NMDS", 
                        distance = "bray", 
                        treatment = interaction(sample_data(ps.16s.2019)$gc, 
                                                sample_data(ps.16s.2019)$year))
    

#NMDS doesn't like NA's, doing this for now until we discuss how to handle them.
NMDS_Data <- na.omit(data)



com = NMDS_Data[,14]
com[is.na(com)] <- 0

env <- NMDS_Data %>% select(year:mass)


#convert com to a matrix
m_com = as.matrix(com)




nmds = metaMDS(com, distance = "bray", k = 3, trymax = 1000)
en = envfit(nmds, env, choices=c(1:3), permutations = 999, na.rm = TRUE)

en


plot(nmds)
plot(en)

#PerMANOVA on data
Community <- NMDS_Data[,31:153]

d.manova <- adonis2(Community ~ NMDS_Data$Urban_Kmeans_Cluster*NMDS_Data$Timepoint, method = "euclidean", data= NMDS_Data)
d.manova



adonis.pw_urban <- pairwise.adonis(x = Community,
                                   factors= NMDS_Data$Urban_Kmeans_Cluster,
                                   sim.method='euclidean',
                                   p.adjust.m='holm')

adonis.pw_timepoint <- pairwise.adonis(x = Community,
                                       factors= NMDS_Data$Timepoint,
                                       sim.method='euclidean',
                                       p.adjust.m='holm')


adonis.pw_urban <- as.data.frame(adonis.pw_urban)
adonis.pw_urban$F.Model <- round(adonis.pw_urban$F.Model, 2)
adonis.pw_urban$R2 <- round(adonis.pw_urban$R2, 2)
adonis.pw_urban

