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
  "emmeans"
)
ipak(packages)


##load augmented PS and spatial data ##
load(
  "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata"
)

#ASV LEVEL====
#load phyloseq objects
ps.f.19 <- ps.its.2019 #fungi 2019
ps.f.20 <- ps.its.2020 #fungi 2020


#Subset soil only compartment
ps.soil.f.19 <- subset_samples(ps.f.19, type == "soil")
ps.soil.f.20 <- subset_samples(ps.f.20, type == "soil")

# Step 1: Extract the OTU table from the phyloseq object
otu_table_data_2019 <- as.data.frame(otu_table(ps.soil.f.19))
otu_table_data_2020 <- as.data.frame(otu_table(ps.soil.f.20))


# Step 2: Transform data as needed
# For Bray-Curtis, Jaccard, and Euclidean, use raw abundance data.
# For Raup-Crick, we need presence-absence data.
presence_absence_data_2019 <- otu_table_data_2019 > 0  # Convert to binary presence/absence
presence_absence_data_2020 <- otu_table_data_2020 > 0  # Convert to binary presence/absence

# Step 3: Calculate Beta Diversity Metrics
# Bray-Curtis Dissimilarity
bray_curtis_2019 <- vegdist(otu_table_data_2019, method = "bray")
bray_curtis_2020 <- vegdist(otu_table_data_2020, method = "bray")


# Jaccard Index
jaccard_2019 <- vegdist(otu_table_data_2019, method = "jaccard")
jaccard_2020 <- vegdist(otu_table_data_2020, method = "jaccard")


#Euclidean Distance
euclidean_2019 <- vegdist(otu_table_data_2019, method = "euclidean")
euclidean_2020 <- vegdist(otu_table_data_2020, method = "euclidean")


# Raup-Crick Dissimilarity (REMINDER - requires presence-absence data)
raup_crick_2019 <- raupcrick(presence_absence_data_2019)
raup_crick_2020 <- raupcrick(presence_absence_data_2020)


# Convert each distance metric to a tibble
bray_curtis_2019 <- as_tibble(bray_curtis_2019)
jaccard_2019 <- as_tibble(jaccard_2019)
raup_crick_2019 <- as_tibble(raup_crick_2019)
euclidean_2019 <- as_tibble(euclidean_2019)

bray_curtis_2020 <- as_tibble(bray_curtis_2020)
jaccard_2020 <- as_tibble(jaccard_2020)
raup_crick_2020 <- as_tibble(raup_crick_2020)
euclidean_2020 <- as_tibble(euclidean_2020)


# Combine them into a data frame
beta_diversity_metrics_2019 <- data.frame(
  bray_curtis_2019 = bray_curtis_2019,
  jaccard_2019 = jaccard_2019,
  raup_crick_2019 = raup_crick_2019,
  euclidean_2019 = euclidean_2019
)


beta_diversity_metrics_2020 <- data.frame(
  bray_curtis_2020 = bray_curtis_2020,
  jaccard_2020 = jaccard_2020,
  raup_crick_2020 = raup_crick_2020,
  euclidean_2020 = euclidean_2020
)


beta_diversity_metrics_2019 <- as.data.frame(beta_diversity_metrics_2019)
beta_diversity_metrics_2019 <- beta_diversity_metrics_2019[1:24, ] #24 corresponds to sample num in 2019 OTU table 


beta_diversity_metrics_2020 <- as.data.frame(beta_diversity_metrics_2020)
beta_diversity_metrics_2020 <- beta_diversity_metrics_2020[1:48, ] #this is to remove the repeating rows that occur for some reason in the calcs. 47 corresponds to sample num in 2020 OTU table







#GENERA LEVEL====
#load phyloseq objects
ps.f.19 <- ps.its.2019.gen #fungi by genera 2019
ps.f.20 <- ps.its.2020.gen #fungi by genera 2020


#Subset soil only compartment
ps.soil.f.19 <- subset_samples(ps.f.19, type == "soil")
ps.soil.f.20 <- subset_samples(ps.f.20, type == "soil")

# Step 1: Extract the OTU table from the phyloseq object
otu_table_data_2019 <- as.data.frame(otu_table(ps.soil.f.19))
otu_table_data_2020 <- as.data.frame(otu_table(ps.soil.f.20))


# Step 2: Transform data as needed
# For Bray-Curtis, Jaccard, and Euclidean, use raw abundance data.
# For Raup-Crick, we need presence-absence data.
presence_absence_data_2019 <- otu_table_data_2019 > 0  # Convert to binary presence/absence
presence_absence_data_2020 <- otu_table_data_2020 > 0  # Convert to binary presence/absence

# Step 3: Calculate Beta Diversity Metrics
# Bray-Curtis Dissimilarity
bray_curtis_2019 <- vegdist(otu_table_data_2019, method = "bray")
bray_curtis_2020 <- vegdist(otu_table_data_2020, method = "bray")


# Jaccard Index
jaccard_2019 <- vegdist(otu_table_data_2019, method = "jaccard")
jaccard_2020 <- vegdist(otu_table_data_2020, method = "jaccard")


#Euclidean Distance
euclidean_2019 <- vegdist(otu_table_data_2019, method = "euclidean")
euclidean_2020 <- vegdist(otu_table_data_2020, method = "euclidean")


# Raup-Crick Dissimilarity (REMINDER - requires presence-absence data)
raup_crick_2019 <- raupcrick(presence_absence_data_2019)
raup_crick_2020 <- raupcrick(presence_absence_data_2020)




# Convert each distance metric to a tibble
bray_curtis_2019 <- as_tibble(bray_curtis_2019)
jaccard_2019 <- as_tibble(jaccard_2019)
raup_crick_2019 <- as_tibble(raup_crick_2019)
euclidean_2019 <- as_tibble(euclidean_2019)

bray_curtis_2020 <- as_tibble(bray_curtis_2020)
jaccard_2020 <- as_tibble(jaccard_2020)
raup_crick_2020 <- as_tibble(raup_crick_2020)
euclidean_2020 <- as_tibble(euclidean_2020)


# Combine them into a data frame
beta_diversity_metrics_2019 <- data.frame(
  bray_curtis_2019 = bray_curtis_2019,
  jaccard_2019 = jaccard_2019,
  raup_crick_2019 = raup_crick_2019,
  euclidean_2019 = euclidean_2019
)


beta_diversity_metrics_2020 <- data.frame(
  bray_curtis_2020 = bray_curtis_2020,
  jaccard_2020 = jaccard_2020,
  raup_crick_2020 = raup_crick_2020,
  euclidean_2020 = euclidean_2020
)


beta_diversity_metrics_2019 <- as.data.frame(beta_diversity_metrics_2019)
beta_diversity_metrics_2019 <- beta_diversity_metrics_2019[1:24, ] #24 corresponds to sample num in 2019 OTU table 


beta_diversity_metrics_2020 <- as.data.frame(beta_diversity_metrics_2020)
beta_diversity_metrics_2020 <- beta_diversity_metrics_2020[1:48, ] #this is to remove the repeating rows that occur for some reason in the calcs. 47 corresponds to sample num in 2020 OTU table

#PHYLUM LEVEL====
#load phyloseq objects
ps.f.19 <- ps.its.2019.phy #fungi by phylum 2019
ps.f.20 <- ps.its.2020.phy #fungi by phylum 2020


#Subset soil only compartment
ps.soil.f.19 <- subset_samples(ps.f.19, type == "soil")
ps.soil.f.20 <- subset_samples(ps.f.20, type == "soil")

# Step 1: Extract the OTU table from the phyloseq object
otu_table_data_2019 <- as.data.frame(otu_table(ps.soil.f.19))
otu_table_data_2020 <- as.data.frame(otu_table(ps.soil.f.20))


# Step 2: Transform data as needed
# For Bray-Curtis, Jaccard, and Euclidean, use raw abundance data.
# For Raup-Crick, we need presence-absence data.
presence_absence_data_2019 <- otu_table_data_2019 > 0  # Convert to binary presence/absence
presence_absence_data_2020 <- otu_table_data_2020 > 0  # Convert to binary presence/absence

# Step 3: Calculate Beta Diversity Metrics
# Bray-Curtis Dissimilarity
bray_curtis_2019 <- vegdist(otu_table_data_2019, method = "bray")
bray_curtis_2020 <- vegdist(otu_table_data_2020, method = "bray")


# Jaccard Index
jaccard_2019 <- vegdist(otu_table_data_2019, method = "jaccard")
jaccard_2020 <- vegdist(otu_table_data_2020, method = "jaccard")


#Euclidean Distance
euclidean_2019 <- vegdist(otu_table_data_2019, method = "euclidean")
euclidean_2020 <- vegdist(otu_table_data_2020, method = "euclidean")


# Raup-Crick Dissimilarity (REMINDER - requires presence-absence data)
raup_crick_2019 <- raupcrick(presence_absence_data_2019)
raup_crick_2020 <- raupcrick(presence_absence_data_2020)



# Convert each distance metric to a tibble
bray_curtis_2019 <- as_tibble(bray_curtis_2019)
jaccard_2019 <- as_tibble(jaccard_2019)
raup_crick_2019 <- as_tibble(raup_crick_2019)
euclidean_2019 <- as_tibble(euclidean_2019)

bray_curtis_2020 <- as_tibble(bray_curtis_2020)
jaccard_2020 <- as_tibble(jaccard_2020)
raup_crick_2020 <- as_tibble(raup_crick_2020)
euclidean_2020 <- as_tibble(euclidean_2020)


# Combine them into a data frame
beta_diversity_metrics_2019 <- data.frame(
  bray_curtis_2019 = bray_curtis_2019,
  jaccard_2019 = jaccard_2019,
  raup_crick_2019 = raup_crick_2019,
  euclidean_2019 = euclidean_2019
)


beta_diversity_metrics_2020 <- data.frame(
  bray_curtis_2020 = bray_curtis_2020,
  jaccard_2020 = jaccard_2020,
  raup_crick_2020 = raup_crick_2020,
  euclidean_2020 = euclidean_2020
)


beta_diversity_metrics_2019 <- as.data.frame(beta_diversity_metrics_2019)
beta_diversity_metrics_2019 <- beta_diversity_metrics_2019[1:24, ] #24 corresponds to sample num in 2019 OTU table 


beta_diversity_metrics_2020 <- as.data.frame(beta_diversity_metrics_2020)
beta_diversity_metrics_2020 <- beta_diversity_metrics_2020[1:48, ] #this is to remove the repeating rows that occur for some reason in the calcs. 47 corresponds to sample num in 2020 OTU table



