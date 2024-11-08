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

#Bac ASV level----
# Extract OTU and Sample tables from the ps.16s.2019 object
otu_table <- as.data.frame(otu_table(ps.16s.2019))
sample_data <- as.data.frame(sample_data(ps.16s.2019))

colnames(sample_data)

# Specify the columns from sample_data that contain treatment information
treatment_columns <- c("year", "Row", "Column", "type", "blk", "gc", "rs", "depth", "mass") 

# Check if the OTU table row names match the sample data row names
if (!identical(rownames(otu_table), rownames(sample_data))) {
  stop("Row names of OTU and sample data do not match.")
}

# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data[, treatment_columns], otu_table)


# View the combined data
head(combined_data)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.16s.2020)) 
sample_data <- as.data.frame(sample_data(ps.16s.2020))

colnames(sample_data)

# Specify the columns from sample_data that contain treatment information
treatment_columns <- c("year", "Row", "Column", "type", "blk", "gc", "rs", "depth", "mass") 

# Check if the OTU table row names match the sample data row names
if (!identical(rownames(otu_table), rownames(sample_data))) {
  stop("Row names of OTU and sample data do not match.")
}

# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Bac.ASV.2019 <- combined_data.2019
Bac.ASV.2020 <- combined_data.2020

write.csv(Bac.ASV.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/Bac.ASV.2019.csv")
write.csv(Bac.ASV.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/Bac.ASV.2020.csv")







#Bac Gen level----
# Extract OTU and Sample tables from the ps.16s.2019 object
otu_table <- as.data.frame(otu_table(ps.16s.2019.gen))
sample_data <- as.data.frame(sample_data(ps.16s.2019.gen))



# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data, otu_table)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.16s.2020.gen)) 
sample_data <- as.data.frame(sample_data(ps.16s.2020.gen))


# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Bac.Gen.2019 <- combined_data.2019
Bac.Gen.2020 <- combined_data.2020

write.csv(Bac.Gen.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Gen.2019.csv")
write.csv(Bac.Gen.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Gen.2020.csv")







#Bac Phy level----
# Extract OTU and Sample tables from the ps.16s.2019 object
otu_table <- as.data.frame(otu_table(ps.16s.2019.phy))
sample_data <- as.data.frame(sample_data(ps.16s.2019.phy))



# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data, otu_table)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.16s.2020.phy)) 
sample_data <- as.data.frame(sample_data(ps.16s.2020.phy))


# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Bac.Phy.2019 <- combined_data.2019
Bac.Phy.2020 <- combined_data.2020

write.csv(Bac.Phy.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Phy.2019.csv")
write.csv(Bac.Phy.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Phy.2020.csv")








#Fungi ASV level----
# Extract OTU and Sample tables from the ITS object
otu_table <- as.data.frame(otu_table(ps.its.2019))
sample_data <- as.data.frame(sample_data(ps.its.2019))

colnames(sample_data)


# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data, otu_table)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.its.2020))
sample_data <- as.data.frame(sample_data(ps.its.2020))


# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Fungi.ASV.2019 <- combined_data.2019
Fungi.ASV.2020 <- combined_data.2020

write.csv(Fungi.ASV.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.ASV.2019.csv")
write.csv(Fungi.ASV.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.ASV.2020.csv")







#Fungi Gen level----
# Extract OTU and Sample tables from the ITS object
otu_table <- as.data.frame(otu_table(ps.its.2019.gen))
sample_data <- as.data.frame(sample_data(ps.its.2019.gen))

colnames(sample_data)


# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data, otu_table)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.its.2020.gen))
sample_data <- as.data.frame(sample_data(ps.its.2020.gen))


# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Fungi.Gen.2019 <- combined_data.2019
Fungi.Gen.2020 <- combined_data.2020

write.csv(Fungi.Gen.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Gen.2019.csv")
write.csv(Fungi.Gen.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Gen.2020.csv")







#Fungi Gen level----
# Extract OTU and Sample tables from the ITS object
otu_table <- as.data.frame(otu_table(ps.its.2019.phy))
sample_data <- as.data.frame(sample_data(ps.its.2019.phy))

colnames(sample_data)


# Merge the OTU table with the selected treatment columns
combined_data.2019 <- cbind(sample_data, otu_table)


# Extract OTU and Sample tables from the ps.16s.2020 object
otu_table <- as.data.frame(otu_table(ps.its.2020.phy))
sample_data <- as.data.frame(sample_data(ps.its.2020.phy))


# Merge the OTU table with the selected treatment columns
combined_data.2020 <- cbind(sample_data, otu_table)


Fungi.Phy.2019 <- combined_data.2019
Fungi.Phy.2020 <- combined_data.2020

write.csv(Fungi.Phy.2019, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Phy.2019.csv")
write.csv(Fungi.Phy.2020, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Phy.2020.csv")






