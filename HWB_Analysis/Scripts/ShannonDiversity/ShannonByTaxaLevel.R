# Load necessary libraries
library(phyloseq)
library(vegan)  # For diversity calculations (if needed)

load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")



# List of phyloseq objects in your environment
phyloseq_list <- list(Bac19 = ps.16s.2019, Bac20 = ps.16s.2020, Fun19 = ps.its.2019, Fun20 = ps.its.2020)  # Add your phyloseq objects here


# Define the taxonomic levels you want to analyze
taxonomic_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# Define the alpha diversity measures you want to calculate (e.g., Shannon, Simpson, Observed richness)
alpha_measures <- c("Shannon", "Simpson", "Observed")
# Create an empty list to store the results
all_alpha_diversity_results <- list()

# Loop through each phyloseq object
for (ps_name in names(phyloseq_list)) {
  
  # Extract the current phyloseq object
  ps <- phyloseq_list[[ps_name]]
  
  # Create a nested list to store results for this phyloseq object
  alpha_diversity_results <- list()
  
  # Loop through each taxonomic level
  for (tax_level in taxonomic_levels) {
    
    # Aggregate taxa at the current taxonomic level
    ps_tax <- tax_glom(ps, taxrank = tax_level)
    
    # Calculate alpha diversity for all specified measures at this level
    alpha_diversity <- estimate_richness(ps_tax, measures = alpha_measures)
    
    # Store the alpha diversity result for this level
    alpha_diversity_results[[tax_level]] <- alpha_diversity
    
    # Optional: Print a summary of the alpha diversity for this level
    cat("\nAlpha Diversity for", ps_name, "at Taxonomic Level:", tax_level, "\n")
    print(summary(alpha_diversity))
  }
  
  # Store the alpha diversity results for this phyloseq object
  all_alpha_diversity_results[[ps_name]] <- alpha_diversity_results
}

# Now all_alpha_diversity_results contains multiple alpha diversity measures for each phyloseq object and taxonomic level

# Example: Accessing Shannon diversity for ps_object1 at the Class level
shannon_ps1_class <- all_alpha_diversity_results[["ps1"]][["Class"]][,"Shannon"]


simpson_ps2_genus <- all_alpha_diversity_results[["Bac19"]][["Genus"]][,"Simpson"]
observed_ps1_family <- all_alpha_diversity_results[["Bac19"]][["Phylum"]][,"Observed"]









# Inspect the taxonomic table to check for completeness
tax_table(ps.16s.2019)

# Example: Aggregate at Genus and Phylum levels
ps_genus <- tax_glom(ps.16s.2019, taxrank = "Genus")
ps_phylum <- tax_glom(ps.16s.2019, taxrank = "Phylum")

# Check the number of taxa at each level
cat("Number of taxa at Genus level: ", ntaxa(ps_genus), "\n")
cat("Number of taxa at Phylum level: ", ntaxa(ps_phylum), "\n")

# Calculate alpha diversity for each level
alpha_genus <- estimate_richness(ps_genus, measures = c("Shannon", "Simpson"))
alpha_phylum <- estimate_richness(ps_phylum, measures = c("Shannon", "Simpson"))

# Compare the results
cat("Alpha diversity at Genus level:\n")
print(summary(alpha_genus))

cat("\nAlpha diversity at Phylum level:\n")
print(summary(alpha_phylum))
