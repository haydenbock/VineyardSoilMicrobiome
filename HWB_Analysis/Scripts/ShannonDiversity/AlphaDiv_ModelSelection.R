#Model selection for Bacteria/Fungi at each taxanomic level

Data <- readxl::read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/Summary_Excel.xlsx")

#Remove a weird column: 
Data <- Data %>% select(-"...1")

# Load necessary libraries
library(nlme)
library(MuMIn)
library(emmeans)

# List of Shannon diversity columns
shannon_metrics <- c("Bac_Shannon_ASV", "Bac_Shannon_Gen", "Bac_Shannon_Phy", 
                     "Fungi_Shannon_ASV", "Fungi_Shannon_Gen", "Fungi_Shannon_Phy")

# Loop over each Shannon diversity metric
for (metric in shannon_metrics) {
  
  # Remove rows with missing values in the current Shannon metric
  data_filtered <- Data[!is.na(Data[[metric]]), ]
  
  # Ensure "depth" is treated as a factor
  data_filtered$depth <- as.factor(data_filtered$depth)
  
  # Define the full model
  full_model <- lme(as.formula(paste(metric, "~ depth * gc * rs * year")), 
                    random = ~ 1 | blk/rs, 
                    data = data_filtered)
  
  # Generate all possible combinations of fixed effects and compare models using AIC
  all_models <- dredge(full_model)
  
  # Extract the best models within 2 AIC units of the top model
  best_models_list <- get.models(all_models, subset = delta < 2)
  
  # Assign each best model to a named object dynamically and calculate emmeans
  for (i in seq_along(best_models_list)) {
    model_name <- paste0("best_model_", metric, "_", i)
    assign(model_name, best_models_list[[i]])
    
    # Extract the fixed effects terms from the model
    fixed_effects_terms <- attr(terms(best_models_list[[i]]), "term.labels")
    
    # Dynamically generate the formula for emmeans using only the fixed effects in the model
    emmeans_formula <- as.formula(paste("~", paste(fixed_effects_terms, collapse = " * ")))
    
    # Calculate mean and standard error using emmeans with the selected terms
    emmeans_results <- tryCatch({
      emmeans(best_models_list[[i]], emmeans_formula)
    }, error = function(e) {
      cat("Error in emmeans for", model_name, ":", e$message, "\n")
      NULL
    })
    
    # Print results if emmeans was successful
    if (!is.null(emmeans_results)) {
      cat("\nResults for", model_name, ":\n")
      print(summary(emmeans_results))
    }
  }
}

# Example to access summary results
# emmeans(best_model_Bac_Shannon_ASV_1, ~ depth * gc * rs * year)


#now plot the results: 

#BACTERIA
#Bac ASV
summary(best_model_Bac_Shannon_ASV_1)
Data %>% 
  ggplot(aes(x = as.character(depth), y = Bac_Shannon_ASV, color = as.character(year))) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)

#Bac Gen
summary(best_model_Bac_Shannon_Gen_1)
Data %>% 
  ggplot(aes(x = as.character(year), y = Bac_Shannon_Gen)) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)

#Bac Phyla
summary(best_model_Bac_Shannon_Phy_1)
Data %>% 
  ggplot(aes(x = as.character(depth), y = Bac_Shannon_Phy, color = as.character(year))) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)



#FUNGI
#Fungi ASV
summary(best_model_Fungi_Shannon_ASV_1) 
Data %>% 
  ggplot(aes(x = as.character(year), y = Fungi_Shannon_ASV)) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)

#Fungi Gen
summary(best_model_Fungi_Shannon_Gen_1)
Data %>% 
  ggplot(aes(x = as.character(year), y = Fungi_Shannon_Gen)) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)

#Fungi Phyla
summary(best_model_Fungi_Shannon_Phy_1)
Data %>% 
  ggplot(aes(x = as.character(depth), y = Fungi_Shannon_Phy, color = as.character(year))) +
  geom_boxplot() + 
  geom_point(position = position_dodge(width = 0.75), size = 2, alpha = 0.7)


