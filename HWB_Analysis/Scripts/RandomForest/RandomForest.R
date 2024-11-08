
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
packages <- c("gridExtra","randomForest","nlme","lme4", "effects","corncob","magrittr","car", "interactions","RVAideMemoire","tibble","doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq",
              "dplyr","MuMIn","readxl","Biostrings","ShortRead","ape","ade4","adegraphics","adespatial","tmap","dada2","jtools","rcompanion","ggnewscale","RColorBrewer","emmeans")
ipak(packages)


##If you happen to need original phyloseq objects
#load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")

RF_Data <- read_excel("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/Summary_Excel.xlsx")

# get rid of weird extra column, and make sure that we are only looking at soil data
RF_Data <- RF_Data %>% select(-`...1`) %>% filter(type == "soil") %>% select(-type)

#Random Forest Approach for Abundance -----
# Make 3 RF's - one for each depth
#Lets start with random seed so the outcome will be repeatable and store train and test data.
set.seed(222)

#rename to readable/simple variables.
TempRF <- RF_Data %>% select(year:mass, Bac_Shannon_Phy, Fungi_Shannon_Phy, Bac_Bray_Curtis_Phy, Fungi_Bray_Curtis_Phy)

#make readable names
TempRF <- TempRF %>% dplyr::rename(Year_of_study = `year`,
         Experimental_block = blk,
         Groundcover = gc,
         Rootstock = rs, 
         Depth = depth,
         Root_mass = mass,
         Bacteria_Alpha = Bac_Shannon_Phy,
         Fungi_Alpha = Fungi_Shannon_Phy,
         Bacteria_Beta = Bac_Bray_Curtis_Phy,
         Fungi_Beta = Fungi_Bray_Curtis_Phy)
         
# Bac Alpha Div=====
#narrow down to the one variable of interest
Bac_Alpha_RF <- TempRF %>% select(-Fungi_Alpha, -Fungi_Beta, -Bacteria_Beta) 

ind <- sample(2, nrow(Bac_Alpha_RF), replace = TRUE, prob = c(0.7, 0.3))
train <- Bac_Alpha_RF[ind==1,]
test <- Bac_Alpha_RF[ind==2,]

train <- as.data.frame(train)


train_BacAlpha <- train 
test_BacAlpha<- test 

#UNCOMMENT IF YOU WANT TO ANALYZE ACROSS DEPTH
#train_shallow <- train %>% filter(Depth == "0") %>% select(-Depth)
#test_shallow<- test %>% filter(Depth == "0") %>% select(-Depth)

#train_moderate <- train %>% filter(Depth == "33") %>% select(-Depth)
#test_moderate<- test %>% filter(Depth == "33") %>% select(-Depth)

#train_deep <- train %>% filter(Depth == "66")%>% select(-Depth)
#test_deep<- test %>% filter(Depth == "66")%>% select(-Depth)


#Random Forest in R
BacAlphaRF<- randomForest(Bacteria_Alpha~., data=train_BacAlpha, ntree=5000, na.action=na.omit, importance=TRUE) 

print(BacAlphaRF)

#extra packages
install.packages("~/Downloads/rfUtilities_2.1-5.tar.gz", repos = NULL, type="source")
library(rfUtilities)

#Extract MSEincrease%
BacAlphaRF_Imp<-importance(BacAlphaRF)  #defaults scale=TRUE, standardizes  
BacAlphaRF_Imp<-data.frame(BacAlphaRF_Imp) #convert to df.




#Fungi Alpha Div ----

#narrow down to the one variable of interest
Fungi_Alpha_RF <- TempRF %>% select(-Bacteria_Alpha, -Fungi_Beta, -Bacteria_Beta) 

ind <- sample(2, nrow(Fungi_Alpha_RF), replace = TRUE, prob = c(0.7, 0.3))
train <- Fungi_Alpha_RF[ind==1,]
test <- Fungi_Alpha_RF[ind==2,]

train <- as.data.frame(train)


train_FungiAlpha <- train 
test_FungiAlpha<- test 

#UNCOMMENT IF YOU WANT TO ANALYZE ACROSS DEPTH
#train_shallow <- train %>% filter(Depth == "0") %>% select(-Depth)
#test_shallow<- test %>% filter(Depth == "0") %>% select(-Depth)

#train_moderate <- train %>% filter(Depth == "33") %>% select(-Depth)
#test_moderate<- test %>% filter(Depth == "33") %>% select(-Depth)

#train_deep <- train %>% filter(Depth == "66")%>% select(-Depth)
#test_deep<- test %>% filter(Depth == "66")%>% select(-Depth)


#Random Forest in R
FungiAlphaRF<- randomForest(Fungi_Alpha~., data=train_FungiAlpha, ntree=5000, na.action=na.omit, importance=TRUE) 

print(FungiAlphaRF)

#extra packages
install.packages("~/Downloads/rfUtilities_2.1-5.tar.gz", repos = NULL, type="source")
library(rfUtilities)

#Extract MSEincrease%
FungiAlphaRF_Imp<-importance(FungiAlphaRF)  #defaults scale=TRUE, standardizes  
FungiAlphaRF_Imp<-data.frame(FungiAlphaRF_Imp) #convert to df.

# Bac Beta Div=====
#narrow down to the one variable of interest
Bac_Beta_RF <- TempRF %>% select(-Fungi_Alpha, -Fungi_Beta, -Bacteria_Alpha) 

ind <- sample(2, nrow(Bac_Beta_RF), replace = TRUE, prob = c(0.7, 0.3))
train <- Bac_Beta_RF[ind==1,]
test <- Bac_Beta_RF[ind==2,]

train <- as.data.frame(train)


train_BacBeta <- train 
test_BacBeta<- test 

#UNCOMMENT IF YOU WANT TO ANALYZE ACROSS DEPTH
#train_shallow <- train %>% filter(Depth == "0") %>% select(-Depth)
#test_shallow<- test %>% filter(Depth == "0") %>% select(-Depth)

#train_moderate <- train %>% filter(Depth == "33") %>% select(-Depth)
#test_moderate<- test %>% filter(Depth == "33") %>% select(-Depth)

#train_deep <- train %>% filter(Depth == "66")%>% select(-Depth)
#test_deep<- test %>% filter(Depth == "66")%>% select(-Depth)


#Random Forest in R
BacBetaRF<- randomForest(Bacteria_Beta~., data=train_BacBeta, ntree=5000, na.action=na.omit, importance=TRUE) 

print(BacBetaRF)

#extra packages
install.packages("~/Downloads/rfUtilities_2.1-5.tar.gz", repos = NULL, type="source")
library(rfUtilities)

#Extract MSEincrease%
BacBetaRF_Imp<-importance(BacBetaRF)  #defaults scale=TRUE, standardizes  
BacBetaRF_Imp<-data.frame(BacBetaRF_Imp) #convert to df.




#Fungi Beta Div ----

#narrow down to the one variable of interest
Fungi_Beta_RF <- TempRF %>% select(-Bacteria_Alpha, -Fungi_Alpha, -Bacteria_Beta) 

ind <- sample(2, nrow(Fungi_Beta_RF), replace = TRUE, prob = c(0.7, 0.3))
train <- Fungi_Beta_RF[ind==1,]
test <- Fungi_Beta_RF[ind==2,]

train <- as.data.frame(train)


train_FungiBeta <- train 
test_FungiBeta<- test 

#UNCOMMENT IF YOU WANT TO ANALYZE ACROSS DEPTH
#train_shallow <- train %>% filter(Depth == "0") %>% select(-Depth)
#test_shallow<- test %>% filter(Depth == "0") %>% select(-Depth)

#train_moderate <- train %>% filter(Depth == "33") %>% select(-Depth)
#test_moderate<- test %>% filter(Depth == "33") %>% select(-Depth)

#train_deep <- train %>% filter(Depth == "66")%>% select(-Depth)
#test_deep<- test %>% filter(Depth == "66")%>% select(-Depth)


#Random Forest in R
FungiBetaRF<- randomForest(Fungi_Beta~., data=train_FungiBeta, ntree=5000, na.action=na.omit, importance=TRUE) 

print(FungiBetaRF)


#Extract MSEincrease%
FungiBetaRF_Imp<-importance(FungiBetaRF)  #defaults scale=TRUE, standardizes  
FungiBetaRF_Imp<-data.frame(FungiBetaRF_Imp) #convert to df.
