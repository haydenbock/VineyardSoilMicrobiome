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

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)

# Fungi.ASV.2019 ====

Fungi.ASV.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.ASV.2019.csv")

#clean + process
Fungi.ASV.2019 <- Fungi.ASV.2019 %>% filter(type == "soil") %>% select( -Row, -Column, -type, -year, -mass, -sample2, -X)

#make community matrix
Fungi.ASV.2019.com = Fungi.ASV.2019[,6:3075]
Fungi.ASV.2019.com[is.na(Fungi.ASV.2019.com)] <- 0
Fungi.ASV.2019.com = as.matrix(Fungi.ASV.2019.com)

#make treatment matrix
Fungi.ASV.2019.env <- Fungi.ASV.2019[,1:5]

#run NMDS
Fungi.ASV.2019.nmds = metaMDS(Fungi.ASV.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Fungi.ASV.2019.nmds, Fungi.ASV.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.ASV.2019.nmds)
plot(en)



full_model <- adonis2(Fungi.ASV.2019.com ~ interaction(Fungi.ASV.2019.env$depth, Fungi.ASV.2019.env$gc, Fungi.ASV.2019.env$rs), method = "bray", data= Fungi.ASV.2019)
reduced_model <- adonis2(Fungi.ASV.2019.com ~ interaction(Fungi.ASV.2019.env$depth, Fungi.ASV.2019.env$gc), method = "bray", data= Fungi.ASV.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.ASV.2019.com,
                             factors= interaction(Fungi.ASV.2019.env$depth, Fungi.ASV.2019.env$gc, Fungi.ASV.2019.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 

Fungi.ASV.2019.bray <- vegdist(Fungi.ASV.2019.com, method="bray")  
Fungi.ASV.2019.betadisper <- betadisper(Fungi.ASV.2019.bray, group=interaction(Fungi.ASV.2019.env$depth, Fungi.ASV.2019.env$gc, Fungi.ASV.2019.env$rs))
anova(Fungi.ASV.2019.betadisper) #p < 0.001
TukeyHSD(Fungi.ASV.2019.betadisper)
Fungi.ASV.2019.betadisper
mean(Fungi.ASV.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.ASV.2019.nmds)$sites)
data.scores$depth = Fungi.ASV.2019.env$depth
data.scores$gc = Fungi.ASV.2019.env$gc
data.scores$rs = Fungi.ASV.2019.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)

vif(en_coord_cont)

data.scores <- data.scores %>% filter(NMDS1 <= 10000)

#plot with ggplot
Fungi.ASV.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)


# Fungi.ASV.2020 ====

Fungi.ASV.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.ASV.2020.csv")

#clean + process
Fungi.ASV.2020 <- Fungi.ASV.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)

#make community matrix
Fungi.ASV.2020.com = Fungi.ASV.2020[,5:5140]
Fungi.ASV.2020.com[is.na(Fungi.ASV.2020.com)] <- 0
Fungi.ASV.2020.com = as.matrix(Fungi.ASV.2020.com)

#make treatment matrix
Fungi.ASV.2020.env <- Fungi.ASV.2020[,1:4]

#run NMDS
Fungi.ASV.2020.nmds = metaMDS(Fungi.ASV.2020.com, distance = "bray", k = 2 , trymax = 1000) 
en = envfit(Fungi.ASV.2020.nmds, Fungi.ASV.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.ASV.2020.nmds)
plot(en)


#PerMANOVA on data
Community <- Fungi.ASV.2020[,5:5140]

full_model <- adonis2(Fungi.ASV.2020.com ~ interaction(Fungi.ASV.2020.env$depth, Fungi.ASV.2020.env$gc, Fungi.ASV.2020.env$rs), method = "bray", data=Fungi.ASV.2020)
reduced_model <- adonis2(Fungi.ASV.2020.com ~ interaction(Fungi.ASV.2020.env$depth, Fungi.ASV.2020.env$gc), method = "bray", data= Fungi.ASV.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.ASV.2020.com,
                             factors= interaction(Fungi.ASV.2020.env$depth, Fungi.ASV.2020.env$gc, Fungi.ASV.2020.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Fungi.ASV.2020.bray <- vegdist(Fungi.ASV.2020.com, method="bray")  
Fungi.ASV.2020.betadisper <- betadisper(Fungi.ASV.2020.bray, group=interaction(Fungi.ASV.2020.env$depth, Fungi.ASV.2020.env$gc, Fungi.ASV.2020.env$rs))
anova(Fungi.ASV.2020.betadisper) #p = 0.0.6528
TukeyHSD(Fungi.ASV.2020.betadisper)
Fungi.ASV.2019.betadisper
mean(Fungi.ASV.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.ASV.2020.nmds)$sites)
data.scores$depth = Fungi.ASV.2020.env$depth
data.scores$gc = Fungi.ASV.2020.env$gc
data.scores$rs = Fungi.ASV.2020.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)



#plot with ggplot
Fungi.asv.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)



# Fungi.Gen.2019 ====

Fungi.Gen.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Gen.2019.csv")

#clean + process
Fungi.Gen.2019 <- Fungi.Gen.2019 %>% filter(type == "soil") %>% select(-X, -Row, -Column, -type, -mass, -sample2, -year)

#make community matrix
Fungi.Gen.2019.com = Fungi.Gen.2019[,5:319]
Fungi.Gen.2019.com[is.na(Fungi.Gen.2019.com)] <- 0
Fungi.Gen.2019.com = as.matrix(Fungi.Gen.2019.com)

#make treatment matrix
Fungi.Gen.2019.env <- Fungi.Gen.2019[,1:4]

#run NMDS
Fungi.Gen.2019.nmds = metaMDS(Fungi.Gen.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Fungi.Gen.2019.nmds, Fungi.Gen.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.Gen.2019.nmds)
plot(en)



full_model <- adonis2(Fungi.Gen.2019.com ~ interaction(Fungi.Gen.2019.env$depth, Fungi.Gen.2019.env$gc, Fungi.Gen.2019.env$rs), method = "bray", data= Fungi.Gen.2019)
reduced_model <- adonis2(Fungi.Gen.2019.com ~ interaction(Fungi.Gen.2019.env$depth, Fungi.Gen.2019.env$gc), method = "bray", data= Fungi.Gen.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.Gen.2019.com,
                             factors= interaction(Fungi.Gen.2019.env$depth, Fungi.Gen.2019.env$gc, Fungi.Gen.2019.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Fungi.Gen.2019.bray <- vegdist(Fungi.Gen.2019.com, method="bray")  
Fungi.Gen.2019.betadisper <- betadisper(Fungi.Gen.2019.bray, group=interaction(Fungi.Gen.2019.env$depth, Fungi.Gen.2019.env$gc, Fungi.Gen.2019.env$rs))
anova(Fungi.Gen.2019.betadisper) #p = 0.0001
TukeyHSD(Fungi.Gen.2019.betadisper)
Fungi.Gen.2019.betadisper
mean(Fungi.Gen.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.Gen.2019.nmds)$sites) 
data.scores$depth = Fungi.Gen.2019.env$depth
data.scores$gc = Fungi.Gen.2019.env$gc
data.scores$rs = Fungi.Gen.2019.env$rs

data.scores <- data.scores %>% filter(NMDS2 <= 1000)

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)

vif(en_coord_cont)

#plot with ggplot
Fungi.asv.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(gc)), size = 2, alpha = 0.5)


# Fungi.Gen.2020 ====

Fungi.Gen.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Gen.2020.csv")

#clean + process
Fungi.Gen.2020 <- Fungi.Gen.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)  

#make community matrix
Fungi.Gen.2020.com = Fungi.Gen.2020[,5:475]
Fungi.Gen.2020.com[is.na(Fungi.Gen.2020.com)] <- 0
Fungi.Gen.2020.com = as.matrix(Fungi.Gen.2020.com)

#make treatment matrix
Fungi.Gen.2020.env <- Fungi.Gen.2020[,1:4]

#run NMDS
Fungi.Gen.2020.nmds = metaMDS(Fungi.Gen.2020.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Fungi.Gen.2020.nmds, Fungi.Gen.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.Gen.2020.nmds)
plot(en)


#PerMANOVA on data
full_model <- adonis2(Fungi.Gen.2020.com ~ interaction(Fungi.Gen.2020.env$depth, Fungi.Gen.2020.env$gc, Fungi.Gen.2020.env$rs), method = "bray", data= Fungi.Gen.2020)
reduced_model <- adonis2(Fungi.Gen.2020.com ~ interaction(Fungi.Gen.2020.env$depth, Fungi.Gen.2020.env$gc), method = "bray", data= Fungi.Gen.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.Gen.2020.com,
                             factors= interaction(Fungi.Gen.2020.env$depth, Fungi.Gen.2020.env$gc, Fungi.Gen.2020.env$rs),
                             sim.method='bray')



adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Fungi.Gen.2020.bray <- vegdist(Fungi.Gen.2020.com, method="bray")  
Fungi.Gen.2020.betadisper <- betadisper(Fungi.Gen.2020.bray, group=interaction(Fungi.Gen.2020.env$depth, Fungi.Gen.2020.env$gc, Fungi.Gen.2020.env$rs))
anova(Fungi.Gen.2020.betadisper) #p = 0.7538
TukeyHSD(Fungi.Gen.2020.betadisper)
Fungi.Gen.2020.betadisper
mean(Fungi.Gen.2020.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.Gen.2020.nmds)$sites)
data.scores$depth = Fungi.Gen.2020.env$depth
data.scores$gc = Fungi.Gen.2020.env$gc
data.scores$rs = Fungi.Gen.2020.env$rs
data.scores <- data.scores %>% filter(NMDS1 < 1)

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)


#plot with ggplot
Fungi.Gen.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)



# Fungi.Phy.2019 ====

Fungi.Phy.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Phy.2019.csv")

#clean + process
Fungi.Phy.2019 <- Fungi.Phy.2019 %>% filter(type == "soil") %>% select(-X, -Row, -Column, -type, -mass, -sample2, -year)

#make community matrix
Fungi.Phy.2019.com = Fungi.Phy.2019[,5:17]
Fungi.Phy.2019.com[is.na(Fungi.Phy.2019.com)] <- 0
Fungi.Phy.2019.com = as.matrix(Fungi.Phy.2019.com)

#make treatment matrix
Fungi.Phy.2019.env <- Fungi.Phy.2019[,1:4]

#run NMDS
Fungi.Phy.2019.nmds = metaMDS(Fungi.Phy.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Fungi.Phy.2019.nmds, Fungi.Phy.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.Phy.2019.nmds)
plot(en)



full_model <- adonis2(Fungi.Phy.2019.com ~ interaction(Fungi.Phy.2019.env$depth, Fungi.Phy.2019.env$gc, Fungi.Phy.2019.env$rs), method = "bray", data= Fungi.Phy.2019)
reduced_model <- adonis2(Fungi.Phy.2019.com ~ interaction(Fungi.Phy.2019.env$depth, Fungi.Phy.2019.env$gc), method = "bray", data= Fungi.Phy.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.Phy.2019.com,
                             factors= interaction(Fungi.Phy.2019.env$depth, Fungi.Phy.2019.env$gc, Fungi.Phy.2019.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Fungi.Phy.2019.bray <- vegdist(Fungi.Phy.2019.com, method="bray")  
Fungi.Phy.2019.betadisper <- betadisper(Fungi.Phy.2019.bray, group=interaction(Fungi.Phy.2019.env$depth, Fungi.Phy.2019.env$gc, Fungi.Phy.2019.env$rs))
anova(Fungi.Phy.2019.betadisper) #p = 0.01821
TukeyHSD(Fungi.Phy.2019.betadisper) 
Fungi.Phy.2019.betadisper
mean(Fungi.Phy.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.Phy.2019.nmds)$sites) 
data.scores$depth = Fungi.Phy.2019.env$depth
data.scores$gc = Fungi.Phy.2019.env$gc
data.scores$rs = Fungi.Phy.2019.env$rs
data.scores <- data.scores %>% filter(NMDS1 < 800)

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)


#plot with ggplot
Fungi.Phy.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(gc)), size = 2, alpha = 0.5)


# Fungi.Phy.2020 ====

Fungi.Phy.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Fungi.Phy.2020.csv")

#clean + process
Fungi.Phy.2020 <- Fungi.Phy.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)  

#make community matrix
Fungi.Phy.2020.com = Fungi.Phy.2020[,5:20]
Fungi.Phy.2020.com[is.na(Fungi.Phy.2020.com)] <- 0
Fungi.Phy.2020.com = as.matrix(Fungi.Phy.2020.com)

#make treatment matrix
Fungi.Phy.2020.env <- Fungi.Phy.2020[,1:4]

#run NMDS
Fungi.Phy.2020.nmds = metaMDS(Fungi.Phy.2020.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Fungi.Phy.2020.nmds, Fungi.Phy.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Fungi.Phy.2020.nmds)
plot(en)


#PerMANOVA on data
full_model <- adonis2(Fungi.Phy.2020.com ~ interaction(Fungi.Phy.2020.env$depth, Fungi.Phy.2020.env$gc, Fungi.Phy.2020.env$rs), method = "bray", data= Fungi.Phy.2020)
reduced_model <- adonis2(Fungi.Phy.2020.com ~ interaction(Fungi.Phy.2020.env$depth, Fungi.Phy.2020.env$gc), method = "bray", data= Fungi.Phy.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Fungi.Phy.2020.com,
                             factors= interaction(Fungi.Phy.2020.env$depth, Fungi.Phy.2020.env$gc, Fungi.Phy.2020.env$rs),
                             sim.method='bray')



adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Fungi.Phy.2020.bray <- vegdist(Fungi.Phy.2020.com, method="bray")  
Fungi.Phy.2020.betadisper <- betadisper(Fungi.Phy.2020.bray, group=interaction(Fungi.Phy.2020.env$depth, Fungi.Phy.2020.env$gc, Fungi.Phy.2020.env$rs))
anova(Fungi.Phy.2020.betadisper) #p = 0.5689
TukeyHSD(Fungi.Phy.2020.betadisper)
Fungi.Phy.2020.betadisper
mean(Fungi.Phy.2020.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Fungi.Phy.2020.nmds)$sites)
data.scores$depth = Fungi.Phy.2020.env$depth
data.scores$gc = Fungi.Phy.2020.env$gc
data.scores$rs = Fungi.Phy.2020.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)


#plot with ggplot
Fungi.Phy.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)

