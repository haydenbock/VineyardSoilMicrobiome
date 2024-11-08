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


# Bac.ASV.2019 ====

Bac.ASV.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.ASV.2019.csv")

#clean + process
Bac.ASV.2019 <- Bac.ASV.2019 %>% filter(type == "soil") %>% select(-X, -Row, -Column, -type, -mass)

#make community matrix
Bac.ASV.2019.com = Bac.ASV.2019[,6:16885]
Bac.ASV.2019.com[is.na(Bac.ASV.2019.com)] <- 0
Bac.ASV.2019.com = as.matrix(Bac.ASV.2019.com)

#make treatment matrix
Bac.ASV.2019.env <- Bac.ASV.2019[,1:5]

#run NMDS
Bac.ASV.2019.nmds = metaMDS(Bac.ASV.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.ASV.2019.nmds, Bac.ASV.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.ASV.2019.nmds)
plot(en)



full_model <- adonis2(Bac.ASV.2019.com ~ interaction(Bac.ASV.2019.env$depth, Bac.ASV.2019.env$gc, Bac.ASV.2019.env$rs), method = "bray", data= Bac.ASV.2019)
reduced_model <- adonis2(Bac.ASV.2019.com ~ interaction(Bac.ASV.2019.env$depth, Bac.ASV.2019.env$gc), method = "bray", data= Bac.ASV.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.ASV.2019.com,
                                   factors= interaction(Bac.ASV.2019.env$depth, Bac.ASV.2019.env$gc, Bac.ASV.2019.env$rs),
                                   sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 

Bac.ASV.2019.bray <- vegdist(Bac.ASV.2019.com, method="bray")  
Bac.ASV.2019.betadisper <- betadisper(Bac.ASV.2019.bray, group=interaction(Bac.ASV.2019.env$depth, Bac.ASV.2019.env$gc, Bac.ASV.2019.env$rs))
anova(Bac.ASV.2019.betadisper) #p = 0.005942
TukeyHSD(Bac.ASV.2019.betadisper)
Bac.ASV.2019.betadisper
mean(Bac.ASV.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.ASV.2019.nmds)$sites)
data.scores$depth = Bac.ASV.2019.env$depth
data.scores$gc = Bac.ASV.2019.env$gc
data.scores$rs = Bac.ASV.2019.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)

vif(en_coord_cont)

#plot with ggplot
Bac.asv.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)


# Bac.ASV.2020 ====

Bac.ASV.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.ASV.2020.csv")

#clean + process
Bac.ASV.2020 <- Bac.ASV.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)

#make community matrix
Bac.ASV.2020.com = Bac.ASV.2020[,5:8108]
Bac.ASV.2020.com[is.na(Bac.ASV.2020.com)] <- 0
Bac.ASV.2020.com = as.matrix(Bac.ASV.2020.com)

#make treatment matrix
Bac.ASV.2020.env <- Bac.ASV.2020[,1:4]

#run NMDS
Bac.ASV.2020.nmds = metaMDS(Bac.ASV.2020.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.ASV.2020.nmds, Bac.ASV.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.ASV.2020.nmds)
plot(en)


#PerMANOVA on data
Community <- Bac.ASV.2020[,5:8108]

full_model <- adonis2(Bac.ASV.2020.com ~ interaction(Bac.ASV.2020.env$depth, Bac.ASV.2020.env$gc, Bac.ASV.2020.env$rs), method = "bray", data= Bac.ASV.2020)
reduced_model <- adonis2(Bac.ASV.2020.com ~ interaction(Bac.ASV.2020.env$depth, Bac.ASV.2020.env$gc), method = "bray", data= Bac.ASV.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.ASV.2020.com,
                             factors= interaction(Bac.ASV.2020.env$depth, Bac.ASV.2020.env$gc, Bac.ASV.2020.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Bac.ASV.2020.bray <- vegdist(Bac.ASV.2020.com, method="bray")  
Bac.ASV.2020.betadisper <- betadisper(Bac.ASV.2020.bray, group=interaction(Bac.ASV.2020.env$depth, Bac.ASV.2020.env$gc, Bac.ASV.2020.env$rs))
anova(Bac.ASV.2020.betadisper) #p = 0.7872
TukeyHSD(Bac.ASV.2020.betadisper)
Bac.ASV.2019.betadisper
mean(Bac.ASV.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.ASV.2020.nmds)$sites)
data.scores$depth = Bac.ASV.2020.env$depth
data.scores$gc = Bac.ASV.2020.env$gc
data.scores$rs = Bac.ASV.2020.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)



#plot with ggplot
Bac.asv.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)



# Bac.Gen.2019 ====

Bac.Gen.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Gen.2019.csv")

#clean + process
Bac.Gen.2019 <- Bac.Gen.2019 %>% filter(type == "soil") %>% select(-X, -Row, -Column, -type, -mass)

#make community matrix
Bac.Gen.2019.com = Bac.Gen.2019[,6:610]
Bac.Gen.2019.com[is.na(Bac.Gen.2019.com)] <- 0
Bac.Gen.2019.com = as.matrix(Bac.Gen.2019.com)

#make treatment matrix
Bac.Gen.2019.env <- Bac.Gen.2019[,1:5]

#run NMDS
Bac.Gen.2019.nmds = metaMDS(Bac.Gen.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.Gen.2019.nmds, Bac.Gen.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.Gen.2019.nmds)
plot(en)



full_model <- adonis2(Bac.Gen.2019.com ~ interaction(Bac.Gen.2019.env$depth, Bac.Gen.2019.env$gc, Bac.Gen.2019.env$rs), method = "bray", data= Bac.Gen.2019)
reduced_model <- adonis2(Bac.Gen.2019.com ~ interaction(Bac.Gen.2019.env$depth, Bac.Gen.2019.env$gc), method = "bray", data= Bac.Gen.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.Gen.2019.com,
                             factors= interaction(Bac.Gen.2019.env$depth, Bac.Gen.2019.env$gc, Bac.Gen.2019.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Bac.Gen.2019.bray <- vegdist(Bac.Gen.2019.com, method="bray")  
Bac.Gen.2019.betadisper <- betadisper(Bac.Gen.2019.bray, group=interaction(Bac.Gen.2019.env$depth, Bac.Gen.2019.env$gc, Bac.Gen.2019.env$rs))
anova(Bac.Gen.2019.betadisper) #p = 0.005942
TukeyHSD(Bac.Gen.2019.betadisper)
Bac.Gen.2019.betadisper
mean(Bac.Gen.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.Gen.2019.nmds)$sites) 
data.scores$depth = Bac.Gen.2019.env$depth
data.scores$gc = Bac.Gen.2019.env$gc
data.scores$rs = Bac.Gen.2019.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)

vif(en_coord_cont)

#plot with ggplot
Bac.asv.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(gc)), size = 2, alpha = 0.5)


# Bac.Gen.2020 ====

Bac.Gen.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Gen.2020.csv")

#clean + process
Bac.Gen.2020 <- Bac.Gen.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)  

#make community matrix
Bac.Gen.2020.com = Bac.Gen.2020[,5:342]
Bac.Gen.2020.com[is.na(Bac.Gen.2020.com)] <- 0
Bac.Gen.2020.com = as.matrix(Bac.Gen.2020.com)

#make treatment matrix
Bac.Gen.2020.env <- Bac.Gen.2020[,1:4]

#run NMDS
Bac.Gen.2020.nmds = metaMDS(Bac.Gen.2020.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.Gen.2020.nmds, Bac.Gen.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.Gen.2020.nmds)
plot(en)


#PerMANOVA on data
full_model <- adonis2(Bac.Gen.2020.com ~ interaction(Bac.Gen.2020.env$depth, Bac.Gen.2020.env$gc, Bac.Gen.2020.env$rs), method = "bray", data= Bac.Gen.2020)
reduced_model <- adonis2(Bac.Gen.2020.com ~ interaction(Bac.Gen.2020.env$depth, Bac.Gen.2020.env$gc), method = "bray", data= Bac.Gen.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.Gen.2020.com,
                             factors= interaction(Bac.Gen.2020.env$depth, Bac.Gen.2020.env$gc, Bac.Gen.2020.env$rs),
                             sim.method='bray')



adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Bac.Gen.2020.bray <- vegdist(Bac.Gen.2020.com, method="bray")  
Bac.Gen.2020.betadisper <- betadisper(Bac.Gen.2020.bray, group=interaction(Bac.Gen.2020.env$depth, Bac.Gen.2020.env$gc, Bac.Gen.2020.env$rs))
anova(Bac.Gen.2020.betadisper) #p = 0.7872
TukeyHSD(Bac.Gen.2020.betadisper)
Bac.Gen.2020.betadisper
mean(Bac.Gen.2020.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.Gen.2020.nmds)$sites)
data.scores$depth = Bac.Gen.2020.env$depth
data.scores$gc = Bac.Gen.2020.env$gc
data.scores$rs = Bac.Gen.2020.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)


#plot with ggplot
Bac.Gen.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5)



# Bac.Phy.2019 ====

Bac.Phy.2019 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Phy.2019.csv")

#clean + process
Bac.Phy.2019 <- Bac.Phy.2019 %>% filter(type == "soil") %>% select(-X, -Row, -Column, -type, -mass)

#make community matrix
Bac.Phy.2019.com = Bac.Phy.2019[,6:43]
Bac.Phy.2019.com[is.na(Bac.Phy.2019.com)] <- 0
Bac.Phy.2019.com = as.matrix(Bac.Phy.2019.com)

#make treatment matrix
Bac.Phy.2019.env <- Bac.Phy.2019[,1:5]

#run NMDS
Bac.Phy.2019.nmds = metaMDS(Bac.Phy.2019.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.Phy.2019.nmds, Bac.Phy.2019.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.Phy.2019.nmds)
plot(en)



full_model <- adonis2(Bac.Phy.2019.com ~ interaction(Bac.Phy.2019.env$depth, Bac.Phy.2019.env$gc, Bac.Phy.2019.env$rs), method = "bray", data= Bac.Phy.2019)
reduced_model <- adonis2(Bac.Phy.2019.com ~ interaction(Bac.Phy.2019.env$depth, Bac.Phy.2019.env$gc), method = "bray", data= Bac.Phy.2019)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.Phy.2019.com,
                             factors= interaction(Bac.Phy.2019.env$depth, Bac.Phy.2019.env$gc, Bac.Phy.2019.env$rs),
                             sim.method='bray')

adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Bac.Phy.2019.bray <- vegdist(Bac.Phy.2019.com, method="bray")  
Bac.Phy.2019.betadisper <- betadisper(Bac.Phy.2019.bray, group=interaction(Bac.Phy.2019.env$depth, Bac.Phy.2019.env$gc, Bac.Phy.2019.env$rs))
anova(Bac.Phy.2019.betadisper) #p = 0.9507
TukeyHSD(Bac.Phy.2019.betadisper) 
Bac.Phy.2019.betadisper
mean(Bac.Phy.2019.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.Phy.2019.nmds)$sites) 
data.scores$depth = Bac.Phy.2019.env$depth
data.scores$gc = Bac.Phy.2019.env$gc
data.scores$rs = Bac.Phy.2019.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)


#plot with ggplot
Bac.Phy.2019 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(gc)), size = 2, alpha = 0.5)


# Bac.Phy.2020 ====

Bac.Phy.2020 <- read.csv("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs/NMDS/Bac.Phy.2020.csv")

#clean + process
Bac.Phy.2020 <- Bac.Phy.2020 %>% filter(type == "soil") %>% select(-X, -newroot, -type, -mass, -year)  

#make community matrix
Bac.Phy.2020.com = Bac.Phy.2020[,5:42]
Bac.Phy.2020.com[is.na(Bac.Phy.2020.com)] <- 0
Bac.Phy.2020.com = as.matrix(Bac.Phy.2020.com)

#make treatment matrix
Bac.Phy.2020.env <- Bac.Phy.2020[,1:4]

#run NMDS
Bac.Phy.2020.nmds = metaMDS(Bac.Phy.2020.com, distance = "bray", k = 2, trymax = 1000)
en = envfit(Bac.Phy.2020.nmds, Bac.Phy.2020.env, choices=c(1:3), permutations = 999, na.rm = TRUE)

plot(Bac.Phy.2020.nmds)
plot(en)


#PerMANOVA on data
full_model <- adonis2(Bac.Phy.2020.com ~ interaction(Bac.Phy.2020.env$depth, Bac.Phy.2020.env$gc, Bac.Phy.2020.env$rs), method = "bray", data= Bac.Phy.2020)
reduced_model <- adonis2(Bac.Phy.2020.com ~ interaction(Bac.Phy.2020.env$depth, Bac.Phy.2020.env$gc), method = "bray", data= Bac.Phy.2020)


#pairwise adonis
adonis.pw <- pairwise.adonis(x = Bac.Phy.2020.com,
                             factors= interaction(Bac.Phy.2020.env$depth, Bac.Phy.2020.env$gc, Bac.Phy.2020.env$rs),
                             sim.method='bray')



adonis.pw <- as.data.frame(adonis.pw)
adonis.pw$F.Model <- round(adonis.pw$F.Model, 3)
adonis.pw$R2 <- round(adonis.pw$R2, 2)
adonis.pw

#Beta dispersion 
Bac.Phy.2020.bray <- vegdist(Bac.Phy.2020.com, method="bray")  
Bac.Phy.2020.betadisper <- betadisper(Bac.Phy.2020.bray, group=interaction(Bac.Phy.2020.env$depth, Bac.Phy.2020.env$gc, Bac.Phy.2020.env$rs))
anova(Bac.Phy.2020.betadisper) #p = 0.7872
TukeyHSD(Bac.Phy.2020.betadisper)
Bac.Phy.2020.betadisper
mean(Bac.Phy.2020.bray)


#Make a nice plot in ggplot
#extract coordinates
data.scores = as.data.frame(scores(Bac.Phy.2020.nmds)$sites)
data.scores$depth = Bac.Phy.2020.env$depth
data.scores$gc = Bac.Phy.2020.env$gc
data.scores$rs = Bac.Phy.2020.env$rs

#extract continous and categorical env 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#manually recode variables for color coding
en_coord_cont <- cbind(Variable = rownames(en_coord_cont), en_coord_cont)
rownames(en_coord_cont) <- 1:nrow(en_coord_cont)


#plot with ggplot
Bac.Phy.2020 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(color = as.character(depth), shape = as.character(rs)), size = 2, alpha = 0.5) + 
  facet_wrap(~gc)

