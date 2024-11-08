#### Setup ####

setwd(
  "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/HaydenOutputs"
)


### Clear workspace ###
rm(list=ls())

### Packages ###
#Corncob Tutorial and more info: https://rdrr.io/github/bryandmartin/corncob/f/vignettes/corncob-intro.Rmd
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
#devtools::install_github("bryandmartin/corncob")
### Load Packages ###
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("corncob","phyloseq", "ggplot2","DESeq2","microbiomeSeq","metagMisc","dplyr","fastDummies","stringr")
ipak(packages)

### Prep Phyloseq Objects ###
load("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/HWB_Analysis/Data&Inputs/Inputs/PS_ccxrs.Rdata")


####CORNCOB Bacteria 2019 Phylum ====

ps<-ps.16s.2019.phy

#select just the soil compartment
ps.sub = ps_b <-subset_samples(ps, type=="soil" )
ps.sub@sam_data$depth<-as.factor(ps.sub@sam_data$depth)



set.seed(1)
da_rootchar <- differentialTest(formula = ~ gc+rs+depth+blk, #the DA formula - if you have a controlling factor it must be included. 
                               phi.formula = ~ gc+rs+depth+blk, #the DV formula - if you have a controlling factor it must be included
                               formula_null = ~ blk, # DA controlling factors
                               phi.formula_null = ~blk,# HB: I changed this code to ~blk, in accordance with https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html; from my understanding, the previous formula would not actually show null dispersion.
                               test = "Wald", boot = FALSE, #wald test is "standard"
                               data = ps.sub, #PS object
                               fdr_cutoff = 0.05) #pval false discovery value
plot(da_rootchar, total = TRUE, B = 50)

#explore significant taxa
l<-as.data.frame(da_rootchar$significant_taxa)
models<-da_rootchar$significant_models
l


#extract means/SE and create new table
p<-plot(da_rootchar)
d<-as.data.frame(p$data)

da_rootchar$p
da_rootchar$significant_models










#
#
# Bacteria 2019
# Each Depth Individually 
#
#
#


##prep phyloseq/subset by depth

ps.sub.top = subset_samples(ps.sub, depth=="0" )
ps.sub.mid = subset_samples(ps.sub, depth=="33" )
ps.sub.bot = subset_samples(ps.sub, depth=="66" )


#create reference table of taxa to ASV
all.tax<-as.data.frame(ps.sub@tax_table@.Data)
phy.gen.all<-cbind(rownames(all.tax), all.tax$sequence, all.tax$Phylum,all.tax$Family,all.tax$Genus)
phy.gen<-as.data.frame(unique(phy.gen.all))
colnames(phy.gen)<-c("asv","seq", "Phylum", "Family","Genus")

##DA Top
###NOTE model overspecified with blk*rs

set.seed(1)
da_top<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.top, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_top$significant_taxa)

#extract means/SE and create new table
p.top<-plot(da_top)
d.top<-as.data.frame(p.top$data)


##DA Mid
###NOTE model overspecified with blk*rs

set.seed(1)
da_mid<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.mid, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_mid$significant_taxa)

#extract means/SE and create new table
p.mid<-plot(da_mid)
d.mid<-as.data.frame(p.mid$data)


##DA bot
###NOTE model overspecified with blk*rs

set.seed(1)
da_bot<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.bot, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_bot$significant_taxa)

#extract means/SE and create new table
p.bot<-plot(da_bot)
d.bot<-as.data.frame(p.bot$data)

d.top
d.mid
d.bot









####CORNCOB Bacteria 2020 Phylum ====

ps<-ps.16s.2020.phy

#select just the soil compartment
ps.sub = ps_b <-subset_samples(ps, type=="soil" )
ps.sub@sam_data$depth<-as.factor(ps.sub@sam_data$depth)



set.seed(1)
da_rootchar <- differentialTest(formula = ~ gc+rs+depth+blk, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ gc+rs+depth+blk, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~ blk, # DA controlling factors
                                phi.formula_null = ~blk,# HB: I changed this code to ~blk, in accordance with https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html; from my understanding, the previous formula would not actually show null dispersion.
                                test = "Wald", boot = FALSE, #wald test is "standard"
                                data = ps.sub, #PS object
                                fdr_cutoff = 0.05) #pval false discovery value
plot(da_rootchar, total = TRUE, B = 50)

#explore significant taxa
l<-as.data.frame(da_rootchar$significant_taxa)
models<-da_rootchar$significant_models
l

#extract means/SE and create new table
p<-plot(da_rootchar)
d<-as.data.frame(p$data)

da_rootchar$p
da_rootchar$significant_models










#
#
# Bacteria 2020
# Each Depth Individually 
#
#
#


##prep phyloseq/subset by depth

ps.sub.top = subset_samples(ps.sub, depth=="0" )
ps.sub.mid = subset_samples(ps.sub, depth=="33" )
ps.sub.bot = subset_samples(ps.sub, depth=="66" )


#create reference table of taxa to ASV
all.tax<-as.data.frame(ps.sub@tax_table@.Data)
phy.gen.all<-cbind(rownames(all.tax), all.tax$sequence, all.tax$Phylum,all.tax$Family,all.tax$Genus)
phy.gen<-as.data.frame(unique(phy.gen.all))
colnames(phy.gen)<-c("asv","seq", "Phylum", "Family","Genus")

##DA Top
###NOTE model overspecified with blk*rs

set.seed(1)
da_top<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.top, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_top$significant_taxa)

#extract means/SE and create new table
p.top<-plot(da_top)
d.top<-as.data.frame(p.top$data)


##DA Mid
###NOTE model overspecified with blk*rs

set.seed(1)
da_mid<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.mid, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_mid$significant_taxa)

#extract means/SE and create new table
p.mid<-plot(da_mid)
d.mid<-as.data.frame(p.mid$data)


##DA bot
###NOTE model overspecified with blk*rs

set.seed(1)
da_bot<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.bot, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_bot$significant_taxa)

#extract means/SE and create new table
p.bot<-plot(da_bot)
d.bot<-as.data.frame(p.bot$data)

d.top
d.mid
d.bot














####CORNCOB Fungi 2019 Phylum ====

ps<-ps.its.2019.phy

#select just the soil compartment
ps.sub = ps_b <-subset_samples(ps, type=="soil" )
ps.sub@sam_data$depth<-as.factor(ps.sub@sam_data$depth)



set.seed(1)
da_rootchar <- differentialTest(formula = ~ gc+rs+depth+blk, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ gc+rs+depth+blk, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~ blk, # DA controlling factors
                                phi.formula_null = ~blk,# HB: I changed this code to ~blk, in accordance with https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html; from my understanding, the previous formula would not actually show null dispersion.
                                test = "Wald", boot = FALSE, #wald test is "standard"
                                data = ps.sub, #PS object
                                fdr_cutoff = 0.05) #pval false discovery value
plot(da_rootchar, total = TRUE, B = 50)

#explore significant taxa
l<-as.data.frame(da_rootchar$significant_taxa)
models<-da_rootchar$significant_models
l

#extract means/SE and create new table
p<-plot(da_rootchar) 
d<-as.data.frame(p$data)

da_rootchar$p
da_rootchar$significant_models










#
#
# Fungi 2019
# Each Depth Individually 
#
#
#


##prep phyloseq/subset by depth

ps.sub.top = subset_samples(ps.sub, depth=="0" )
ps.sub.mid = subset_samples(ps.sub, depth=="33" )
ps.sub.bot = subset_samples(ps.sub, depth=="66" )


#create reference table of taxa to ASV
all.tax<-as.data.frame(ps.sub@tax_table@.Data)
phy.gen.all<-cbind(rownames(all.tax), all.tax$sequence, all.tax$Phylum,all.tax$Family,all.tax$Genus)
phy.gen<-as.data.frame(unique(phy.gen.all))
colnames(phy.gen)<-c("asv","seq", "Phylum", "Family","Genus")

##DA Top
###NOTE model overspecified with rs included

set.seed(1)
da_top<- differentialTest(formula = ~ gc+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~ blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.top, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_top$significant_taxa)

#extract means/SE and create new table
p.top<-plot(da_top)
d.top<-as.data.frame(p.top$data)


##DA Mid
###NOTE model overspecified with blk*rs

set.seed(1)
da_mid<- differentialTest(formula = ~ blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.mid, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_mid$significant_taxa)

#extract means/SE and create new table
p.mid<-plot(da_mid)
d.mid<-as.data.frame(p.mid$data)


##DA bot
###NOTE model overspecified with blk*rs

set.seed(1)
da_bot<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.bot, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_bot$significant_taxa)

#extract means/SE and create new table
p.bot<-plot(da_bot)
d.bot<-as.data.frame(p.bot$data)

d.top
d.mid
d.bot









####CORNCOB Fungi 2020 Phylum ====

ps<-ps.its.2020.phy

#select just the soil compartment
ps.sub = ps_b <-subset_samples(ps, type=="soil" )
ps.sub@sam_data$depth<-as.factor(ps.sub@sam_data$depth)



set.seed(1)
da_rootchar <- differentialTest(formula = ~ gc+rs+depth+blk, #the DA formula - if you have a controlling factor it must be included. 
                                phi.formula = ~ gc+rs+depth+blk, #the DV formula - if you have a controlling factor it must be included
                                formula_null = ~ blk, # DA controlling factors
                                phi.formula_null = ~blk,# HB: I changed this code to ~blk, in accordance with https://rdrr.io/github/bryandmartin/corncob/man/differentialTest.html; from my understanding, the previous formula would not actually show null dispersion.
                                test = "Wald", boot = FALSE, #wald test is "standard"
                                data = ps.sub, #PS object
                                fdr_cutoff = 0.05) #pval false discovery value
plot(da_rootchar, total = TRUE, B = 50)

#explore significant taxa
l<-as.data.frame(da_rootchar$significant_taxa)
models<-da_rootchar$significant_models
l

#extract means/SE and create new table
p<-plot(da_rootchar)
d<-as.data.frame(p$data)

da_rootchar$p
da_rootchar$significant_models










#
#
# Fungi 2020
# Each Depth Individually 
#
#
#


##prep phyloseq/subset by depth

ps.sub.top = subset_samples(ps.sub, depth=="0" )
ps.sub.mid = subset_samples(ps.sub, depth=="33" )
ps.sub.bot = subset_samples(ps.sub, depth=="66" )


#create reference table of taxa to ASV
all.tax<-as.data.frame(ps.sub@tax_table@.Data)
phy.gen.all<-cbind(rownames(all.tax), all.tax$sequence, all.tax$Phylum,all.tax$Family,all.tax$Genus)
phy.gen<-as.data.frame(unique(phy.gen.all))
colnames(phy.gen)<-c("asv","seq", "Phylum", "Family","Genus")

##DA Top
###NOTE model overspecified with blk*rs

set.seed(1)
da_top<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.top, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_top$significant_taxa)

#extract means/SE and create new table
p.top<-plot(da_top)
d.top<-as.data.frame(p.top$data)


##DA Mid
###NOTE model overspecified with blk*rs

set.seed(1)
da_mid<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.mid, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_mid$significant_taxa)

#extract means/SE and create new table
p.mid<-plot(da_mid)
d.mid<-as.data.frame(p.mid$data)


##DA bot
###NOTE model overspecified with blk*rs

set.seed(1)
da_bot<- differentialTest(formula = ~ gc+rs+blk, #the DA formula - if you have a controlling factor it must be included. 
                          phi.formula = ~ gc+rs+blk, #the DV formula - if you have a controlling factor it must be included
                          formula_null = ~ blk, # DA controlling factors
                          phi.formula_null = ~blk, #DV controlling factors
                          test = "Wald", boot = FALSE, #wald test is "standard"
                          data = ps.sub.bot, #PS object
                          fdr_cutoff = 0.05) #pval false discovery value
#see significant taxa
as.data.frame(da_bot$significant_taxa)

#extract means/SE and create new table
p.bot<-plot(da_bot)
d.bot<-as.data.frame(p.bot$data)

d.top
d.mid
d.bot

















