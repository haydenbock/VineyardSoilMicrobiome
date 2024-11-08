#Package Management
BiocManager::install("dada2", version = "3.19")
library(dada2); packageVersion("dada2")
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(readr)

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
packages <- c("nlme","lme4", "effects","corncob","magrittr","car", "interactions","RVAideMemoire","tibble","doParallel","phyloseq", "ggplot2","vegan","microbiome","microbiomeSeq",
              "dplyr","Biostrings","ShortRead","tmap","dada2","jtools","ggnewscale","RColorBrewer","limma")
ipak(packages)


#Import Phyloseq object(s)


#calculations 

counts<-as.data.frame(ps_obj@otu_table)
meta<-as.data.frame(ps_obj@sam_data)

#### Calculate alphadiversity ####
rich <- estimate_richness(ps_obj,measures="Shannon")
richtable<-cbind(meta,rich)


