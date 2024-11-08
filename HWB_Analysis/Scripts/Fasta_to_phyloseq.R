#09/27/2024
#Basic Overview (no coding) - https://www.youtube.com/watch?v=TnWeVUx5ERU
#Phyloseq- https://www.youtube.com/watch?v=YRZC3henVPE&t=40s
#https://microbiome.netlify.app
#https://bioconductor.org/install/


library(phyloseq)

#Step 1 ----
#https://www.youtube.com/watch?v=YRZC3henVPE&t=40s


### Data structure of phyloseq object
data(GlobalPatterns) # Loading included example data
# The DADA2 sequence inference method is reference-free

###### Part I: Creating a phyloseq object
### 1. Prepare otu_table from raw reading fastq files
BiocManager::install("dada2", version = "3.19")
library(dada2); packageVersion("dada2")


#2x250 Illumina Miseq amplicon sequencing of the V4 region of the 16S rRNA from mouse gut
# http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar, unzip the data
path <- "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/Hayden_Learning/StabilityNoMetaG/"
list.files(path)

miseq_path <- file.path("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/Hayden_Learning/StabilityNoMetaG/")
list.files(miseq_path)
filt_path <- file.path("~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/Hayden_Learning/StabilityNoMetaG/filtered/")
list.files(filt_path) #in your miseq folder, you need to make a blank folder called "filtered" for this to work

# Forward and reverse fastq file format: SAMPLENAME_R1_001.fastq & SAMPLENAME_R2_001.fastq
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# Inspect read quality profiles
plotQualityProfile(fnFs[1:2]) # visualize the quality profiles of the forward reads
plotQualityProfile(fnRs[1:2]) # visualize the quality profiles of the reverse reads

### Trim and Filter
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(245, 160), 
                    maxN=0, maxEE=2, truncQ=2, compress=TRUE)
}


#NOTES: On line 42 in the trunc() function, the 245 comes from the forward reads where the quality drops (visualized throught the figure), and the 
#       160 comes from where the quality drops on the reverse reads. 

#Step 2: ----
#https://www.youtube.com/watch?v=FMXa_oYLAXM
#stopped at about 4 minutes into the video, about line 70

## Infer sequence variants
# Dereplicate the sequences to remove redundancy
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Name the dereplicated objects by their sample names
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Estimate the error rates from a (sufficiently large) subset of the data 
# to distinguish sequencing errors from real biological variation
ddF <- dada(derepFs[1:40], err=NULL, selfConsist=TRUE)
plotErrors(ddF, nominalQ=TRUE)

ddR <- dada(derepRs[1:40], err=NULL, selfConsist=TRUE)
plotErrors(ddR, nominalQ=TRUE)

# Perform inference from the pooled sequencing reads from all samples (pool=TRUE)
dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE) # take a while to run
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE) # take a while to run

# merge the inferred forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

## Construct sequence table 
otu_table.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

## remove chimeric sequences
otu_table <- removeBimeraDenovo(otu_table.all)



#Creating a phyloseq object-----
#Video: https://www.youtube.com/watch?v=WQvkG_s3G0c

### 2. Prepare tax_table
# A training set of classified sequences
# https://www.dropbox.com/sh/mfcivbudmc21cqt/AAB1l-AUM5uKvjrR33ct-cTXa?dl=0
ref_fasta <- "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/Hayden_Learning/rdp_train_set_14.fa.gz"

# Assign taxonomy: comparing sequence variants to rdp_train_set_14 classified sequences, 
taxtab <- assignTaxonomy(otu_table, refFasta = ref_fasta)

### 3. Prepare sample_data
# https://web.stanford.edu/class/bios221/data/MIMARKS_Data_combined.csv
sampledf <- read.csv("https://web.stanford.edu/class/bios221/data/MIMARKS_Data_combined.csv", header=TRUE)

#correcting the sample_id column
sampledf$SampleID <- paste0(gsub("00", "", sampledf$host_subject_id), "D", sampledf$age-21)
sampledf <- sampledf[!duplicated(sampledf$SampleID),] # Remove duplicate entries for reverse reads
rownames(otu_table) <- gsub("124", "125", rownames(otu_table)) # Fixing an odd discrepancy
all(rownames(otu_table) %in% sampledf$SampleID) # TRUE

rownames(sampledf) <- sampledf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
               "host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
               "diet", "family_relationship", "genotype", "SampleID") # 14 columns
sampledf <- sampledf[rownames(otu_table), keep.cols]
   
### 4. Construct phylogenetic tree
# The DADA2 sequence inference method is reference-free, so we need to construct 
# the phylogenetic tree in order to relate the inferred sequence variants

## align the sequence using DECIPHER R package
#BiocManager::install("DECIPHER")
library(DECIPHER) 
seqs <- getSequences(otu_table)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

## construct a phylogenetic tree using phangorn R package
#BiocManager::install("phanghorn")
library(phangorn)

# first construct a neighbor-joining tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# next fit a GTR maximum likelihood tree
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0))

### 5. Combine data into a phyloseq object
ps_obj <- phyloseq(otu_table(otu_table, taxa_are_rows = FALSE),
                   tax_table(taxtab), 
                   sample_data(sampledf),
                   phy_tree(fitGTR$tree))

saveRDS(ps_obj, "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/shared_RootAgroEco/research_projects/vineyard_soil_microbiome/Hayden_Learning/ps_obj.rds")



#Analyzing Phyloseq Object-----
