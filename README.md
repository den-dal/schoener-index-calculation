# FILTER data by YEAR
````
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
#library(ape)
library(Biostrings)

setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/rbcL")
list.files()

### convert csv with sequences into fasta if not done already
library(tidyverse)
library(readr)
list.files()
csv = read_csv("rbcL_296zOTUs_sequences&taxonomy.csv")
writeLines(paste0(">", csv$zOTU, "\n", csv$sequence), "D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/rbcL/rbcL_zOTU_sequences&taxonomy.fasta")

### 1. Import data ####
otu_mat    <- read.table("OTUtable_rbcL_min100readspersample.txt", 
                         header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(otu_mat)
meta_df   <- read.table("metadata_rbcL_min100readspersample.txt", 
                        header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(meta_df)
head(meta_df)
tax_df      <- read.table("taxonomy_rbcL_296zOTUs_15082025.txt",   # your taxonomy file
                          header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(tax_df)
head(tax_df)
seqs <- readDNAStringSet("rbcL_zOTU_sequences.fasta", format="fasta")

#to analyze data only from year
meta_df <- subset(meta_df, Year == 2022)
dim(meta_df)
# Filter the zOTUs table accordingly
otu_mat <- otu_mat[, rownames(meta_df), drop = FALSE]
dim(otu_mat)
otu_mat <- otu_mat[rowSums(otu_mat) >0, , drop = FALSE]

### Collapse OTUs to Genus level ----
# Match taxonomy to OTU table by OTU IDs
tax_df_sub <- tax_df[rownames(otu_mat), ]
dim(tax_df_sub)

# Use the Genus column (change name if needed)
genus_vec <- tax_df_sub$Genus

# Make sure matrix is numeric
otu_mat_num <- as.matrix(otu_mat)
mode(otu_mat_num) <- "numeric"

# Sum OTUs by Genus (rows = genera, columns = samples)
genus_mat <- rowsum(otu_mat_num, group = genus_vec)

### 3. Collapse Genus-level reads to Population level ----
# Match metadata to samples (columns in genus_mat)
meta_sub <- meta_df[colnames(genus_mat), ]
pop_vec  <- meta_sub$Population

# Sum samples within each Population (rows = genera, cols = Populations)
genus_pop_mat <- t(rowsum(t(genus_mat), group = pop_vec))

### Save table reads per genus per population ----
write.table(genus_pop_mat,
             file = "2022rbcL_genus_reads_by_Population_min100readspersample.txt",
             sep = "\t", quote = FALSE, col.names = NA)
````
### Schoener index between Populations using spaa ----
````
library(spaa)

### Remove populations with zero total reads ----
genus_pop_mat <- genus_pop_mat[, colSums(genus_pop_mat) > 0, drop = FALSE]

### Convert to proportional read abundances per Population ----
# columns = Populations, each column sums to 1
prop_mat <- sweep(genus_pop_mat, 2, colSums(genus_pop_mat), "/")

### Schoener niche overlap using spaa ----
# spaa expects: rows = resources, columns = species
# here: rows = Genera (resources), columns = Populations (species) -> OK
schoener_dist <- niche.overlap(prop_mat, method = "schoener")

### 4. Convert 'dist' object to full Population × Population matrix ----
schoener_mat <- as.matrix(schoener_dist)

### 5. (optional) Save to file ----
write.table(schoener_mat,
           file = "Schoener_genus_by_Population2022_26112025.txt",
           sep = "\t", quote = FALSE, col.names = NA)
````
