# schoener-index-calculation
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
````
