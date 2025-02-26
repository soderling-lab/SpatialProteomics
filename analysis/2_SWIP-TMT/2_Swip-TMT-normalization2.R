#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw

# Initialize renv if necessary
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::init()

# Necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}
if (!requireNamespace("UniProt.ws", quietly = TRUE)) {
  BiocManager::install("UniProt.ws")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load required packages
library(UniProt.ws)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(dplyr)

# Project's root directory
root <- "~/Documents/SoderlingLab/SpatialProteomics"

# Input data file
input_data <- "Transformed_LOPIT_LRRK2_young.csv"

# Project directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "transformeddata")

# Create project output directories if necessary
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

# Load peptide data
myfile <- file.path(downdir, input_data)
peptides <- data.table::fread(myfile)

# Map all Uniprot Protein numbers to stable entrez IDs
message("\nCreating gene identifier map.")

# Remove unwanted proteins
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750","P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(!Protein %in% ig_prots)
peptides <- peptides[!grepl("Keratin", peptides$Function), ]

# Collect all Uniprot IDs
uniprot <- sub(";.*$|__.*$", "", peptides$Protein)
uniprot <- unique(uniprot)

# Map Uniprot IDs to Entrez using UniProt.ws
up <- UniProt.ws(taxId=9606)
entrezd <- AnnotationDbi::select(up, keys=uniprot, columns=c("xref_geneid"), keytype="UniProtKB")
entrezd$GeneID <- gsub("[;\"]", "", entrezd$GeneID)
entrezd$GeneID  <- as.integer(entrezd$GeneID)
entrez <- setNames(rep(NA_integer_, length(uniprot)), uniprot)
entrez[entrezd$From] <- entrezd$GeneID

# Handle missing Entrez IDs
prots_naEntrez <- entrez[is.na(entrez)]
num_na_entrez <- sum(is.na(entrez))
message("The following length of ids are in peptides (raw) but Entrez ID was not found: ", paste(num_na_entrez, collapse = ", "))

missEntrezd <- AnnotationDbi::select(up, keys=names(prots_naEntrez), columns=c("xref_geneid"), keytype="UniProtKB")
missEntrezd$GeneID <- sapply(missEntrezd$GeneID, function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    first_id <- strsplit(as.character(x), ";")[[1]][1]
    return(as.integer(first_id))
  }
})
missEntrezd$GeneID  <- as.integer(missEntrezd$GeneID)
entrez[missEntrezd$From] <- missEntrezd$GeneID

# Manually map remaining missing IDs
mapped_by_hand <- c(Q6A1A2 = 5170, A0A0J9YX94 = 105373377, A0A3B3IU46 = 353267, A2A3N6 = 266971, A6NKF1 = 29901, A8CG34 = 9883, A8MPP1 = 100302090, C4AMC7 = 100287171, E9PAV3 = 4666, L0R819 = 110599588, O75420 = 2887, P01859 = 3501, P0CG12 = 113455421, P20742 = 7040, P49674 = 324, P51784 = 5971, Q14CW9 = 56970, Q5T1J5 = 51142, Q5T3I0 = 54865, Q63ZY6 = 260294, Q68D20 = 441194, Q6ZN08 = 7617, Q8NI35 = 6331, Q96PV7 = 54540, Q9NQA3 = 100287171, Q9UF83 = 101060157, Q9UJX3 = 999, Q9UPS6 = 23067, Q9Y4D8 = 283450, P00761 = 100302368)
entrez[names(mapped_by_hand)] <- mapped_by_hand

message("Fixed some! Now this # are still missing: ", paste(sum(is.na(entrez)), collapse = ", "))

# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys=as.character(entrez), column="SYMBOL", keytype="ENTREZID", multiVals="first")

# Create gene identifier mapping data.table
gene_map <- data.table(uniprot = names(entrez), entrez = entrez, symbol = gene_symbols)
gene_map$id <- paste(gene_map$symbol, gene_map$uniprot, sep="|")

# Preprocess peptide data
peptides <- peptides %>%
  mutate(Genotype = Genotype) %>%
#   mutate(Gene = Gene) %>%
  mutate(Function = Function) %>%
  mutate(Protein = Protein) %>%
  mutate(Mixture = Mixture) %>%
  mutate(BioFraction = BioFraction) %>%
  mutate(Intensity = Intensity) %>%
  mutate(Abundance = Abundance) %>%
  mutate(Condition = interaction(Genotype, BioFraction)) %>%
  mutate(Relative_Intensity = Relative_Intensity) %>%
  dplyr::select(Protein, Mixture, Genotype, BioFraction, Condition, Intensity, Abundance, Relative_Intensity) %>%
  group_by(Protein)

# Save key results
gene_name <- sub("^[^_]*_[^_]+_([^_]+).*$", "\\1", input_data)

myfile <- file.path(datadir, paste0(gene_name, "_tmt.rda"))
save(peptides, file=myfile, version=2)
message("saved: ", myfile)

myfile <- file.path(datadir, paste0(gene_name, "_gene_map.rda"))
save(gene_map, file=myfile, version=2)
message("saved: ", myfile)