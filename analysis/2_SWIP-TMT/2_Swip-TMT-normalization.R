#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw
renv::init()
# 
# ################## What you need installed ###########################
# install.packages("lmerTest")
# install.packages("dplyr")
# install.packages("data.table")
# install.packages("reshape2")
# install.packages("BiocManager")
# BiocManager::install("AnnotationDbi")
# install.packages("rvest")
# install.packages("httr")
# install.packages("xml2")
# # 
# # # ## Things to load from github
# devtools::install_github("soderling-lab/tidyProt")
# devtools::install_github("soderling-lab/getPPIs")
# devtools::install_github("soderling-lab/geneLists")
# # 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("UniProt.ws")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("UniProt.ws")
# BiocManager::install("org.Mm.eg.db")
# install.packages("GenomeInfoDbData")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")

# library(org.Hs.eg.db)
### ---------------------------------- ###
### 2023
## This file has now been modified to adapt to a new experiment on Vps35,
## LRRK2, SNCA. You may also change the input (genes/protein data) and analyze your data specifically.
## 
## Please head to the original program for just this file as it has
## been heavily modified to deal with 
## updated data from proteomics core. 
## If you have PSMs, this can be a guide but check the original implementation.
##
## ---- input

# project's root directory
root = "~/Documents/SoderlingLab/SpatialProteomics"

# INPUT data is zipped up in root/data/
# zip_file = "TMT.zip"
# input_meta = "TMT-samples.csv"
# input_data = "TMT-raw-peptide.csv"

# input_data is already in format of output 'SWIP_TMT'. Running this for gene_map
input_data = "10415/Transformed_KinSub10415_unbiased_normalized_041124_AnnotatedHumanKinome.csv"
# input_data = "data_JBK10297/10297_SupplementalData_Proteome_122723_Transformed.csv"
# ---- functions

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.


## ---- prepare the R workspace
# prepare the R workspace for the analysis
library(devtools)
library(UniProt.ws)
library(RCurl)

library(UniProt.ws)
library(AnnotationDbi)
########-----------------##############
# If you cannot load any of the below packages please install them,
# if they do not exist. Head to the soderling-lab github repository and install the github repo to your R console.
# You must download :

#.    among the remaining indicated by the git repo.
########-------------------#######

# load required packages and functions
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(getPPIs) # soderling-lab/getPPIs for mouse PPIs
	library(geneLists) # soderling-lab/geneLists for gene mapping
	library(data.table) # for working with tables
	# library(doParallel) # for parallel processing
#   library(biomaRt)
})

# load project specific functions and data
devtools::load_all(root, quiet=TRUE)

# project directories:
datadir <- file.path(root, "data") # Key pieces of data saved as rda
rdatdir <- file.path(root, "rdata") # Temporary data files
tabsdir <- file.path(root, "tables") # Output tables saved as excel files
downdir <- file.path(root, "transformeddata") # Misc downloads/temporary files

# create project output directories if necessary
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

# extract the raw TMT data from zipped file
myfile <- file.path(downdir, input_data) # all peptide information
peptides <- data.table::fread(myfile)

## ---- map all Uniprot Protein numbers to stable entrez IDs

message("\nCreating gene identifier map.")

# first, remove any non-mouse proteins from the data
# peptides <- peptides %>% filter(grepl("Mus musculus",Origin))

# remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(!Protein %in% ig_prots)

# also remove Keratins
peptides <- peptides[!grepl("Keratin", peptides$Function), ]

# collect all Uniprot IDs
uniprot <- sub(";.*$|__.*$", "", peptides$Protein)
uniprot <- unique(uniprot)
# map Uniprot IDs to Entrez using online MGI batch query function
## Mouse
# entrez <- geneLists::queryMGI(uniprot) ## MOUSE
# names(entrez) <- uniprot

## Human
up <- UniProt.ws(taxId=9606) #mouse 10009

entrezd <- AnnotationDbi::select(up, keys=uniprot, columns=c("xref_geneid"), # xref_geneid
                  keytype="UniProtKB")
entrezd$GeneID <- gsub("[;\"]", "", entrezd$GeneID)
entrezd$GeneID  <- as.integer(entrezd$GeneID)
entrez <- setNames(rep(NA_integer_, length(uniprot)), uniprot) # initialize vector w all uniprots
# Using vectorized indexing to update 'entrez'
entrez[entrezd$From] <- entrezd$GeneID

## check all Protein ids in peptides have been mapped appropriately to entrez.
# peptides_in_entrez <- all(peptides$Protein %in% mgi_map$UniProt)
# missing_ids <- c(setdiff(peptides$Protein, mgi_map$UniProt))
prots_naEntrez <- entrez[is.na(entrez)]
num_na_entrez <- sum(is.na(entrez))

message("The following length of ids are in peptides (raw) but Entrez ID was not found: ",paste(num_na_entrez, collapse = ", "))

missEntrezd <- select(up, keys=names(prots_naEntrez), columns=c("xref_geneid"), # xref_geneid
                  keytype="UniProtKB")
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

na_entrezX2sum <- sum(is.na(entrez))
na_entrezX2 <- entrez[is.na(entrez)]
message("Fixed some! Now this # are still missing",paste(na_entrezX2sum, collapse = ", "))
# map any remaining missing IDs by hand
message("Mapping missing IDs by hand.\n")
mapped_by_hand <- c(Q6A1A2 = 5170,
A0A0J9YX94 = 105373377,
A0A3B3IU46 = 353267,
A2A3N6 = 266971,
A6NKF1 = 29901,
A8CG34 = 9883,
A8MPP1 = 100302090,
C4AMC7 = 100287171,
E9PAV3 = 4666,
L0R819 = 110599588,
O75420 = 2887,
P01859 = 3501,
P0CG12 = 113455421,
P20742 = 7040,
P49674 = 324,
P51784 = 5971,
Q14CW9 = 56970,
Q5T1J5 = 51142,
Q5T3I0 = 54865,
Q63ZY6 = 260294,
Q68D20 = 441194,
Q6ZN08 = 7617,
Q8NI35 = 6331,
Q96PV7 = 54540,
Q9NQA3 = 100287171,
Q9UF83 = 101060157,
Q9UJX3 = 999,
Q9UPS6 = 23067,
Q9Y4D8 = 283450,
P00761 = 100302368)
entrez[names(mapped_by_hand)] <- mapped_by_hand


# entrez[names(mapped_by_hand)] <- mapped_by_hand
print(paste0("how many uniprots werent mapped? ",sum(is.na(entrez))))
## remove blanks (what else to do?)
# entrez<- na.omit(entrez)

# check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
# if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }

# map entrez ids to gene symbols using twesleyb/geneLists
# NOTE: getIDs is just an easier-to-use wrapper around AnnotationDbi::mapIDs
# NOTE: you need the org.Mm.eg.db package from Bioconductor for mapping mouse genes
# can only do getIDs for mouse. 
gene_symbols <- mapIds(org.Hs.eg.db, keys=as.character(entrez), 
			column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# gene_symbols <- select(up, keys=names(entrez), columns=c(""), # xref_geneid
#                   keytype="GeneID")
# create gene identifier mapping data.table
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = gene_symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")


####################################################################
## 2023. PP--- WE DO NOT DO. DATA IS ALREADY TRANSFORMED IN 0_prepData_transform.py 
## Deleted all unnecessary calls, please refer to Tyler's original implementation for fruther details.
####################################################################

peptides <- peptides %>%
	# other munge
	mutate(Genotype = Genotype) %>%
  mutate(Gene = Gene) %>%
  mutate(Function = Function) %>%
	mutate(Protein = Protein) %>%
	mutate(Mixture = Mixture)%>%
	mutate(BioFraction = BioFraction) %>%
  mutate(Intensity = Intensity)%>%
	mutate(Abundance = Abundance) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
  mutate(Relative_Intensity = Relative_Intensity) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,
	              Intensity, Abundance, Relative_Intensity) %>%
	# calculate relative_Intensity (sum normalization)str
  group_by(Protein)

## ----  save key results

# # final normalized protein in tidy format as rda object
#
# names(peptides)[names(peptides) == "Protein"] <- "Protein"

# peptides <- peptides %>% 
#   dplyr::select(Protein,Gene, Origin, Function, Mixture,Genotype,BioFraction,Abundance)

# peptides <- data.frame(peptides)
# peptides <- peptides %>%
# group_by(Protein)


gene_name <- sub("^[^_]*_([^_]+)_.*$", "\\1", input_data)

# peptides$Protein <- peptides$Protein

myfile <- file.path(datadir, paste0(gene_name, "_tmt.rda"))
save(peptides,file=myfile,version=2)
message("saved: ", myfile)

# save gene_map
myfile <- file.path(datadir,paste0(gene_name, "_gene_map.rda"))
save(gene_map,file=myfile,version=2)
message("saved: ", myfile)
