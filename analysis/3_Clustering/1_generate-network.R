#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: generate protein co-variation (correlation) network and perform
#   network enhancement
## Run 2_SWIP/2_SWIP-TMT-normalization.R first with transformed data to get needed data tables!!
## downloads
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
library(devtools)
if (!requireNamespace("igraph", quietly = TRUE)) {
    install.packages("igraph" )
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2" )
}
if (!requireNamespace("neten", quietly = TRUE)) {
    devtools::install_github("soderling-lab/neten")
}
## ---- Input:
root <- "~/Documents/SoderlingLab/SpatialProteomics"
data_dir <- file.path(root, "data")
# setting working directory manually in the event session restarts
setwd("~/Documents/SoderlingLab/SpatialProteomics/analysis/3_Clustering")

## Change to your gene name
gene_name <- "" # change
load(file.path(data_dir, "KinSub10415_gene_map.rda"))
load(file.path(data_dir, "KinSub10415_tmt.rda"))                                                                                                                                             
## ---- Output:

# * adjm.rda
# * ne_adjm.rda

# * ne_adjm.csv --> for leidenalg clustering!

# NOTE: large ouput files (>100mb) saved in root/rdata bc too big to be tracked by git


## ---- prepare the working environment

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

#################################################
####### You must change all below to your specific gene needs ##############

# load data in root/data


################################################################
### If

input_data = "10415/Transformed_KinSub10415_unbiased_normalized_041124_AnnotatedHumanKinome.csv"
downdir = "../../transformeddata"  # Step back twice and then go to transformeddata
myfile <- file.path(downdir, input_data) # all peptide information
novel <- data.table::fread(myfile)

gene_name <- sub("^[^_]*_([^_]+)_.*$", "\\1", input_data)

### Where 
# myfile <- file.path("../../data", paste0(gene_name, "_tmt.rda"))
# save(peptides,file=myfile,version=2)
# message("saved: ", myfile)
# ########################

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(neten) # twesleyb/neten
  library(igraph)
  library(data.table)
})


## ---- create covariation network

message("Generating covariation network...")

# peptides <- peptides %>%
#   subset(select = -c(V1, Gene, Function, Origin)) %>%3
#   ungroup()

# no median summarization of bioreplicates
# network is constructed from log2(Intensity) ~ Abundance
## adjust for duplicate accessions&fractions
peptides <- peptides %>%
  group_by(Protein, Mixture, Genotype, BioFraction, Condition) %>%
  summarize(
    Intensity = sum(Intensity, na.rm = TRUE), ## MEAN?
    Abundance = sum(Abundance, na.rm = TRUE),
    Relative_Intensity = sum(Relative_Intensity, na.rm = TRUE), # sum checked, and correct.
    .groups = 'drop'
  )
  
dm <- peptides %>% ## need to adjust for sum
  reshape2::dcast(Protein ~ Mixture + Genotype + BioFraction, value.var = "Abundance") %>%
  as.data.table() %>%
  as.matrix(rownames="Protein")

# there are a small number proteins with some missing vals
# e.g. Q9QUN7 = low abundance only quantified in 4/7 fractions
# remove these proteins
idx <- apply(dm,1,function(x) any(is.na(x)) | any(x <0) ) ## also removing negative abundances ....
warning(sum(idx)," proteins with any missing or negative values are removed.")
filt_dm <- dm[!idx,] # picks idx that are Not True! not true= not missing
print(!any(filt_dm<0))
stopifnot(!any(filt_dm<0)) # checks if any values are less than 0


## calculate correlation matrix
adjm <- cor(t(filt_dm), method="pearson",use="complete.obs")

# save to file
currentdir = getwd()
print(currentdir)
enclosingfolder = dirname(currentdir)
root = dirname(enclosingfolder)

# stop()
## ---- network enhancement

# Wang et al., 2018 (Nature Communications; PMID:30082777)

message("Performing network enhancement...")
folder = "rdata"
if (!dir.exists(file.path(root, folder))) {
  dir.create(file.path(root, folder))
}
print(paste("created directory", file.path(root, folder)))
# FIXME: is there any room for speed improvements here? Probably...
ne_adjm <- neten::neten(adjm) # result is robust to neten parameters
Sys.time()

## ---- save networks as csv in rdata
# coerce to data.table and save adjm.csv
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root, folder, paste0(gene_name, "_adjm.csv"))
print(paste("Writing adjm_dt to", myfile))
data.table::fwrite(adjm_dt, myfile)
message("saved: ", myfile)

# coerce to data.table and save ne_adjm.csv
ne_adjm_dt <- as.data.table(ne_adjm,keep.rownames="Protein")
myfile <- file.path(root, folder,paste0(gene_name, "_ne_adjm.csv"))
print(paste("Writing ne_adjm_dt to", myfile))
data.table::fwrite(ne_adjm_dt, myfile)
message("saved: ", myfile)


## ---- save data as rda

# NOTE: data are saved in root/rdata
# adjm could be make smaller by melting to edge list
# but it is still too big at ~ 150 mb
# ne_adjm is smaller bc it is sparse, but still to large to be easily tracked by
# git

# adjm

myfile <- file.path(root, folder, paste0(gene_name, "_adjm.rda"))
print(paste("Writing adjm vs2 to", myfile))
save(adjm, file=myfile,version=2)
message("saved: ", myfile)

# ne adjm
myfile <- file.path(root, folder, paste0(gene_name,"_ne_adjm.rda"))
save(ne_adjm, file=myfile,version=2)
message("saved: ", myfile)
print(paste("Writing ne_adjm vs2 to", myfile))
