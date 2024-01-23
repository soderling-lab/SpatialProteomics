#!/usr/bin/env Rscript 

# title: SwipProteomics
# description: generate protein co-variation (correlation) network
# author: twab

library(stringr)


## ---- Input:
root <- "~/Documents/SoderlingLab/SpatialProteomics"
inputfile = "Vps35_adjm.rda"

input_adjm <- file.path(root,"rdata",inputfile)
gene_name <- strsplit(inputfile, "_")[[1]][1]

print(gene_name)
`%notin%` <- function(x, y) !(x %in% y)


# species from which to collect PPIs for PPI graph
# these are taxonomy ids. human, norway rat,house mouse.. do we need human?
os_keep = c(9606, 10116, 10090)


## ---- Output:
# ppi_adjm.rda

# NOTE: large ouput files saved in root/rdata bc too big to be tracked by git

stopifnot(file.exists(input_adjm))


## ---- prepare the working environment

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load data in root/data
data(Vps35_gene_map) ## this is not mandatory, but calculated in analysis/2_SWIP-TMT or analysis/1_WASH-iBioID

# load data in root/rdata
load(input_adjm)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(getPPIs) # twesleyb/getPPIs
  library(data.table)
  # library(impute) not in twesleyb
})


# load mouse PPIs compiled from HitPredict
data(musInteractome)


## ---- create ppi network

# NOTE: PPIs are NOT used to identify communities

# map uniprot to entrez
uniprot <- colnames(adjm)
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx] # accessing entrez column in gene_map with the necessary ids

# given entrez, collect ppis from musInteractome
ppi_df <- musInteractome %>% 
	 subset(osEntrezA %in% entrez & osEntrezB %in% entrez)  %>% 
	subset(Interactor_A_Taxonomy %in% os_keep) %>%
	subset(Interactor_B_Taxonomy %in% os_keep) %>%
	dplyr::select(osEntrezA, osEntrezB)

# map back to uniprot and cast to matrix for igraph
idx <- match(ppi_df$osEntrezA,gene_map$entrez) # matches/outputs index of osEntreza in genemap entrez
idy <- match(ppi_df$osEntrezB,gene_map$entrez)
ppi_dm <- ppi_df %>% # creates new dataframe with cols ProtA and ProtB and 
                    # contains the uniprot values based on the indexes
	dplyr::mutate(ProtA = gene_map$uniprot[idx], ProtB = gene_map$uniprot[idy]) %>%
	dplyr::select(ProtA,ProtB) %>% 
	as.matrix() # without this it is a dataframe, now a matrix

# create igraph graph, directed as in does not have an arrow
g <- igraph::graph_from_edgelist(ppi_dm, directed=FALSE) # can onlytake matrix

# simplify (weight=0,1) and get the adjacency matrix
ppi_adjm <- as.matrix(igraph::as_adjacency_matrix(igraph::simplify(g)))

# collect proteins that are missing
missing_prots <- colnames(adjm)[colnames(adjm) %notin% colnames(ppi_adjm)]
# these proteins are unconnected^, but we include them in the ppi_adjm so the
# networks have matching vertex sets

# add missing cols
tmp_cols <- matrix(0, nrow=nrow(ppi_adjm),ncol=length(missing_prots))
colnames(tmp_cols) <- missing_prots
rownames(tmp_cols) <- rownames(ppi_adjm)
tmp_dm <- cbind(ppi_adjm,tmp_cols)

# add missing rows
tmp_rows <- matrix(0, nrow=length(missing_prots),ncol=ncol(tmp_dm))
colnames(tmp_rows) <- colnames(tmp_dm)
rownames(tmp_rows) <- missing_prots

# full(ppi)_adjm
full_adjm <- rbind(tmp_dm,tmp_rows)

# sort rows and cols to match adjm
ppi_adjm <- full_adjm[colnames(adjm),rownames(adjm)]


## ---- save networks

# coerce to data.table and write to file
ppi_dt <- as.data.table(ppi_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata", paste0(gene_name, "_ppi_adjm.csv"))
data.table::fwrite(ppi_dt, myfile)
message("saved: ", myfile)


## ---- save as rda

# ppi adjm
myfile <- file.path(root,"rdata", paste0(gene_name, "_ppi_adjm.rda"))
save(ppi_adjm, file=myfile,version=2)
message("saved: ", myfile)
