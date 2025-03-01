#!/usr/bin/env Rscript

# title: SwipProteomics
# description: module-level analysis with mixed models
# author: Tyler W Bradshaw

### Doing with abundance  termed log2(rel_intensity)
## ---- Inputs

# Input data in root/data/
root <- "~/Documents/SoderlingLab/SpatialProteomics"


## ---- Prepare the R environment

devtools::load_all(root, quiet=TRUE)

gene_name = "Vps35"
# load the data
data(Vps35_tmt) # change to your gene 
data(Vps35_gene_map)# change to your gene 
data(Vps35_ne_adjm_partition) # change it to your gene 

# imports
suppressPackageStartupMessages({
  devtools::install_github("soderling-lab/tidyProt") # NEED THIS FOR getContrast 
  library(lmerTest)
	library(dplyr)
	library(data.table)
	library(doParallel)
  library(tidyProt)
  library(openxlsx)
  library(writexl) ## you can install using install.packages("writexl")
})


## ---- Function

# fit mixed-model to log2 relative (scaled to sum) Intensity
# 1 indicates only random intercept is being evaluated with "Protein"
# Condition is independent variable.
fx <-  log2(Rel_Intensity) ~ 0 + Condition + (1|Protein)

fitModule <- function(prots, tidy_prot, fx) {
  # build list of input args for lmerTest
  lmer_args <- list()
  lmer_args[["formula"]] <- fx
  # only looking at WASH proteins in this data key
  lmer_args[["data"]] <- tidy_prot %>% subset(Protein %in% prots)
  # fit the model with some lmer control
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- do.call(lmerTest::lmer, lmer_args)
  # assess overall contrast and collect results
  LT <- tidyProt::getContrast(fm,"tg","wt") # which column to look at
  result <- lmerTestContrast(fm,LT) %>%
	  mutate(Contrast='tg-wt') %>% unique() %>%
	  mutate('nProts'=length(prots))
  return(result)
} #EOF


## ---- main

modules <- split(names(partition), partition) #  partition at end will also work., if not names, rownames
names(modules) <- paste0("M",names(modules))

message("k Modules: ", length(modules))

## ---- loop to fit module-level models ad assess contrast

# now looking at gene TMT NOT PARTITIONED
# scale Intensity,... looking at abundance.
tidy_prot <- peptides %>%
  group_by(Protein) %>%
  mutate(Rel_Intensity=Intensity/sum(Intensity))

# examine fit of wash complex proteins
# USER: you may do the following, but you must change the analysis to your gene,
# Change wherever you see 'Wash'

# washc_prots <- gene_map$uniprot[grepl("Washc*", gene_map$symbol)] # "Wash*"
# # get only washc prots in tidy_prot
# prots = washc_prots[washc_prots %in% tidy_prot$Protein]
# prots only looks at washc4, washc2, washc5, washc1
# passing in all data into prots(4 wash complexes) and then fitting with lmer

# CHECK <- tidy_prot %>%  subset(Protein %in% prots)

# fm = lmerTest::lmer(fx, tidy_prot)
fm = lmerTest::lmer(fx, tidy_prot)

print(fm)


# why is this being done for just WASH complexes
LT = getContrast(fm,"tg","wt")

lmerTestContrast(fm,LT) %>%
	mutate(Contrast = 'tg-wt') %>%
	mutate(Pvalue = formatC(Pvalue)) %>%
	mutate(nProts = nrow(tidy_prot)) %>%
	unique() %>% knitr::kable()


# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to do module-level analysis
results_list <- foreach(module = names(modules))  %do% { #dopar
  fitModule(modules[[module]], tidy_prot, fx)
} # EOL

names(results_list) <- names(modules)

message("running: ")

## collect results
results_df <- bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>%
	arrange(Pvalue)


## ---- save results

# drop singular col
results_df$isSingular <- NULL

# re-arrange column order
results_df <- results_df %>%
dplyr::select(Module, nProts, Contrast, log2FC,
		percentControl, SE, Tstatistic,
		Pvalue, FDR, Padjust, DF, S2)

# annotate candidate sig modules (Bonferroni padjust < 0.05)
results_df <- results_df %>%
	mutate(candidate = Padjust < 0.05) %>% # significant candidate if p value is 
  # less than 0.05.... contrast is signfiicant and not by chance.
	arrange(desc(candidate))

# summary
message("n Sig modules: ", sum(results_df$candidate))

# list of results

# data.frame describing network partition
idx <- match(names(partition),gene_map$uniprot)
df <-  data.table(UniProt = names(partition),
    Entrez = gene_map$entrez[idx],
    Symbol = gene_map$symbol[idx],
    Membership = partition)

df[, Symbol := as.character(Symbol)]
# results list:
results_list <- list()
results_list[["Partition"]] <- df %>% arrange(Membership)
results_list[["Module Results"]] <- results_df


# save in root/tables
myfile <- file.path(root,"tables", paste0(gene_name, "-TMT-Module-Results.xlsx")) ## This is ranked by log2(Relative_Intensity)
write_xlsx(results_list, myfile)
message("saved :", myfile)

# save results as rda in root/data
# MODULE RESULTS IS CHANGING HERE TO RESULTSDF
module_results <- results_df
myfile <- file.path(root,"data", paste0(gene_name, "module_results.rda"))
save(module_results, file=myfile, version=2)
message("saved :", myfile)

# save module results to csv file
myfile = file.path(root,"data", paste0(gene_name, "module_results.csv"))
module_results_table <- as.data.table(module_results)
data.table::fwrite(module_results_table, myfile)
message("saved :", myfile)


# save sig modules
sig_modules <- module_results$Module[module_results$candidate]
myfile <- file.path(root,"data", paste0(gene_name, "sig_modules.rda"))
save(sig_modules, file=myfile, version=2)
message("saved :", myfile)

