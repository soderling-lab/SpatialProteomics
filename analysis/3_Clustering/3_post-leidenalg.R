#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## project root dir
root <- "~/Documents/SoderlingLab/SpatialProteomics"

## input in root/rdata
input_part <- "KinSub10415_partition.csv" # change

## output is [input_part].rda
output_part <- paste0(tools::file_path_sans_ext(input_part),".rda")
output_csv <-  paste0(tools::file_path_sans_ext(input_part),".csv")

## ---- prepare the R env

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

## ---- load partition


myfile <- file.path(root,"rdata",input_part)
df <- data.table::fread(myfile,drop=1)

partition_list <- unlist(apply(df,1,function(x) list(x)),
			 recursive=FALSE, use.names=FALSE)
print(length(partition_list))
# add 1 bc python is 0-based
part <- partition_list[[1]] + 1

# set small modules to 0
modules <- split(names(part),part)
too_small <- as.numeric(names(which(sapply(modules,length) < 5)))
part[part %in% too_small] <- 0


## status

message("Total number of modules: ", length(unique(part))-1)


## ---- save as rda

# save partition
partition <- data.frame(Module = part)
colnames(partition)[colnames(partition) == "x"] <- "Module"

partition <- arrange(partition, Module)


myfile <- file.path(root,"data", output_part)
print(output_part)
save(partition, file = myfile, version = 2)
message("saved: ", myfile)
# save as csv
myfile <- file.path(root, "data", output_csv)
write.csv(partition, file = myfile, row.names = TRUE) # row.names = FALSE is often desirable to prevent row names being saved as an extra column
message("saved: ", myfile)
