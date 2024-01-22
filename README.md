# SpatialProteomics

Taken from source [SWIP-Proteomics by T. Wesley & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file)

Please refer to the SWIP-proteomics github for the foundation of this program. Below are the specific details on how to get variation analysis for each protein individually, and then protein modules.

Refer to the following paper for methods and further conceptual understanding of the pipeline. [Genetic disruption of WASHC4 drives endo-lysosomal dysfunction and cognitive-movement impairments in mice and humans](https://elifesciences.org/articles/61590) 


**In the event you run into environment/packaging issues:**

<u> R: <u>
If you are missing a package:

`install.packages("rvest")`

`BiocManager::install("AnnotationDbi")`

If you cannot find the package, it is likely a github repository. Please check soderling-lab repositories, the package will likely be in there.

`devtools::install_github("user/{githubrepo}")`


## Necessary Pre-processing
Given the Tandem Mass Taggged (TMT) Mass Spectrometry data from the Proteomics Core, please use the Normalized data sheet moving further.

Please run '0_prepData_transform.py' from 'analysis/2_SWIP-TMT' or 'analysis/3_Clustering' with your target dataset, and prep the data to start at '1_generate-network.R', '1_MSstatsTMT-analysis.R'

This will transform the data to a long format so 

## MSstatsTMT
### Assess differential protein abundance for intrafraction comparisons between WT and MUT

Navigate into and run: 
- analysis/2_SWIP-TMT/1_MSstatsTMT-analysis.R
    - output: adjacency matrix, neten adjacency matrix, TMT protein data
- analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R

### Protein clustering using Leiden Algorithm and Linear Mixed Effect Models
Navigate into analysis/3_Clustering and run:
- analysis/3_Clustering/1_generate-network.R
- analysis/3_Clustering/2_leidenalg-clustering.py
- analysis/3_Clustering/3_post-leidenalg.R
- analysis/3_Clustering/4_generate-PPI-network.R

You can load these packages using `library(PKG)`
