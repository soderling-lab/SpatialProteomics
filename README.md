# SpatialProteomics

Taken from source [SWIP-Proteomics by T. WA Bradshaw & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file)

Please refer to the SWIP-proteomics github for the foundation of this program. Below are the specific details on how to get variation analysis for each protein individually, and then protein modules.

All programming in analysis/ sourced from Tyler WA Bradshaw. Adapted in certain areas to fit current needs/updates.

Refer to the following paper for methods and further conceptuals of the pipeline. [Genetic disruption of WASHC4 drives endo-lysosomal dysfunction and cognitive-movement impairments in mice and humans](https://elifesciences.org/articles/61590)


## Pipeline Overview
Tandem Mass Tag (TMT) is a chemical label that facilitates sample multiplexing in mass spectrometry for the quantification and identification of proteins, often used in large-scale proteomic studies. Each sample mixture contains several fractions (X in number), each labeled with a unique chemical barcode. Each fraction corresponds to different cellular components with different densities. Typically, $\geq 2$ mixtures are analyzed, providing various replicates of WT (Wild Type) & MUT (Mutant) fractions for detailed analysis. In an unperturbed/unmutated environment, every fraction is analyzed, and then a specific gene is knocked out, creating a corresponding MUT fraction. Thus, for X WT fractions in a mixture, there are also X MUT fractions, totaling 2X fractions per mixture.

In this study, 7 WT and 7 MUT fractions were analyzed in each mixture. The TMT proteomic data allows for examining how proteins cluster in terms of abundance or regulation relative to others. Any correlation in abundance change, whether negative or positive, is tracked. An adjacency matrix is created to show how each protein's abundance changes across fractions (WT & MUT) and how this change differs from every other protein's behavior. This matrix is signed positive (+).

The Leiden algorithm is used for community detection, analyzing inter-module and intra-module connections. The 'Surprise' optimizer is chosen for its ability to consider both node weightage and edge count in a module.

With the modules identified, we use a linear mixed model to determine which modules have proteins that differ significantly in WT fractions compared to MUT fractions. A module is deemed significant if the p-value is $\leq 0.05$ between WT and MUT fractions. 

Given these modules, we find which pathway (gene ontology) the proteins in each specific module are enriched in. 

We also plot the averaged mixtures for each module, analyzing how the proteins vary across fractions demonstrating changes within WT, within MUT and WT and MUT compared. 

**Example data files are provided under example_data/ for all steps. Please refer to this to get a sense of the program**

## Environments
### R:
This project uses `renv` for dependency management. To set up the R environment for this project:
1. Open R with the current repo (SpatialProteomics)
    - if you are in bash, type `R` to enter R console
2. Set the working directory to where SpatialProteomics/ was downloaded. `setwd('path/to/SpatialProteomics')`
2. Run `install.packages("renv")`
3. Run `renv::restore(lockfile = "analysis/renv.lock")`. This restores the R environment from renv.lock in `analysis/`.

If you are missing a package:
`install.packages("PKG")`

It may be a package from BiocManger specifically. In that event:
`BiocManager::install("AnnotationDbi")`


If you cannot find the package, it is likely a github repository. Please check soderling-lab repositories, the package will likely be in there.

`devtools::install_github("{github-user}/{githubrepo}")`

You can load these packages using `library(PKG)`

### Python

For running `analysis/3_Clustering/2_leidenalg-clustering.py` you will need a specific python environment. 

In the `~/analysis/` of this repository `SpatialProteomics/`, please run:

`conda env create -f environment.yml` 

and activate:
`conda activate spatial_env`

In the event you run into missing packages, please install with `conda install {pkg}` or `pip install {pkg}`


**NOTE: For the singular python file in the module analysis pipeline. You must activate the environment, and then run it from the command line.**

`python analysis/3_Clustering/2_leidenalg-clustering.py`

#### Gene Ontology env
When you start with the gene ontology pipeline (after finishing gathering all the module data) you must deactivate the `spatial_env` environment `conda deactivate`. 

In the `~/geneontologies/` dir please run:
`conda env create -f go_environment.yml`

and activate:
`conda activate GOenv`

#### ShinyGO environment
If you decide to webscrape with ShinyGo to get pathways for each module. Please deactivate whatever env is active `conda deactivate` and run `conda env create shinygo_environment.yml`

#### plot_MUTvsWT environment
If you would like to create plots comparing mutant vs wildtype fractions, then deactivate any current conda environment `conda deactivate` then run `conda env create -f plot_env.yml`

## Program Instructions

### Necessary Pre-processing
Given the Tandem Mass Taggged (TMT) Mass Spectrometry data from the Proteomics Core, please use the Normalized data sheet (labeled in excel file) moving further.

Please run your target dataset through `0_prepData_transform.py` from `analysis/2_SWIP-TMT/` or `analysis/3_Clustering/` to prep the data to start at `1_generate-network.R` or `1_MSstatsTMT-analysis.R`

In  `0_prepData_transform.py`, navigate into the `main()` function and change the filename and worksheet accordingly. This will transform the data to a long format so the fractions and mixtures are structured vertically rather than horizontally.

**In the event you are starting from the PSMs, please run 0_PD-data-preprocess.R** 
Please refer to the source code [SWIP-Proteomics by T. Wesley & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file) if you need to start here. 

### MSstatsTMT
**Assess differential protein abundance between fractions between WT and MUT models**

Navigate into and run: 


- analysis/2_SWIP-TMT/1_MSstatsTMT-analysis.R
    - output: adjacency matrix, neten adjacency matrix, TMT protein data
     - You can only run this if you have the PSMs. Reach out the proteomics core in the event you need this. Some examples are under `~/PSM`

NOTE: notice that 2_SWIP-TMT-normalization is before 1_MSstatsTMT-analysis. You must get all the output data tables from 2_ before running 1_. Additionally, you will be able to get individual protein data (mut vs wt differences in abundance and significance of change) given you have the Peptide-Spectrum Match (PSMs) from the proteomics core. In the event you do, start with 'analysis/2_SWIP-TMT/0_PD-data-preprocess.R' to analyze the PSMs and get the necessary datasets (pd_psm, pd_annotation, mut_vs_control) etc. 

### Leiden Algorithm
**Spatial Protein clustering in terms of abundance regulation**

First to normalize the data, run `analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R` with the long, transformed data structure.
    - Start with this to get the entire pipeline below for protein clustering
    - saved files will be printed in console

Navigate into analysis/3_Clustering to change variables and run:
- analysis/3_Clustering/1_generate-network.R
- analysis/3_Clustering/2_leidenalg-clustering.py
    - please activate the python environment using environment.yml file and run this from the command line after making your edits. `python 2_leidenalg-clustering.py`
- analysis/3_Clustering/3_post-leidenalg.R
- analysis/3_Clustering/4_generate-PPI-network.R

Given the output from the four above programs, you may now find out which modules are significant in that WT and MUT protein abundances vary within the module.

### Identify significant clusters (modules) with Linear Mixed Models
Navigate into analysis/4_Module-Analysis and run:
- analysis/4_Module-Analysis/1_module-lmerTest-analysis.R

You have found the clustered modules in your dataset! Please navigate into ~/tables, the protein modules are saved under `{YOUR_GENE}_TMT-Module-Results.xlsx`
 - Each protein is assigned to a module. You can also analyze how each protein does individually using MSStatsTMT procedure detailed above. MSStatsTMT will yield how much the abundance of **each** protein changed relative to WT and MUT fraction variances.

### Gene Ontologies
If you want to remain on the pipeline and look for gene ontologies,
you may run the following file the original author used for WASH proteins (changes needed for your GOI):
- analysis/4_Module-Analysis/2_module-GSEA.R

**OR**

You can scrape [ShinyGo](http://bioinformatics.sdstate.edu/go/) using `~/ShinyGo/ShinyGO_analysis/1_analyze_modules.ipynb`. This works in an automated loop but sometimes the website is spontaneous so in intervals it will er.. 
- ie. if you have 50 modules, it may err on the 25th module, but all until 25 will save.

You can create a heatmap with the identified pathways modules using `ShinyGo/ShinyGO_analysis/2_cluster_analysis.py`

**OR** if you have a specific pathway in mind (Parkinson's, Autism..)

Navigate into `~/geneontologies` and run:
 - `geneontologies/get_KEGGnums.ipynb `
    - input: `~/tables/{KOGENE}-TMT-Module-Results.xlsx`
    - gets KEGG ortholog nums to get pathways in next file

 - `geneontologies/get_pathways.ipynb`
    - input: `~/tables/{KOGENE}_moduleResults_KEGGortholog.csv`
    - Change the Pathway of Interest (POI) from Parkinson's to what you are looking for!
    - This will calculate gene_ontologies with your POI using a hypergeometric test. Each module is analyzed at a time, and the overlap between genes in that module vs genes in POI will indicate p-value. (p-val < 0.05 signifies the module is enriched in POI) 
        - P-value is adjusted with benjamin-hochberg for multiple testing (FDR) since multiple memberships are analyzed at once.
        - output: `~/enrichments_DIY/{YOURGENE_KO}_{YOUR_POI}_hypergeometricDIY.csv`

    - GSEApy is an inbuilt wrapper for Enrichr (allows list of human/mouse genes to compare against numerous biological libraries- pathways, diseases, genesets). This analysis is also provided as a benchmark using inbuilt python/enrichment analysis.
        - documentation: https://gseapy.readthedocs.io/en/latest/introduction.html
        - output: `~/enrichments_GSEA/{YOURGENE_KO}_{YOUR_POI}_GSEApy.csv`

## Plot WildType vs Mutant fractions
This program provides a pipeline to plot all proteins in each mixture individually or averaged for each module. The average of all protein variations for WT and MUT are highlighted in the plot, and individual protein changes are shown in a more transparent color to provide context.

Please navigate into `~/compare_MutvsWTfractions/WT_Mut_fractionPlot.ipynb` and change the variables in the second cell to your data. Please work through the notebook, changing the input to what you need specifically!

You can just run main, at the function definition there are calls (commented out). If you need to debug, uncomment for convenience.
