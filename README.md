# SpatialProteomics

Taken from source [SWIP-Proteomics by T. Wesley & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file)

Please refer to the SWIP-proteomics github for the foundation of this program. Below are the specific details on how to get variation analysis for each protein individually, and then protein modules.

Refer to the following paper for methods and further conceptuals of the pipeline. [Genetic disruption of WASHC4 drives endo-lysosomal dysfunction and cognitive-movement impairments in mice and humans](https://elifesciences.org/articles/61590) 


## Environment
### R:
This project uses `renv` for dependency management. To set up the R environment for this project:
1. Open R with the current repo (SpatialProteomics)
    - if you are in bash, type `R` to enter R console
2. Run `install.packages("renv")`
3. Run `renv::restore()`. This restores the R environment from renv.lock in 'analysis'.

If you are missing a package:
`install.packages("PKG")`

It may be a package from BiocManger specifically. In that event:
`BiocManager::install("AnnotationDbi")`


If you cannot find the package, it is likely a github repository. Please check soderling-lab repositories, the package will likely be in there.

`devtools::install_github("{github-user}/{githubrepo}")`

You can load these packages using `library(PKG)`

### Python

For running 'analysis/3_Clustering/2_leidenalg-clustering.py' you will need a specific python environment. 

At the home directory of this repository 'SpatialProteomics', please run:

`conda env create -f analysis/environment.yml` 

and activate:
`conda activate spatial_env`

In the event you run into missing packages, please install with `conda install {pkg}` or `pip install {pkg}`

**NOTE: For the singular python file in the module analysis pipeline. You must activate the environment, and then run it from the command line.**

`python analysis/3_Clustering/2_leidenalg-clustering.py`

## Necessary Pre-processing
Given the Tandem Mass Taggged (TMT) Mass Spectrometry data from the Proteomics Core, please use the Normalized data sheet (labeled in xlsx file) moving further.

Please run '0_prepData_transform.py' from 'analysis/2_SWIP-TMT' or 'analysis/3_Clustering' with your target dataset, to prep the data to start at '1_generate-network.R', '1_MSstatsTMT-analysis.R'

This will transform the data to a long format so the fractions and mixtures are structured vertically rather than horizontally.

**In the event you are starting from the PSMs, please run 0_PD-data-preprocess.R** 
Please refer to the source code [SWIP-Proteomics by T. Wesley & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file) if you need to start here. 

### MSstatsTMT
### Assess differential protein abundance for intrafraction comparisons between WT and MUT

Navigate into and run: 
- analysis/2_SWIP-TMT/1_MSstatsTMT-analysis.R
    - output: adjacency matrix, neten adjacency matrix, TMT protein data
     - You can only run this if you have the PSMs. Reach out the proteomics core in the event you need this. Some examples are under ~/PSM

- analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R
    - Start with this to get the entire pipeline below for protein clustering
    - saved files will be printed in console

### Protein clustering into modules using Leiden Algorithm

Please ensure you have run ~/analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R

Navigate into analysis/3_Clustering and run:
- analysis/3_Clustering/1_generate-network.R
- analysis/3_Clustering/2_leidenalg-clustering.py
    - please activate the python environment using environment.yml file and run this from the command line after making your edits. `python 2_leidenalg-clustering.py`
- analysis/3_Clustering/3_post-leidenalg.R
- analysis/3_Clustering/4_generate-PPI-network.R

Given the output from the four above programs, you may now find out which modules are significant in that WT and MUT protein abundances vary within the module.

### Identify significant modules with Linear Mixed Models
Navigate into analysis/4_Module-Analysis and run:
- analysis/4_Module-Analysis/1_module-lmerTest-analysis.R

You have found the clustered modules in your dataset! Please navigate into ~/tables, the protein modules are saved under '{YOUR_GENE}_TMT-Module-Results.xlsx'
 - Each protein is assigned to a module. You can also analyze how each protein does individually using MSStatsTMT procedure detailed above. MSStatsTMT will yield how much the abundance of each protein changed as a result of the Gene Of Interest (GOI) being knocked out.

### Gene Ontologies
If you want to look at Gene Ontologies, you may run the following file the original author used for WASH proteins:
- analysis/4_Module-Analysis/2_module-GSEA.R

**OR**

You can scrape ShinyGo using '~/ShinyGo/ShinyGO_analysis/1_analyze_modules.ipynb'. This works for in an automated loop but sometimes the websites is spontaneous so in intervals it will er.. 
- ie. if you have 50 modules, it may err on the 25th module.
You can create a heatmap with the identified pathways and each module using 'ShinyGo/ShinyGO_analysis/2_cluster_analysis.py'

**OR** if you have a specific pathway in mind..

Navigate into ~/geneontologies and run:
 - geneontologies/get_KEGGnums.ipynb 
    - input: ~/tables/{KOGENE}-TMT-Module-Results.xlsx
 - geneontologies/get_pathways.ipynb
    - input: ~/tables/{KOGENE}_moduleResults_KEGGortholog.csv
    - This will calculate gene_ontologies with your pathway of interest (POI) using a hypergeometric test. Each module is analyzed at a time, and the overlap between genes in that module vs genes in POI will indicate p-value. (p-val <0.05 == module is enriched in that pathway) 
        - P-value is adjusted with benjamin-hochberg for multiple testing (FDR) since multiple memberships are analyzed at once.
        - output: ~/enrichments_DIY/{YOURGENE_KO}_{YOUR_POI}_hypergeometricDIY.csv

    - GSEApy is an inbuilt wrapper for Enrichr (allows list of human/mouse genes to compare against numerous biological libraries- pathways, diseases, genesets). This analysis is also provided as a benchmark using inbuilt python/enrichment analysis.
        - documentation: https://gseapy.readthedocs.io/en/latest/introduction.html
        - output: ~/enrichments_GSEA/{YOURGENE_KO}_{YOUR_POI}_GSEApy.csv
