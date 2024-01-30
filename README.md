# SpatialProteomics

Taken from source [SWIP-Proteomics by T. WA Bradshaw & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file)

Please refer to the soderling-lab/SwipProteomics github for the foundation of this program. Below are the specific details on how to get variation analysis for each protein individually, and then protein modules.

All programming in `analysis/` sourced from Tyler WA Bradshaw. Adapted in certain areas to fit current needs/updates.

Refer to the following paper for methods and further conceptuals of the pipeline. [Genetic disruption of WASHC4 drives endo-lysosomal dysfunction and cognitive-movement impairments in mice and humans](https://elifesciences.org/articles/61590)


## Pipeline Overview
Tandem Mass Tag (TMT) is a labeling technique used in mass spectrometry for protein quantification and identification in large-scale proteomic studies. In this study, we utilized TMT to analyze fractionated samples from WildType (WT) and transgenic disease models (MUT).

Biological samples are taken from both WT and MUT models, and each is centrifuged into seven distinct cellular fractions.  Each cellular fraction is labeled with a TMT barcode, differentiating it from the others. A plasma control is taken alongside each biological sample $$(7 \text{ WT TMT-fractions} + 7 \text{ MUT TMT-fractions} + 2 {\text{ Control sample}})$$ per replicate.
One replicate has 16 TMT-fractions representing the WT, MUT mouse models. This experiment is done three times (triplicate) with different biological samples taken every time, yielding three mixtures. The triplicate is significant to diminish any bias or contamination when coming to analysis.

The fractions are processed through a mass spectrometer, and the proteome of all three mixtures is analyzed. The control fractions in each mixture are excluded from analysis. Each protein is correlated against itself across all fractions and mixtures for both WT and MUT models.

The resulting data forms an adjacency matrix of the protein network, which is then [network enhanced](https://github.com/soderling-lab/neten) to improve the signal-to-noise ratio [[1]](#references). Given the adjaceny matrix, a simplistic network is constructed. We then implement the [leiden algorithm](https://github.com/soderling-lab/leiden) to further cluster with intra-module and inter-module analysis, optimizing the Surprise function for community detection. [[1]](#references)   [[2]](#references). During module identification, variations between WT and MUT conditions are not considered. Post clustering, these variations are reintegrated, and we fit a linear mixed model to assess contrasts in protein abundance between WT and MUT. A module is considered significant if the contrast indicates that the transgenic condition substantially affects protein expression levels within that module, determined by a $\text{p-value} \leq 0.05$.
Each module, along with differences in protein expression between WT and MUT conditions, is visually represented. Proteins are plotted individually, and their averages are used to highlight distinctions between WT and MUT expressions.

Following module detection, we explore gene ontologies to identify enriched pathways using tools like ShinyGo [[4]](#references), GSEApy [[5]](#references), and hypergeometric tests. Additionally, for specific pathways (e.g., Parkinson's Disease), we determine which modules have proteins closely related to the relevant gene set.

**Example data files are provided under example_data/ for all steps. Please refer to this to get a sense of the program**

## Environments
### R:

This project uses `renv` for dependency management. To set up the R environment for this project:
1. Open R with the current repo (SpatialProteomics)
    - if you are in bash, type `R` to enter R console
2. Set the working directory to where `SpatialProteomics/` was downloaded. `setwd('path/to/SpatialProteomics')`
2. Run `install.packages("renv")`
3. Run `renv::restore(lockfile = "analysis/renv.lock")`. This restores the R environment from [renv.lock](analysis/renv.lock) in `analysis/`.

If you are missing a package:
`install.packages("PKG")`

It may be a package from [BiocManger](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) specifically. In that event:
`BiocManager::install("AnnotationDbi")`


If you cannot find the package, it is likely a github repository. Please check soderling-lab repositories, the package will likely be in there.

`devtools::install_github("{github-user}/{githubrepo}")`

You can load these packages using `library(PKG)`

### Python

**NOTE:** For notebooks (`.ipynb` files), you have to push the environment to the jupyter kernel so you can select it as a kernel. Please run 
`python -m ipykernel install --user --name <env> --display-name "displayname"`

<u> Spatial Proteomics Environment</u>

To run [`analysis/3_Clustering/2_leidenalg-clustering.py`](analysis/3_Clustering/2_leidenalg-clustering.py) you will need a specific python environment. 

In the `~/analysis/` of this repository `SpatialProteomics/`, please run:

`conda env create -f`[` environment_spatial.yml`](analysis/environment_spatial.yml)

and activate:
`conda activate spatial_env`

In the event you run into missing packages, please install with `conda install {pkg}` or `pip install {pkg}`

**NOTE: For the singular python file in the module analysis pipeline. You must activate the environment, and then run it from the command line:** `python analysis/3_Clustering/2_leidenalg-clustering.py`

<u> Gene Ontology Environment </u>

When you start with the gene ontology pipeline (after gathering all the module data) you must deactivate the `spatial_env` environment `conda deactivate`.

In the [~/geneontologies/](geneontologies) dir please run:        
`conda env create -f`[`environment_GO.yml`](geneontologies/environment_GO.yml)

and activate:
`conda activate GOenv`

#### ShinyGO environment
If you decide to webscrape with ShinyGo to get pathways for each module. Please deactivate whatever env is active `conda deactivate` and run

 `conda env create -f `[`environment_ShinyGo.yml`](ShinyGo/environment_ShinyGo.yml)

#### compare_MutvsWTfractions environment
If you would like to create plots comparing mutant vs wildtype fractions, then deactivate any current conda environment `conda deactivate` then run `conda env create -f ` [`environment_plot.yml`](compare_MutvsWTfractions/environment_plot.yml)

## Program Instructions

### Necessary Pre-processing
Given the Tandem Mass Taggged (TMT) Mass Spectrometry data from the Proteomics Core, please use the Normalized data sheet (labeled in excel file) moving further.

1. Please run your target dataset through [0_prepData_transform.py](analysis/2_SWIP-TMT/0_PD-data-preprocess.R) from `analysis/2_SWIP-TMT/` or `analysis/3_Clustering/` to prep the data to start at `1_generate-network.R` or `1_MSstatsTMT-analysis.R`
    - In  `0_prepData_transform.py`, navigate into the `main()` function and change the filename and worksheet accordingly. This will transform the data to a long format so the fractions and mixtures are structured vertically rather than horizontally.

**In the event you are starting from the PSMs, please run 0_PD-data-preprocess.R** 
Please refer to the source code [SWIP-Proteomics by T. Wesley & S. H. Soderling](https://github.com/soderling-lab/SwipProteomics?tab=readme-ov-file) if you need to start here. 

### MSstatsTMT
**Assess differential protein abundance between fractions between WT and MUT models**

Navigate into and run:
- [analysis/2_SWIP-TMT/0_PD-data-preprocess.R](analysis/2_SWIP-TMT/0_PD-data-preprocess.R)
    - input: PSM report, PSM samples 
    - You can only run this if you have the PSMs. Reach out the proteomics core in the event you need this. Some examples are under `~/PSM`
    - Variables and Gene strings must be changed to your case

- [analysis/2_SWIP-TMT/1_MSstatsTMT-analysis](/analysis/2_SWIP-TMT/1_MSstatsTMT-analysis.R)
    - input: data tables from `0_PD-data-preprocess.R`
    - output: adjacency matrix, neten adjacency matrix, TMT protein data

### Leiden Algorithm
**Spatial Protein clustering in terms of abundance regulation**

First to normalize the data, run [analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R](analysis/2_SWIP-TMT/2_Swip-TMT-normalization.R) with the long, transformed data structure.
- Start with this to get the entire pipeline below for protein clustering
- Saved files will be printed in console

Navigate into [analysis/3_Clustering](analysis/3_Clustering) to change variables and run:
- [1_generate-network.R](analysis/3_Clustering/1_generate-network.R)
- [2_leidenalg-clustering.py](analysis/3_Clustering/2_leidenalg-clustering.py)
    - please activate the python environment using environment.yml file and run this from the command line after making your edits. `python 2_leidenalg-clustering.py`
- [3_post-leidenalg.R](analysis/3_Clustering/3_post-leidenalg.R)
- [4_generate-PPI-network.R](analysis/3_Clustering/4_generate-PPI-network.R)

Given the output from the four above programs, you may now find out which modules are significant when MUT model introduces different protein expression levels. 

### Identify significant clusters (modules) with Linear Mixed Models
Navigate into [analysis/4_Module-Analysis](analysis/4_Module-Analysis) and run:
- [1_module-lmerTest-analysis.R](analysis/4_Module-Analysis/1_module-lmerTest-analysis.R)

You have found the clustered modules in your dataset! Please navigate into `~/tables`, the protein modules are saved under `{YOUR_GENE}_TMT-Module-Results.xlsx`
 - Each protein is assigned to a module. You can also analyze how each protein does individually using MSStatsTMT procedure detailed above. MSStatsTMT will yield how much the abundance of **each** protein changed relative to WT and MUT variances.

### Gene Ontologies
If you want to remain on the original pipeline and look for gene ontologies,
you may run the following file the original author used for WASH proteins (changes needed for your GOI):
- [analysis/4_Module-Analysis/2_module-GSEA.R](analysis/4_Module-Analysis/2_module-GSEA.R)

**OR**

You can scrape [ShinyGo](http://bioinformatics.sdstate.edu/go/) using [~/ShinyGo/ShinyGO_analysis/module_pathways_1.ipynb](ShinyGo/ShinyGO_analysis/module_pathways_1.ipynb). This works in an automated loop but sometimes the website is spontaneous so in intervals it will er.. 
- ie. if you have 50 modules, it may err on the 25th module, but all until 25 will save.

You can create a heatmap with the identified pathways modules using [ShinyGo/ShinyGO_analysis/cluster_analysis_2.py](ShinyGo/ShinyGO_analysis/cluster_analysis_2.py)

**OR** if you have a specific pathway in mind (Parkinson's, Autism..)

Navigate into [~/geneontologies](geneontologies/get_pathways.ipynb) and run:
 - [get_KEGGnums.ipynb](geneontologies/get_KEGGnums.ipynb)
    - input: `~/tables/{KOGENE}-TMT-Module-Results.xlsx`
    - gets KEGG ortholog nums to get pathways in next file

 - [get_pathways.ipynb](geneontologies/get_KEGGnums.ipynb)
    - input: `~/tables/{KOGENE}_moduleResults_KEGGortholog.csv`
    - Change the Pathway of Interest (POI) from Parkinson's to what you are looking for!
    - This will calculate gene_ontologies with your POI using a hypergeometric test. Each module is analyzed at a time, and the overlap between genes in that module vs genes in POI will indicate p-value. (p-val < 0.05 signifies the module is enriched in POI) 
        - P-value is adjusted with benjamin-hochberg for multiple testing (FDR) since multiple memberships are analyzed at once.
        - output: `~/enrichments_DIY/{YOURGENE_KO}_{YOUR_POI}_hypergeometricDIY.csv`

    - [GSEApy](https://gseapy.readthedocs.io/en/latest/introduction.html) is an inbuilt wrapper for Enrichr (allows list of human/mouse genes to compare against numerous biological libraries- pathways, diseases, genesets). This analysis is also provided as a benchmark using inbuilt python/enrichment analysis.
        - documentation: https://gseapy.readthedocs.io/en/latest/introduction.html
        - output: `~/enrichments_GSEA/{YOURGENE_KO}_{YOUR_POI}_GSEApy.csv`

## Plot WildType vs Mutant fractions
This program provides a pipeline to plot all proteins in each mixture individually or averaged for each module. The average of all protein variations for WT and MUT are highlighted in the plot, and individual protein changes are shown in a more transparent color to provide context.

Please navigate into [~/compare_MutvsWTfractions/WT_Mut_fractionPlot.ipynb](compare_MutvsWTfractions/WT_Mut_fractionPlot.ipynb) and change the variables in the second cell to your data. Please work through the notebook, changing the input to what you need specifically!

You can just run main, at the function definition there are calls (commented out). If you need to debug, uncomment for convenience.


## References
[1] Genetic Disruption of WASHC4 Drives Endo-lysosomal Dysfunction and Cognitive-Movement Impairments in Mice and Humans.
Courtland J.L., Bradshaw T.W.A., Waitt G., Soderblom E., Ho T., Rajab A., Vancini R., Kim I.H., Soderling S.H. (2021). eLife; 10:e61590 [doi: 10.7554/eLife.61590](https://elifesciences.org/articles/61590)

[2] From Louvain to Leiden: guaranteeing well-connected communities.
Traag, V.A., Waltman. L., Van Eck, N.-J. (2018). Scientific reports, 9(1), 5233. [10.1038/s41598-019-41695-z](https://www.nature.com/articles/s41598-019-41695-z)

[3] Network Enhancement as a general method to denoise weighted biological networks.
Wang B., Pourshafeie A., Zitnik M., Zhu J., Bustamante C.D., Batzoglou S., Leskovec J. (2018). Nature Communications, 9, 3108. [10.1038/s41467-018-05469-x](https://www.nature.com/articles/s41467-018-05469-x)

[4] Steven Xijin Ge, Dongmin Jung, Runan Yao, ShinyGO: a graphical gene-set enrichment tool for animals and plants, Bioinformatics, Volume 36, Issue 8, April 2020, Pages 2628–2629, https://doi.org/10.1093/bioinformatics/btz931

[5] Zhuoqing Fang, Xinyuan Liu, Gary Peltz, GSEApy: a comprehensive package for performing gene set enrichment analysis in Python,
Bioinformatics, 2022;, btac757, https://doi.org/10.1093/bioinformatics/btac757