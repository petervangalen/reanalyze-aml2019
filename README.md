# Reanalyze AML 2019

This repository contains scripts to analyze gene expression in healthy and malignant cells from the single-cell RNA-seq data on AML patient bone marrow samples that was published in 2019 ([van Galen, Hovestadt _et al._, Cell 2019](https://pubmed.ncbi.nlm.nih.gov/30827681/)).

### Set up

First, you need to [clone or download](https://www.freecodecamp.org/news/10-important-git-commands-that-every-developer-should-know/) this repository to your local disk. Then, you need to download the gene expression data to its top level folder (it is too large to share on Github). The gene expression data is saved as a Seurat object (Seurat_AML.rds), which is available from Figshare:

[https://doi.org/10.6084/m9.figshare.30581066.v1](https://doi.org/10.6084/m9.figshare.30581066.v1)

### Bar and sina plots

Use the script `230411_Bar_and_sina_plots.R` to look at the expression of genes. For example, you can generate these plots, or similar for any other gene by changing the variable `mygene`:

![alt text](/images/HOXA9_bar_and_sina.png "HOXA9 expression")

The bar plot shows the mean expression across cell types in healthy donors (green) and malignant cells in AML patients at diagnosis (red). The sina/violin plot shows expression in every individual cell (symbols).

You can also split the analysis by donors or any other variable in `aml@meta.data`.

### Heatmaps

Use the script `230824_Heatmaps.R` to visualize expression of several genes at once. This script generates heatmaps showing expression of genes (rows) in every cell (columns), or summary heatmaps showing average expression of genes (rows) in cell types (columns). For example:

![alt text](/images/translation_factors.png "Translation initiation")

### Folders

Folders contain scripts to reproduce specific published visualizations:
* 2024_Talvard-Balland → Talvard-Balland _et al._, J Clin Invest, 2024. [PMID: 38916965](https://pubmed.ncbi.nlm.nih.gov/38916965/)
* 2024_Lee → Lee _et al_., bioRxiv, 2024. [doi: https://doi.org/10.1101/2024.12.20.629809](https://doi.org/10.1101/2024.12.20.629809)

**Other resources:** With the original publication, the raw data was deposited on GEO (accession [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256)); you don't need the raw data for the analysis described here. The original code for cell type classification etc. by Volker Hovestadt is found on the Bernstein lab [aml2019](https://github.com/BernsteinLab/aml2019) repository. 
