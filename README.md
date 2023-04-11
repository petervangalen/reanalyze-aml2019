# Reanalyze AML 2019

This repository contains scripts to reanalyze the single-cell RNA-seq data on AML patient samples that was originally published in 2019 ([van Galen et al, Cell 2019](https://pubmed.ncbi.nlm.nih.gov/30827681/)). The analyses here mostly pertain to assessing gene expression in normal and malignant bone marrow cell types.

With the original publication, the raw data was deposited on GEO (accession [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256)); you don't need this for any of the analysis described here. The code by Volker Hovestadt is found on the Bernstein lab [aml2019](https://github.com/BernsteinLab/aml2019) repository. 

First, you need to download the data to the top level folder of this repository. To do this, navigate to this folder in the terminal, and run the following command:
`wget https://www.dropbox.com/s/399x045zc57fiut/Seurat_AML.rds`

Then, you can run `Visualize_expression.R` to look at the expression of genes. For example, you can generate these plots, or similar for any other gene by changing the variable `mygene`:
![alt text](/images/HOXA9_bar_and_sina.png "HOXA9 expression")
The bar plots show the mean expression across cell types in healthy donors (green) and malignant cells in AML patients at diagnossis (red). The sina/violin plots show expression in every individual cell (symbol).

You can also split the analysis by donors or any other variable in aml@meta.data
