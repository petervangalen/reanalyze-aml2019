# Peter van Galen, 240719
# Create heatmap of translation initiation factors from single-cell AML paper (van Galen, Hovestadt et al, Cell 2019) for Supplemental Figure 1

# Load required libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(data.table)
library(ComplexHeatmap)

# Set working directory
# fmt: skip
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/Github/reanalyze-aml2019/2024_Lee")
# VM
setwd("/home/unix/vangalen/reanalyze-aml2019/2024_Lee")

# Start with a clean slate
rm(list = ls())

# Functions & colors
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}
popcol_tib <- read_excel("../Single-cell_AML_colors.xlsx")
cell_colors <- popcol_tib$mycol
names(cell_colors) <- popcol_tib$pop
heat_colors <- read_excel(
  "../Single-cell_AML_colors.xlsx",
  sheet = 2,
  col_names = "heat"
)$heat

# Load AML paper data -----------------------------------------------------------------------------

# Load AML paper data. This first needs to be downloaded as described in README.md
aml <- readRDS("../Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))

# Scale data (scaled values should be used for heatmap)
aml <- ScaleData(aml, scale.max = 4)


# Create & wrangle tibble with metadata -----------------------------------------------------------

# Extract metadata from Seurat object
metadata_tib <- as_tibble(aml@meta.data, rownames = "cell")

# Filter for AML cells at Dx or BM cells. Also remove healthy cells from AML patients since their gene expression may be aberrant
metadata_filter <- metadata_tib %>%
  filter(
    grepl("AML.*D0", orig.ident) &
      grepl("-like", CellType) |
      grepl("BM", orig.ident)
  ) %>%
  mutate(
    Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")
  ) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))

# Arrange
metadata_filter <- metadata_filter %>% arrange(Donor, CellType)


# Filter by genes and cell types ------------------------------------------------------------------

#fmt: skip
genes <- c("EIF1", "EIF1AX", "EIF1AY", "EIF2S1", "EIF2S2", "EIF2S3",
  "EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E", "EIF3F", "EIF3G", "EIF3H",
  "EIF3I", "EIF3J", "EIF3K", "EIF3CL", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4E",
  "EIF4E1B", "EIF4EBP1", "EIF4G1", "EIF4G2", "EIF4G3", "DDX3X", "EIF4B",
  "EIF4H", "EIF5", "EIF5B"
)
#fmt: skip
metadata_filter <- metadata_filter %>%
  filter(
    CellType %in%
      c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC", "HSC-like",
        "Prog-like", "GMP-like", "ProMono-like", "Mono-like", "cDC-like"
      )
  )
anno_fac <- metadata_filter %>%
  group_by(CellType, Donor) %>%
  summarize() %>%
  .$Donor


# Create expression matrix and summary expression matrix ------------------------------------------

# Extract expression matrix from Seurat object
expr_mat <- as.matrix(LayerData(aml, layer = "scale.data"))[genes, ]

# Create expression matrix for Heatmap. First, convert the matrix to a data table & add metadata (make sure it joins by "cell")
plot_mat <- expr_mat[, metadata_filter$cell]
plot_dt <- as.data.table(as.table(plot_mat))
colnames(plot_dt) <- c("gene", "cell", "expr")
plot_dt <- left_join(plot_dt, metadata_filter)

# Summarize the data table to get the average for each group
summarized_dt <- plot_dt[, .(avg_expr = mean(expr)), by = .(gene, CellType)]

# Convert the summarized data table back to a matrix
summarized_mat <- matrix(
  summarized_dt$avg_expr,
  nrow = length(unique(summarized_dt$gene)),
  ncol = length(unique(summarized_dt$CellType))
)
rownames(summarized_mat) <- unique(summarized_dt$gene)
colnames(summarized_mat) <- unique(summarized_dt$CellType)
summarized_mat


# Plot Heatmap (every column is a cell type) ------------------------------------------------------

# Define annotation objects
top_anno.ha <- HeatmapAnnotation(
  Split = anno_fac,
  annotation_name_gp = gpar(fontsize = 10),
  border = T
)

# Order by expression in primitive cells
ordered_genes <- names(sort(rowMeans(summarized_mat[, 7:9]), decreasing = T))
summarized_mat <- summarized_mat[ordered_genes, ]

# Create Heatmap object
hm2 <- Heatmap(
  summarized_mat,
  col = heat_colors,
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 10),
  show_column_names = T,
  column_split = anno_fac,
  top_annotation = top_anno.ha,
  name = "Expr",
  column_title_gp = gpar(fontsize = 10),
  border = T
)

pdf("240719_eIF4F-heatmap", width = 6, height = 6)
print(hm2)
dev.off()
