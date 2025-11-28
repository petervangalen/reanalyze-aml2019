# Peter van Galen, 230830
# Create gene expression heatmap from single-cell AML paper (van Galen, Hovestadt et al, Cell 2019)

# Load required libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(data.table)
library(ComplexHeatmap)

# Set working directory
# fmt: skip
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/Github/reanalyze-aml2019/2024_Talvard-Balland")

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

# This first needs to be downloaded as described in README.md
aml <- readRDS("../Seurat_AML.rds")

# Scale data (scaled values should be used for heatmap)
aml <- ScaleData(aml, scale.max = 4)


# Create & wrangle tibble with metadata -----------------------------------------------------------

# Extract metadata from Seurat object
metadata_tib <- as_tibble(aml@meta.data, rownames = "cell")

# Filter for AML cells at Dx or BM cells
metadata_filter <- metadata_tib %>%
  filter(grepl("AML.*D0", orig.ident) | grepl("BM", orig.ident)) %>%
  mutate(
    Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")
  ) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))

# Subset for cell types of interest
# fmt: skip
metadata_filter <- metadata_filter %>% mutate(CellType_Donor = paste0(CellType, "_", Donor)) %>% filter(CellType_Donor %in%
  c("HSC_Healthy", "Prog_Healthy", "GMP_Healthy", "ProMono_Healthy", "Mono_Healthy", "cDC_Healthy",
    "pDC_Healthy", "earlyEry_Healthy", "lateEry_Healthy", "ProB_Healthy", "B_Healthy", "Plasma_Healthy", "T_Healthy", "CTL_Healthy", "NK_Healthy",
    "HSC-like_AML", "Prog-like_AML", "GMP-like_AML", "ProMono-like_AML", "Mono-like_AML", "cDC-like_AML",
    "pDC_AML", "earlyEry_AML", "lateEry_AML", "ProB_AML", "B_AML", "Plasma_AML", "T_AML", "CTL_AML", "NK_AML")) %>%
  mutate(CellType_Donor = factor(CellType_Donor, levels =
  c("HSC_Healthy", "Prog_Healthy", "GMP_Healthy", "ProMono_Healthy", "Mono_Healthy", "cDC_Healthy",
    "pDC_Healthy", "earlyEry_Healthy", "lateEry_Healthy", "ProB_Healthy", "B_Healthy", "Plasma_Healthy", "T_Healthy", "CTL_Healthy", "NK_Healthy",
    "HSC-like_AML", "Prog-like_AML", "GMP-like_AML", "ProMono-like_AML", "Mono-like_AML", "cDC-like_AML",
    "pDC_AML", "earlyEry_AML", "lateEry_AML", "ProB_AML", "B_AML", "Plasma_AML", "T_AML", "CTL_AML", "NK_AML")))

metadata_filter <- metadata_filter %>% arrange(CellType_Donor)

# Filter by genes and cell types ------------------------------------------------------------------

# Select genes to visualize & extract matrix with expression values
# We will not include CEACAM1 because it is only detected in 0.4% of cells:
#mean(LayerData(aml, layer = "counts")["CEACAM1",])*100
genes <- c("HAVCR2", "LGALS9", "HMGB1")
expr_mat <- as.matrix(LayerData(aml, layer = "scale.data"))[genes, ]

# Create expression matrix for Heatmap
plot_mat <- expr_mat[, metadata_filter$cell]
plot_mat[, 1:3]

# Define annotation objects
top_anno.ha <- HeatmapAnnotation(
  CellType = metadata_filter$CellType,
  Donor = metadata_filter$Donor,
  col = list(
    CellType = cell_colors[unique(as.character(metadata_filter$CellType))]
  ),
  annotation_name_gp = gpar(fontsize = 10),
  border = T
)

# Create Heatmap object
hm1 <- Heatmap(
  plot_mat,
  col = heat_colors[3:11],
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 10),
  show_column_names = F,
  column_split = metadata_filter$Donor,
  top_annotation = top_anno.ha,
  name = "Expr",
  column_title_gp = gpar(fontsize = 10),
  border = T,
  raster_quality = 10
)

pdf("230830_Zeiser_Heatmap.pdf", width = 10, height = 2)
print(hm1)
dev.off()

# To annotate in Illustrator
metadata_filter$Donor %>% tabyl
# P-values
plot_mat[, 1:3]
t.test(
  x = plot_mat["HAVCR2", grepl("BM", colnames(plot_mat))],
  y = plot_mat["HAVCR2", !grepl("BM", colnames(plot_mat))]
)
t.test(
  x = plot_mat["LGALS9", grepl("BM", colnames(plot_mat))],
  y = plot_mat["LGALS9", !grepl("BM", colnames(plot_mat))]
)
t.test(
  x = plot_mat["HMGB1", grepl("BM", colnames(plot_mat))],
  y = plot_mat["HMGB1", !grepl("BM", colnames(plot_mat))]
)
# Out of curiosity
as_tibble(t(plot_mat), rownames = "cell") %>%
  left_join(metadata_filter) %>%
  group_by(CellType_Donor) %>%
  summarize(mean_TIM3 = mean(HAVCR2)) %>%
  print(n = Inf)
