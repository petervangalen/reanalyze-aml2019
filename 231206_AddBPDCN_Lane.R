# Peter van Galen, 230830
# Integrate data from AML paper (2019) and BPDCN paper (2023)
# Links to papers: https://doi.org/10.1016/j.cell.2019.01.031, https://doi.org/10.1038/s41586-023-06156-8
# Note: there are unresolved errors below and the script is not very well polished

# Load required libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(janitor)
library(data.table)
library(ComplexHeatmap)
library(ggrastr)

# Start with a clean slate
rm(list=ls())

# Functions & colors
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))
popcol_tib <- read_excel("Single-cell_AML_colors.xlsx")
cell_colors_aml <- popcol_tib$mycol
names(cell_colors_aml) <- popcol_tib$pop
heat_colors <- read_excel("Single-cell_AML_colors.xlsx", sheet = 2, col_names = "heat")$heat

# LOAD AML PAPER DATA -----------------------------------------------------------------------------

# This first needs to be downloaded as described in README.md
aml <- readRDS("Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))

# Scale data (scaled values should be used for heatmap)
aml <- ScaleData(aml, scale.max = 4)


# LOAD BPDCN PAPER DATA -----------------------------------------------------------------------------
# Griffin et al, Nature 2023

# Load functions, colors, and data
# This is for Peter's local file system - the same files can be found on https://github.com/petervangalen/Single-cell_BPDCN
popcol.tib <- read_excel("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/Single-cell_BPDCN_colors.xlsx")
cell_colors_bpdcn <- popcol.tib$hex[1:21]
names(cell_colors_bpdcn) <- popcol.tib$pop[1:21]
sample_colors <- popcol.tib$hex[23:40]
names(sample_colors) <- popcol.tib$pop[23:40]
group_colors <- popcol.tib$hex[41:43]
names(group_colors) <- popcol.tib$pop[41:43]
mut_colors <- popcol.tib$hex[44:46]
names(mut_colors) <- popcol.tib$pop[44:46]
malignant_stage <- popcol.tib$hex[47:50]
names(malignant_stage) <- popcol.tib$pop[47:50]

# Load Seurat objects. These can be downloaded from https://vangalenlab.bwh.harvard.edu/resources/bpdcn-uv/
seurat_files <- list.files("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/04_XV-seq",
                           pattern = "*.rds", full.names = T)
bpdcn_ls <- lapply(seurat_files, function(x) readRDS(x))
names(bpdcn_ls) <- cutf(basename(seurat_files), d = "_")
bpdcn <- merge(bpdcn_ls[[1]], bpdcn_ls[2:length(bpdcn_ls)])
bpdcn$orig.ident2 <- ifelse(grepl("BM", bpdcn$orig.ident), yes = cutf(bpdcn$replicate, d = "\\."), no = bpdcn$orig.ident)

# Add metadata from 11.3_Classify_malignant_BPDCN.R
MalignantCalls_df <- read.table("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/11_pDC_expr/11.3_MalignantCalls_Final.txt", header = T, row.names = "cell")
all(rownames(bpdcn@meta.data) == rownames(MalignantCalls_df))
bpdcn$is_malignant <- MalignantCalls_df$is_malignant

# Load genotyping information
genotyping_tables.tib <- read_excel("~/DropboxMGB/Projects/Single-cell_BPDCN/AnalysisPeter/scBPDCN-analysis/04_XV-seq/XV-seq_overview.xlsx")
# Replace different MTAP entries with one, just as in 4.1_Add_GoT-XV_to_Seurat.R
genotyping_tables.tib$Mutation <- gsub("MTAP.rearr.*", "MTAP.rearr", genotyping_tables.tib$Mutation)

# Ensure nice ordering in visualizations
bpdcn$orig.ident <- factor(bpdcn$orig.ident, levels = c("BM", "Pt1Dx", "Pt1Rem", "Pt5Dx", "Pt9Dx", "Pt10Dx", "Pt10Rel",
                                                        "Pt12Dx", "Pt12Rel", "Pt14Dx", "Pt15Dx", "Pt16Dx"))
bpdcn$CellType <- factor(bpdcn$CellType, levels = names(cell_colors_bpdcn))

# UMAP of normal cell  types
p1 <- as_tibble(bpdcn@meta.data) %>% filter(orig.ident == "BM") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point_rast(size = 0.1) +
  scale_color_manual(values = cell_colors_bpdcn) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
p1

# Scale data (scaled values should be used for heatmap)
bpdcn <- ScaleData(bpdcn, scale.max = 4)


# MAKE HEATMAPS OF PIK3R5 AND SPI1 IN 2019 HEALTHY DONOR DATA -------------------------------------

# Create & wrangle tibble with metadata ---------
metadata_tib <- as_tibble(aml@meta.data, rownames = "cell")

# Filter for AML cells at Dx or BM cells. Also remove healthy cells from AML patients since their gene expression may be aberrant
metadata_filter <- metadata_tib %>% filter(grepl("AML.*D0", orig.ident) & grepl("-like", CellType) | grepl("BM", orig.ident)) %>%
  mutate(Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))

# Bar plot showing cell type proportions in AML patients
metadata_filter %>% filter(Donor == "AML") %>% select(orig.ident, CellType) %>% group_by(orig.ident, CellType) %>% count() %>% group_by(orig.ident) %>% mutate(prop = n/sum(n)) %>% ggplot(aes(x = orig.ident, y = prop, fill = CellType)) + geom_bar(stat = "identity") + scale_fill_manual(values = cell_colors_aml) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())

# Arrange
metadata_filter <- metadata_filter %>% arrange(Donor, CellType)

# Filter by genes and cell types ----------------

# Select genes to visualize & extract matrix with expression values
genes <- c("PIK3R5", "SPI1")
expr_mat <- as.matrix(LayerData(aml, layer = "scale.data"))[genes,]

# First, subset expression matrix. Then, convert the matrix to a data table & add metadata (make sure it joins by "cell")
plot_mat <- expr_mat[,metadata_filter$cell]
plot_mat[,1:3]
plot_dt <- as.data.table(as.table(plot_mat))
colnames(plot_dt) <- c("gene", "cell", "expr")
plot_dt <- left_join(plot_dt, metadata_filter)

# Split by healthy/AML
plot_healthy_dt <- filter(plot_dt, Donor == "Healthy")
plot_aml_dt <- filter(plot_dt, Donor == "AML") # USED IN NEXT SECTION
# Check
plot_healthy_dt$CellType %>% tabyl
plot_aml_dt$CellType %>% tabyl

# Calculate average expression per cell type; create matrices
dt2mat <- function(x) {
  # Summarize the data table to get the average for each group
  summarized_dt <- x[, .(avg_expr = mean(expr)), by = .(gene, CellType)]
  
  # Convert the summarized data table back to a matrix
  summarized_mat <- matrix(summarized_dt$avg_expr, 
                           nrow = length(unique(summarized_dt$gene)),
                           ncol = length(unique(summarized_dt$CellType)))
  rownames(summarized_mat) <- unique(summarized_dt$gene)
  colnames(summarized_mat) <- unique(summarized_dt$CellType)
  return(summarized_mat)
}

summarize_healthy_mat <- dt2mat(plot_healthy_dt)
summarize_aml_mat <- dt2mat(plot_aml_dt) # USED IN NEXT SECTION

# Create Heatmap object
hm1 <- Heatmap(summarize_healthy_mat[1,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)
hm2 <- Heatmap(summarize_healthy_mat[2,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)

draw(hm1 %v% hm2)


# MAKE HEATMAPS OF PIK3R5 AND SPI1 IN 2023 HEALTHY DONOR DATA -------------------------------------

# Create & wrangle tibble with metadata ---------
metadata_tib <- as_tibble(bpdcn@meta.data, rownames = "cell")

# Filter for cells from healthy donors
metadata_filter <- metadata_tib %>% filter(orig.ident == "BM" | orig.ident != "BM" & is_malignant == "Malignant" & CellType == "pDC") %>% arrange(CellType)

# Select genes to visualize & extract matrix with expression values
genes <- c("PIK3R5", "SPI1")
expr_mat <- as.matrix(LayerData(bpdcn, layer = "scale.data"))[genes,]

# First, subset expression matrix. Then, convert the matrix to a data table & add metadata (make sure it joins by "cell")
plot_mat <- expr_mat[,metadata_filter$cell]
plot_mat[,1:3]
plot_dt <- as.data.table(as.table(plot_mat))
colnames(plot_dt) <- c("gene", "cell", "expr")
plot_dt <- left_join(plot_dt, metadata_filter)

# Split by healthy/BPDCN
plot_healthy_dt <- filter(plot_dt, orig.ident == "BM")
plot_bpdcn_dt <- filter(plot_dt, orig.ident != "BM")
# Check
plot_healthy_dt$CellType %>% tabyl
plot_bpdcn_dt$CellType %>% tabyl

# Calculate average expression per cell type; create matrices
summarize_healthy_mat <- dt2mat(plot_healthy_dt)
summarize_bpdcn_mat <- dt2mat(plot_bpdcn_dt)

# Create Heatmap object
hm3 <- Heatmap(summarize_healthy_mat[1,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)
hm4 <- Heatmap(summarize_healthy_mat[2,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)

draw(hm3 %v% hm4)


# MAKE HEATMAP OF PIK3R5 AND SPI1 IN LEUKEMIA DATA ------------------------------------------------

# The following is questionable at best, since the two matrices were scaled differently
summarize_leukemia_mat <- cbind(summarize_aml_mat, summarize_bpdcn_mat)

hm5 <- Heatmap(summarize_leukemia_mat[1,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)
hm6 <- Heatmap(summarize_leukemia_mat[2,,drop=F],
               col = heat_colors,
               cluster_rows = F,
               cluster_columns = F,
               row_names_gp = gpar(fontsize = 10),
               show_column_names = T,
               column_title_gp = gpar(fontsize = 10),
               border = T)

draw(hm5 %v% hm6)


# 240124: Quick and dirty analysis to compare better (this failed miserably)

common_genes <- intersect(rownames(aml), rownames(bpdcn))

aml_geneSubset <- subset(aml, features = common_genes)
bpdcn_geneSubset <- subset(bpdcn, features = common_genes)

# Remove AML layers
aml_small1 <- DietSeurat(aml_geneSubset)
aml_small2 <- DietSeurat(aml_geneSubset, layers = "counts")
aml_small3 <- DietSeurat(aml_geneSubset, layers = list(counts = "counts"))

Layers(aml_small1) # none of this worked
Layers(aml_small2) # none of this worked

aml_geneSubset[["RNA"]]$data <- NULL # Ok finally
aml_geneSubset[["RNA"]]$scale.data <- NULL # Ok finally
Layers(aml_geneSubset) # Ok finally

# Remove BPDCN layers
bpdcn_small1 <- DietSeurat(bpdcn_geneSubset)
bpdcn_small2 <- DietSeurat(bpdcn_geneSubset, layers = "counts")
bpdcn_small3 <- DietSeurat(bpdcn_geneSubset, layers = list(counts = "counts"))

Layers(bpdcn_small1) # none of this worked
Layers(bpdcn_small2) # none of this worked

bpdcn_geneSubset[["RNA"]]$data <- NULL
bpdcn_geneSubset[["RNA"]]$scale.data <- NULL

Layers(bpdcn_geneSubset) # nope

LayerData(bpdcn_geneSubset, "data") <- NULL

Layers(bpdcn_geneSubset) # still nope

# Maybe this is the issue?
Version(aml)
Version(bpdcn)

# As you can read here, I am encountering numerous issues: https://github.com/satijalab/seurat/issues/8054#issuecomment-1907373226

# What I would like to do is merge, if I had been able to remove the data and scale.data layers:
all <- merge(aml_geneSubset, bpdcn_geneSubset)
