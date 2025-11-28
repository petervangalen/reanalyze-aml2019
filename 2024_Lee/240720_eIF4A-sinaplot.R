# Peter van Galen, 240720
# Generate plots for translation initation factor genes from the single-cell AML data

# Load required libraries
library(tidyverse)
library(Seurat)
library(ggpubr)
library(ggforce)

# Start with a clean slate
rm(list = ls())

# Frequently used function
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Set working directory
# fmt: skip
setwd("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/Github/reanalyze-aml2019/2024_Lee")
# VM
setwd("/home/unix/vangalen/reanalyze-aml2019/2024_Lee")

# Load AML paper data. This first needs to be downloaded as described in README.md
aml <- readRDS("../Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))

# Generate a metadata tibble and add gene expression (choose one gene)
metadata <- as_tibble(aml@meta.data, rownames = "cell")
mygene <- "EIF4A1"
metadata$mygene <- LayerData(aml, layer = "data")[mygene, ]

# Filter for AML cells at Dx or BM cells. Also remove "healthy" cells from AML patients since their gene expression may be aberrant
metadata.filter <- metadata %>%
  filter(
    grepl("AML.*D0", orig.ident) &
      grepl("-like", CellType) |
      grepl("BM", orig.ident)
  ) %>%
  mutate(
    Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")
  ) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))


# Three violins (240719) ------------------------------------------------------

# Calculate P-values for all genes, then extract EIF4A1
HSC_Pvals <- FindMarkers(
  aml,
  ident.1 = filter(metadata.filter, CellType == "HSC")$cell,
  ident.2 = filter(metadata.filter, CellType == "HSC-like")$cell,
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.feature = 0,
  return.thresh = 1.01
)
Prog_Pvals <- FindMarkers(
  aml,
  ident.1 = filter(metadata.filter, CellType == "Prog")$cell,
  ident.2 = filter(metadata.filter, CellType == "Prog-like")$cell,
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.feature = 0,
  return.thresh = 1.01
)
GMP_Pvals <- FindMarkers(
  aml,
  ident.1 = filter(metadata.filter, CellType == "GMP")$cell,
  ident.2 = filter(metadata.filter, CellType == "GMP-like")$cell,
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.feature = 0,
  return.thresh = 1.01
)
#fmt: skip
plot_P <- tribble(~x, ~y, ~P,
  1.5, 3.5, HSC_Pvals["EIF4A1", ]$p_val_adj,
  3.5, 3.5, Prog_Pvals["EIF4A1", ]$p_val_adj,
  5.5, 3.5, GMP_Pvals["EIF4A1", ]$p_val_adj
)

# The sina/violin plot shows expression in every individual cell (symbol)
p1 <- metadata.filter %>%
  filter(
    CellType %in% c("HSC", "Prog", "GMP", "HSC-like", "Prog-like", "GMP-like")
  ) %>%
  mutate(
    CellType = factor(
      CellType,
      levels = c("HSC", "HSC-like", "Prog", "Prog-like", "GMP", "GMP-like")
    )
  ) %>%
  ggplot(aes(x = CellType, y = mygene, color = Donor)) +
  geom_sina(scale = "width", size = 0.2) +
  geom_violin(
    scale = "width",
    draw_quantiles = 0.5, # line location changes with ggplot 4.0.0
    color = "black",
    fill = NA
  ) +
  geom_text(
    data = plot_P,
    aes(x = x, y = y, label = paste("p = ", scales::scientific(P, digits = 2))),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    color = "black"
  ) +
  ylab("Normalized expression, log(TP10K+1)") +
  ggtitle(mygene) +
  scale_x_discrete(drop = F) +
  scale_color_manual(values = c(Healthy = "#66cdaa", AML = "#cd5c5c")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid = element_blank()
  )

# Visualize the plot
pdf("240720_eIF4A-sinaplot.pdf", width = 4, height = 4)
p1
dev.off()

# Look at P-values for other eIF4F complex members
HSC_Pvals[c("EIF4A1", "EIF4A2", "EIF4E", "EIF4G1"), ]
Prog_Pvals[c("EIF4A1", "EIF4A2", "EIF4E", "EIF4G1"), ]
GMP_Pvals[c("EIF4A1", "EIF4A2", "EIF4E", "EIF4G1"), ]


# Merge violins to simplify figure (250927) -----------------------------------

# Calculate P-values for all genes, then extract EIF4A1
HSPC_Pvals <- FindMarkers(
  aml,
  ident.1 = filter(metadata.filter, CellType %in% c("HSC", "Prog"))$cell,
  ident.2 = filter(
    metadata.filter,
    CellType %in% c("HSC-like", "Prog-like")
  )$cell,
  logfc.threshold = 0,
  min.pct = 0,
  min.cells.feature = 0,
  return.thresh = 1.01
)

# Plot
metadata.filter %>%
  filter(
    CellType %in% c("HSC", "Prog", "HSC-like", "Prog-like")
  ) %>%
  ggplot(aes(x = Donor, y = mygene, color = Donor)) +
  geom_sina(scale = "width", size = 0.1) +
  geom_violin(
    scale = "width",
    color = "grey",
    fill = NA
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    color = "black",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 1.4,
    y = 3.5,
    label = paste(
      "p =",
      scales::scientific(HSPC_Pvals["EIF4A1", ]$p_val_adj, digits = 2)
    ),
    size = 5
  ) +
  ylab("Normalized expression, log(TP10K+1)") +
  ggtitle(paste0(mygene, " in healthy vs. AML HSPCs")) +
  scale_x_discrete(drop = F) +
  scale_color_manual(values = c(Healthy = "#66cdaa", AML = "#FF0066")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave("240720_eIF4A-sinaplot-merge.pdf", height = 5, width = 5)

# This last figure is not actually used in the paper
