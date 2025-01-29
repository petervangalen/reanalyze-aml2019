# Peter van Galen, 230411
# Check expression of genes of interest in data from single-cell AML paper (van Galen, Hovestadt et al, Cell 2019)

# Load required libraries
library(tidyverse)
library(Seurat)
library(ggforce)
library(cowplot)

# Start with a clean slate
rm(list=ls())

# Frequently used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load AML paper data. This first needs to be downloaded as described in README.md
aml <- readRDS("Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))


# Plot single gene --------------------------------------------------------------------------------

# Generate a metadata tibble and add gene expression (choose one gene)
metadata <- as_tibble(aml@meta.data, rownames = "cell")
mygene <- "HOXA9"
mygene <- "HAVCR2"
mygene <- "LGALS9"
mygene <- "CEACAM1"
mygene <- "HMGB1"
mygene <- "AKT1"
mygene <- "AKT2"
mygene <- "AKT3"
mygene <- "MTOR"
mygene <- "EIF4A1"
mygene <- "EIF4A2"
mygene <- "EIF4A3"
mygene <- "EIF2S1"
mygene <- "EIF2S2"
mygene <- "EIF2S3"
mygene <- "PIK3R5"
mygene <- "SF3B1"
mygene <- "SRSF2"
mygene <- "JAK2"
mygene <- "GAPDH"
mygene <- "CD200"
mygene <- "BCOR"
# The following has changed with Seurat version 5. Update if you get an error.
metadata$mygene <- LayerData(aml, layer = "data")[mygene,]

# Filter for AML cells at Dx or BM cells. Also remove "healthy" cells from AML patients since their gene expression may be aberrant
metadata.filter <- metadata %>% filter(grepl("AML.*D0", orig.ident) & grepl("-like", CellType) | grepl("BM", orig.ident)) %>%
  mutate(Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))

# The bar plot shows the mean expression across cell types in healthy donors (green) and malignant cells in AML patients at diagnossis (red)
p1 <- metadata.filter %>% group_by(CellType) %>%
  summarize(n = n(), mean_mygene = mean(mygene), Donor = unique(Donor)) %>%
  ggplot(aes(x = CellType, y = mean_mygene, fill = Donor)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(Healthy = "#66cdaa", AML = "#cd5c5c")) +
  ylab("Mean normalized expression") +
  ggtitle(mygene) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

# Visualize the plot 
p1
                                       
# The sina/violin plot shows expression in every individual cell (symbol)
p2 <- metadata.filter %>%
  ggplot(aes(x = CellType, y = mygene, color = Donor)) +
  geom_sina(scale = "width") +
  geom_violin(scale = "width", draw_quantiles = 0.5, color = "black", fill = NA) +
  ylab("Normalized expression, log(TP10K+1)") +
  ggtitle(mygene) +
  scale_x_discrete(drop = F) +
  scale_color_manual(values = c(Healthy = "#66cdaa", AML = "#cd5c5c")) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

# Visualize the plot 
p2
                                       
# Save pdf
pdf(paste0(mygene, "_plots.pdf"), width = 9, height= 9)
plot_grid(p1, p2, ncol = 1)
dev.off()

# Note: you cannot back-calculate the number of transcripts per 10K based on the mean because the mean of log-transformed values differs from the log of the mean. Example for GAPDH:
cells.seu <- subset(aml, CellType == "HSC")
cells.id <- colnames(cells.seu)[grepl("BM", colnames(cells.seu))]
counts.mat <- LayerData(subset(aml, cells = cells.id), layer = "counts")
norm_counts <- sweep(counts.mat, 2, colSums(counts.mat), FUN = "/") * 10000
mean(norm_counts["GAPDH",])
mean( log1p(norm_counts["GAPDH",]) )
# For HSC, there are 3.43 GAPDH transcripts per 10K transcripts, and the mean normalized value is 1.0
# For Mono, there are 10 GAPDH transcripts per 10K transcripts, and the mean normalized value is 2.0
# For T, there are 1.71 GAPDH transcripts per 10K transcripts, and the mean normalized value is 0.46

# Plot signature ----------------------------------------------------------------------------------

# Define signature name and components
mysig <- "EIF4F_complex"
mygenes <- c("EIF4E", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4G1", "EIF4G2", "EIF4G3", "EIF4B", "EIF4H")
aml_sign <- AddModuleScore(aml, features = list(mygenes), name = mysig)

# Extract metadata
metadata <- as_tibble(aml_sign@meta.data, rownames = "cell") %>%
  rename_with(~ mysig, paste0(mysig, "1"))

# Like above, filter for AML cells at Dx or BM cells. Also remove "healthy" cells from AML patients since their gene expression may be aberrant
metadata.filter <- metadata %>% filter(grepl("AML.*D0", orig.ident) & grepl("-like", CellType) | grepl("BM", orig.ident)) %>%
  mutate(Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))
            
# The sina/violin plot shows expression in every individual cell (symbol)
p3 <- metadata.filter %>%
  filter(CellType %in% c("HSC", "Prog", "GMP", "HSC-like", "Prog-like", "GMP-like")) %>%
  mutate(CellType = factor(CellType, levels = c("HSC", "HSC-like", "Prog", "Prog-like", "GMP", "GMP-like"))) %>%
  ggplot(aes(x = CellType, y = !!sym(mysig), color = Donor)) +
  geom_sina(scale = "width") +
  geom_violin(scale = "width", draw_quantiles = 0.5, color = "black", fill = NA) +
  #stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black") +
  ylab(paste(mysig, "Signature Score")) +
  scale_x_discrete(drop = F) +
  scale_color_manual(values = c(Healthy = "#66cdaa", AML = "#cd5c5c")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))
            
# Visualize the plot 
p3
