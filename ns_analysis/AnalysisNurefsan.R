# 230411,Nurefsan Sariipek
# Check expression of genes of interest in data from single-cell AML paper (van Galen et al, Cell 2019)

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggforce)
library(ggplot2)
library(cowplot)
library(randomcoloR)
library(readxl)
library(data.table)
library(janitor)

# Clean the enviroment
rm(list=ls())

# Frequently used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load AML paper data. See https://github.com/petervangalen/reanalyze-aml2019 first for how to obtain this file
# For Nurefsan: setwd("/Users/dz855/Documents/R/AML2019/")
# For Peter: setwd("~/DropboxMGB/GitHub/reanalyze-aml2019/")
aml <- readRDS(file = "Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))

# Add a variable to metadata to identify the FLT3 mutation across the samples
normal <- c(paste0("BM", 1:4), "BM5.34p", "BM5.34p38n")
flt3_mutated <- c("AML210A", "AML419A", "AML997", "AML329", "AML328")
flt3_wt <- c("AML1012", "AML314", "AML371", "AML420B", "AML475", "AML556", "AML707B", "AML722B", "AML870", "AML916", "AML921A")
flt3_wt_cell_line <- c("MUTZ3", "OCI.AML3")
aml$FLT3 <- case_when(grepl(paste(flt3_mutated, collapse = "|"), aml$orig.ident) ~ "FLT3_mutated_AML",
                      grepl(paste(flt3_wt, collapse = "|"), aml$orig.ident) ~ "FLT3_wildtype_AML",
                      grepl(paste(flt3_wt_cell_line, collapse = "|"), aml$orig.ident) ~ "FLT3_wildtype_CellLine",
                      grepl(paste(normal, collapse = "|"), aml$orig.ident) ~ "Normal")

# Add gene expression as metadata
metadata <- as_tibble(aml@meta.data, rownames = "cell")
# Select your gene of interest (One of these)
mygene <- "HOXA9"
mygene <- "HAVCR2"
mygene <- "LGALS9"
mygene <- "CEACAM1"
mygene <- "HMGB1"
# Add expression data to metadata
metadata$mygene <- LayerData(aml, layer = "data")[mygene,]

# Remove "normal" cells from AML patients since their gene expression may be aberrant 
metadata.filter <- metadata %>% filter(grepl("AML|MUTZ", orig.ident) & grepl("-like", CellType) | grepl("BM", orig.ident))
# Compare cell number before and after filtering (check that it makes sense)
metadata %>% group_by(orig.ident) %>% count() %>%
  left_join(count(group_by(metadata.filter, orig.ident)), by = "orig.ident") %>% print(n = 100)
# Filter out samples with <10 cells
metadata.filter <- metadata.filter %>% group_by(orig.ident) %>% filter(n() >= 10)

# Wrangle data into a better order for visualization later
metadata.filter <- metadata.filter %>%
  mutate(FLT3 = factor(FLT3, levels = c("Normal", "FLT3_wildtype_AML", "FLT3_mutated_AML", "FLT3_wildtype_CellLine"))) %>%
  arrange(FLT3)
orig.order <- unique(metadata.filter$orig.ident)
metadata.filter$orig.ident <- factor(metadata.filter$orig.ident, levels = orig.order)

# Plot the mean gene expression across the samples, split by samples. Filter samples with <10 cells.
p1 <- metadata.filter %>% group_by(orig.ident) %>%
  summarize(n = n(), mean_mygene = mean(mygene), FLT3 = unique(FLT3)) %>%
  ggplot(aes(x = orig.ident, y = mean_mygene, fill = FLT3)) +
  ylab("Mean normalized expression log(TP10K+1)") +
  geom_bar(stat="identity") +
  scale_x_discrete(drop=F) +
  ggtitle(mygene) +
  theme_bw() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

# Visualize the plot
p1

# Sina plot
p2 <- metadata.filter %>%
  ggplot(aes(x = orig.ident, y = mygene, color = FLT3)) +
  geom_violin(scale = "width") +
  geom_sina(scale = "width") +
  ylab("Normalized expression log(TP10K+1)") +
  ggtitle(mygene) +
  scale_x_discrete(drop = F) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 8, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

# Visualize the plot        
p2
