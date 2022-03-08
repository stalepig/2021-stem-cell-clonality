setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
library(dplyr)
library(scales)
library(umap)
library(gplots)
source("camMonocleHelper.R")

# d <- read.csv(file = "table_stem_repro_dynamics.csv",header = T)
# m <- data.matrix(d[,3:9])
# row.names(m) <- d$gene_short_name

# hm <- heatmap.2(x = m,Colv = F,dendrogram = "none",trace = "none",cexRow = 0.5,lhei = c(2,10))
# print(hm)

## makes stem cell trajectory
# load(file = "cds_epi_stem.rds")
# cds.stem <- order_cells(cds = cds.stem,verbose = T)
# colData(cds.stem)$time.bin <- cut(x = pseudotime(cds.stem),breaks = 5)
# cds.stem <- add.gene.sum.column(cds = cds.stem,genes.of.interest = c("Ly6a"),col.name = "idx.Ly6a",use.raw=T)
# colData(cds.stem)$Ly6a.status <- ifelse(log10(colData(cds.stem)$idx.Ly6a+0.1) > 0.783,"hi","lo")
# colData(cds.stem)$Ly6a.status <- factor(colData(cds.stem)$Ly6a.status,levels=c("lo","hi"))
# goi <- c("Clu","Ly6a","Lgr5","Ascl2","Anxa2","Sulf2")
# goi.id <- map.gene.short.names(cds = cds.stem,gene.short.names = goi)
# g <- plot_genes_violin(cds_subset = cds.stem[goi.id,],group_cells_by = "time.bin",log_scale = F,ncol = 2) +
#   theme_classic(base_size = 16)
# print(g)
# ggsave(filename = "plot_epi_stem_pseudotime_genes.svg",plot = g)

## plots Clu on stem cell trajectory
h <- plot_cells(cds = cds.stem,genes = c("Clu"),cell_size = 3,trajectory_graph_segment_size = 3,trajectory_graph_color = "red")
ggsave(filename = "plot_epi_stem_pseudotime_Clu.png",plot = h)