setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
library(gplots)
source("camMonocleHelper.R")

# load(file = "cds_full.rds")
# d <- read.csv(file = "table_stem_repro_dynamics.csv",header = T)


# cds.full.subset <- cds.full[map.gene.short.names(cds = cds.full,gene.short.names = d$gene_short_name),]
# colData(cds.full.subset)$condition <- with(colData(cds.full.subset),paste(annotation,sample.name,sep = "-"))
# agg <- aggregate.expression.by.factor(cds = cds.full.subset,grouping = "condition",do.print = T)

## computes mean of each cell type to screen changing cell types
# means <- colMeans2(x = agg)
# sds <- colSds(x = agg) / sqrt(dim(agg)[[1]])
# names(means) <- colnames(agg)
# names(sds) <- colnames(agg)

## prepares the plot
# comb <- data.frame(condition = names(means),mean = means,sd = sds)
# comb$cell.type <- sapply(X = comb$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = "-")[[1]][1]
# })
# comb$sample <- sapply(X = comb$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = "-")[[1]][2]
# })
# to.focus <- c("Epithelium","Endothelium","Fibroblast","Myofibroblast","T cell","Lymphatic")
# comb$focus <- ifelse(test = comb$cell.type %in% to.focus,1,0.2)

## plots
# g <- ggplot(data = comb,mapping = aes(x = sample,y = mean)) +
#   geom_line(mapping = aes(group = cell.type,colour=cell.type,alpha=focus),size=3) +
#   geom_ribbon(mapping = aes(ymin = mean-sd,ymax = mean+sd,fill = cell.type,group=cell.type),alpha=0.1) +
#   theme_classic(base_size = 24)
# print(g)