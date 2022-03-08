setwd("~/Documents/work_adj/05_23_20_crypt_single_cell/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
source("camMonocleHelper.R")

# load(file = "cds_full.rds",verbose = T)

# cds.full <- preprocess_cds(cds = cds.full,num_dim = 100,verbose = T)
# cds.full <- reduce_dimension(cds = cds.full,verbose = T)
# cds.full <- cluster_cells(cds = cds.full,verbose = T)
# cds.full <- learn_graph(cds = cds.full,verbose = T)

# marker_test_res <- top_markers(cds = cds.full,reference_cells = 5000,verbose = T)
# spl <- split(x = marker_test_res,f = list(marker_test_res$cell_group))
# top_specific_markers <- do.call("rbind",lapply(X = spl,FUN = function(df) {
#   df10 <- df[order(df$marker_score,decreasing = T),]
#   df10[1:3,]
# }))
# top_specific_markers_ids <- unique(top_specific_markers$gene_id)
# g <- plot_genes_by_group(cds = cds.full,markers = top_specific_markers_ids,group_cells_by = "cluster",ordering_type = "cluster_row_col",max.size = 3)

# assignments <- read.csv(file = "cluster_assignments_full.csv",header = T)
# cds.full <- assign.clusters(cds = cds.full,assignments = assignments)
# cds.full <- assign.clusters(cds = cds.full,assignments = assignments,column = 3,new.col.name = "annotation.2")

# sample.names <- c("D00","D06","D09","D12","D16","D19","D23","IL10")
# injury.status <- c("None","DSS","DSS","DSS","DSS","DSS","DSS","IL10")
# colData(cds.full)$sample.name <- sapply(X = colData(cds.full)$sample,FUN = function(x) {
#   return(sample.names[x])
# })
# colData(cds.full)$injury <- sapply(X = colData(cds.full)$sample,FUN = function(x) {
#   return(injury.status[x])
# })

# barcodes <- data.frame(colData(cds.full))
# breakdown <- table(barcodes$sample.name,barcodes$annotation)
# breakdown10 <- t(apply(X = breakdown,MARGIN = 1,FUN = function(x) {
#   return(x/sum(x)*100)
# }))

# d <- data.frame(breakdown10)
# d$sample <- row.names(d)
# d10 <- melt(data = d,id.var = "sample")
# colnames(d10) <- c("sample","celltype","abundance")
# d10$celltype <- factor(d10$celltype,levels = c("Epithelium","Fibroblast","Myofibroblast","Endothelium","Lymphatic","Glia","T.cell","B.cell","Macrophage","Neutrophil","Dendritic.cell"))
# pal <- c(rgb(1.0,0,0),
#          rgb(0,1.0,0),rgb(0,0.5,0),
#          rgb(1.0,1.0,0),rgb(0.67,0.67,0),rgb(0.33,0.33,0),
#          rgb(0,0,1.0),rgb(0,0,0.5),
#          rgb(0,1.0,1.0),rgb(0,0.67,0.8),rgb(0,0.33,0.6))
# h <- ggplot(data = d10,mapping = aes(x = sample,y = abundance)) +
#   geom_bar(stat = "identity",mapping = aes(fill = celltype),position = "stack") +
#   scale_fill_manual(values = pal) +
#   theme_classic(base_size = 16)
# print(h)

# l <- plot_cells(cds = cds.full,color_cells_by = "annotation",show_trajectory_graph = F,group_label_size = 4,cell_size = 1) +
#   scale_color_calc()
# print(l)