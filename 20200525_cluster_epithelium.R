setwd("~/Documents/work_adj/05_23_20_crypt_single_cell/monocle3/")

library(monocle3)
library(magrittr)
library(dplyr)
library(ggthemes)
source("camMonocleHelper.R")

# load(file = "cds_full.rds")
# barcodes <- colData(cds.full) %>% data.frame
# cds.epi <- cds.full[,subset(barcodes,barcodes$annotation=="Epithelium") %>% row.names]
# cds.epi <- recluster.cds(cds.epi)

# marker_test_res <- top_markers(cds = cds.epi,reference_cells = 5000,verbose = T)
# top_specific_markers_ids <- filter.marker.results(marker_test_res,num.to.keep = 6,criterion = "marker_score")
# g <- plot_genes_by_group(cds = cds.epi,markers = top_specific_markers_ids,group_cells_by = "cluster",ordering_type = "cluster_row_col",max.size = 4)
# print(g)

# colData(cds.epi)$injury <- factor(colData(cds.epi)$injury,levels=c("None","DSS","IL10"))
# h <- plot_cells(cds = cds.epi,color_cells_by = "injury",show_trajectory_graph = F,label_cell_groups = F) +
  # scale_color_colorblind()
# print(h)

## Based on analysis of day 00, one of the clusters is full of dying cells and should be excluded
# barcodes.epi <- colData(cds.epi) %>% data.frame
# cds.epi.d00 <- cds.epi[,subset(barcodes.epi,barcodes.epi$sample.name=="D00") %>% row.names]
# cds.epi.d00 <- cluster_cells(cds = cds.epi.d00,verbose = T,cluster_method = "leiden")
# cds.epi.d00 <- learn_graph(cds = cds.epi.d00,verbose = T)
# marker_test_res <- top_markers(cds = cds.epi.d00,verbose = T)
# top_specific_markers_ids <- filter.marker.results(marker_test_res,num.to.keep = 20,criterion = "marker_score")
# g <- plot_genes_by_group(cds = cds.epi.d00,markers = top_specific_markers_ids,group_cells_by = "cluster",ordering_type = "cluster_row_col",max.size = 4)
# print(g)

## remove cluster 1 and re-cluster and filter
# colData(cds.epi)$cell.group <- clusters(cds.epi)
# barcodes <- colData(cds.epi) %>% data.frame
# cds.epi.fil <- cds.epi[,subset(barcodes,!barcodes$cell.group %in% c(1)) %>% row.names]
# cds.epi.fil <- recluster.cds(cds = cds.epi.fil)
# marker_test_res <- top_markers(cds = cds.epi.fil,verbose = T)
# top_specific_markers_ids <- filter.marker.results(marker_test_res,num.to.keep = 8,criterion = "marker_score")
# g <- plot_genes_by_group(cds = cds.epi.fil,markers = top_specific_markers_ids,group_cells_by = "cluster",ordering_type = "cluster_row_col",max.size = 4)
# print(g)

## remove cluster 4 (myeloid-epithelial doublets) and re-cluster and filter
# colData(cds.epi.fil)$cell.group <- clusters(cds.epi.fil)
# barcodes.fil <- colData(cds.epi.fil) %>% data.frame
# cds.epi.fil.fil <- cds.epi.fil[,subset(barcodes.fil,!barcodes.fil$cell.group %in% c(4)) %>% row.names]
# cds.epi.fil.fil <- recluster.cds(cds = cds.epi.fil.fil)
# cds.epi.fil.fil <- cluster_cells(cds = cds.epi.fil.fil,cluster_method = "leiden")
# marker_test_res <- top_markers(cds = cds.epi.fil.fil,verbose = T)
# top_specific_markers_ids <- filter.marker.results(marker_test_res,num.to.keep = 8,criterion = "specificity")
# g <- plot_genes_by_group(cds = cds.epi.fil.fil,markers = top_specific_markers_ids,group_cells_by = "cluster",ordering_type = "cluster_row_col",max.size = 4)
# print(g)
# h <- plot_cells(cds = cds.epi.fil.fil,color_cells_by = "injury",show_trajectory_graph = F,label_cell_groups = F) +
#   geom_point(mapping=aes(colour=injury,size=injury,alpha=injury)) +
#   scale_size_manual(values = c(2,0.5,0.5)) +
#   scale_alpha_manual(values = c(0.8,0.3,0.3)) +
#   scale_color_colorblind()
# print(h)
# colData(cds.epi.fil.fil)$cell.group <- clusters(cds.epi.fil.fil)
# l <- plot_cells(cds = cds.epi.fil.fil,show_trajectory_graph = F) +
#   facet_wrap(facets = .~sample.name) +
#   geom_point(mapping=aes(colour=cell.group),size=1)
# print(l)
