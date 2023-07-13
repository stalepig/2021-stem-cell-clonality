setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/notacliu_shared/data/2022/07_01_22_Regev_UC_reanalysis/")

library(Matrix)
library(monocle3)
library(dplyr)
source("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/camMonocleHelper.R")

load(file = "cds_full.rds",verbose = T)

# cds.full <- align_cds(cds = cds.full,alignment_group = "Subject",verbose = T)
# cds.full <- reduce_dimension(cds = cds.full,verbose = T)
# cds.full <- cluster_cells(cds = cds.full,verbose = T)
# cds.full <- learn_graph(cds = cds.full,verbose = T)

## Repeats the clustering to get fewer number of clusters
# cds.full <- cluster_cells(cds = cds.full,resolution = 1e-5,verbose = T)

## finds markers
# mtr <- top_markers(cds = cds.full,verbose = T)
# tm <- filter.marker.results(marker_test_res = mtr,num.to.keep = 5,criterion = "specificity")
# print(plot_genes_by_group(cds = cds.full,markers = tm,max.size = 5))

## filters out clusters 2, 4, 6 based on reduced expression of KLF4/5 and EPCAM transcripts
## removes clusters 14 and 11 as these are likely fusions
# colData(cds.full)$cell.group <- clusters(cds.full)
# barcodes <- colData(cds.full) %>% data.frame
# cds.fil <- cds.full[,subset(barcodes,!barcodes$cell.group %in% c(2,4,6,11,14)) %>% row.names]

## resets the analysis
# cds.fil <- preprocess_cds(cds = cds.fil,num_dim = 100,verbose = T)
# cds.fil <- align_cds(cds = cds.fil,alignment_group = "Subject",verbose = T)
# cds.fil <- reduce_dimension(cds = cds.fil,verbose = T)
# cds.fil <- cluster_cells(cds = cds.fil,resolution = 1e-5,verbose = T)

## based on some diagnostics, the filtering was overkill, and should not have removed clusters 2, 4, and 6 as these seem to be somewhat intermediary phenotypes?
## the filtered UMAP looks like crap