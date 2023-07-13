setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/notacliu_shared/data/2022/07_01_22_Regev_UC_reanalysis/")

library(monocle3)
library(dplyr)
library(readxl)
library(ggthemes)
source("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/camMonocleHelper.R")

# load(file = "cds_full.rds",verbose = T)

## finds markers
# mtr <- top_markers(cds = cds.full,verbose = T)
# write.table(x = mtr,file = "table_cluster_markers.csv",sep = ",",row.names = F)
# tm <- filter.marker.results(marker_test_res = mtr,num.to.keep = 5,criterion = "specificity")
# h <- plot_genes_by_group(cds = cds.full,markers = tm,max.size = 5)
# ggsave(filename = "plot_cluster_markers.svg",plot = h)

## make cluster assignments
# assignments <- read_xlsx(path = "cluster_assignments.xlsx",sheet = 1,col_names = T) %>% data.frame
# cds.full <- assign.clusters(cds = cds.full,assignments = assignments,column = 2,new.col.name = "Cell.Type")
# colData(cds.full)$cell.group <- clusters(cds.full)
# save(cds.full,file = "cds_full.rds")

## filters out doublet clusters
# barcodes <- colData(cds.full) %>% data.frame
# cds.fil <- cds.full[,subset(barcodes,!barcodes$Cell.Type %in% c("Doublets - Myeloid","Doublets - Mesenchyme")) %>% row.names]
# save(cds.fil,file = "cds_filtered.rds")
# g <- plot_cells(cds = cds.fil,color_cells_by = "Cell.Type",show_trajectory_graph = F,group_label_size = 6,label_cell_groups = F)
# print(g)
# ggsave(filename = "plot_UMAP.png",plot = g,dpi = 300)

## checks abundance of different biopsy types in each cluster
# tab <- table(colData(cds.fil)$Health,colData(cds.fil)$Cell.Type)
# tab.norm.tot <- sweep(x = tab,MARGIN = 1,STATS = rowSums(tab),FUN = "/")
# tabsure <- tab[,-7]
# tabsure.norm.tot <- sweep(x = tabsure,MARGIN = 1,STATS = rowSums(tabsure),FUN = "/")
# tabdf <- data.frame(tab.norm.tot)
# colnames(tabdf) <- c("Health","Cell.Type","Freq")
# tabdf$Health <- factor(tabdf$Health,levels=c("Healthy","Non-inflamed","Inflamed"))
# k <- ggplot(data = tabdf,mapping = aes(x = Cell.Type,y = Freq*100)) +
#   geom_col(mapping = aes(fill = Health),position = position_dodge()) +
#   scale_fill_colorblind() +
#   ylab("Abundance (% of recovered cells)") +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(k)
# ggsave(filename = "plot_cell_type_abundance.svg",plot = k)

## Are samples labeled Non-inflamed really that different from Inflamed? Are annotations informative?
## Studies the concordance of genes by analyzing top markers by health status
# mtra <- top_markers(cds = cds.fil,group_cells_by = "Health",genes_to_test_per_group = 500,verbose = T)
# tm <- filter.marker.results(marker_test_res = mtra,num.to.keep = 20,criterion = "specificity")
# print(plot_genes_by_group(cds = cds.fil,group_cells_by = "Health",markers = tm))
# write.table(x = mtra,file = "table_markers_by_health.csv",sep = ",",row.names = F)

## analyzes top markers 

## plots expression of select genes
# goi <- c("CAPG","APRT","SULF2","WFDC2")
# cds.sub <- cds.fil[goi,]
# print(plot.expression.with.dots(cds.subset = cds.sub,grouping = "Cell.Type"))
