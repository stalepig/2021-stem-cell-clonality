setwd("/Volumes/Macintosh HD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)
library(readxl)
library(magrittr)
library(tibble)
library(ggplot2)
source(file = "camMonocleHelper.R")

# load(file = "cds_all_aligned.rds",verbose = T)

# goi <- c("Ptprc","Krt8","Epcam","Krt14","Cdh1","Lyve1","Col6a1","Fn1","Krt7","Chga","Cd3g","Jchain","Plp1","Lyz2","Cd19","Gata3","Cd74","Dclk1","Vil1","Igkc","Vwf","Pecam1","S100a8","S100a9","Acta2","Hba-a2","Hba-a1")
# g <- plot_genes_by_group(cds = cds.allnn,markers = goi,max.size = 6)
# ggsave(filename = "plot_dots_global_markers.png",plot = g)

# assignments <- read_xlsx(path = "assignments_all.xlsx",col_names = T) %>% data.frame
# cds.allnn <- assign.clusters(cds = cds.allnn,assignments = assignments)

# tm <- top_markers(cds = cds.allnn,group_cells_by = "annotation",verbose = T)
# g10 <- plot_genes_by_group(cds = cds.allnn,markers = filter.marker.results(tm,5,"pseudo_R2"),group_cells_by = "annotation")
# ggsave(filename = "plot_dots_global_markers_unbiased.png",plot = g10)
# write.table(x = tm,file = "table_global_markers_unbiased.csv",sep = ",",row.names = F)
# save(cds.allnn,file = "cds_all_aligned.rds")

## now subclassifies each of the original clusters
# barcodes <- colData(cds.allnn) %>% data.frame
# cds.sub <- cds.allnn[,subset(barcodes,barcodes$annotation == "Immune") %>% row.names]
# mtr <- top_markers(cds = cds.sub,verbose = T)
# print(plot_genes_by_group(cds = cds.sub,markers = filter.marker.results(mtr,5,"marker_score"),group_cells_by = "cluster"))
# cds.sub <- cds.allnn[,subset(barcodes,barcodes$annotation == "Mesenchyme") %>% row.names]
# mtr <- top_markers(cds = cds.sub,verbose = T)
# print(plot_genes_by_group(cds = cds.sub,markers = filter.marker.results(mtr,5,"pseudo_R2"),group_cells_by = "cluster"))
# cds.sub <- cds.allnn[,subset(barcodes,barcodes$annotation == "Endothelium") %>% row.names]
# mtr <- top_markers(cds = cds.sub,verbose = T)
# print(plot_genes_by_group(cds = cds.sub,markers = filter.marker.results(mtr,5,"marker_score"),group_cells_by = "cluster"))
# cds.sub <- cds.allnn[,subset(barcodes,barcodes$annotation == "Epithelium") %>% row.names]
# mtr <- top_markers(cds = cds.sub,verbose = T)
# print(plot_genes_by_group(cds = cds.sub,markers = filter.marker.results(mtr,5,"specificity"),group_cells_by = "cluster"))

## adds in additional annotations
# assignments <- read_xlsx(path = "assignments_all.xlsx",col_names = T) %>% data.frame
# cds.allnn <- assign.clusters(cds = cds.allnn,assignments = assignments,column = 3,new.col.name = "annotation.2")
# cds.allnn <- assign.clusters(cds = cds.allnn,assignments = assignments,column = 4,new.col.name = "annotation.3")
# cds.allnn <- assign.clusters(cds = cds.allnn,assignments = assignments,column = 5,new.col.name = "filter")

## filters out keratinocytes and doublet/dead cells
# barcodes <- colData(cds.allnn) %>% data.frame
# cds.allnn.fil <- cds.allnn[,subset(barcodes,is.na(barcodes$filter)) %>% row.names]
# barcodes <- colData(cds.allnn.fil) %>% data.frame
# cds.allnn.fil2 <- cds.allnn.fil[,subset(barcodes,!barcodes$annotation == "Keratinocyte") %>% row.names]
# save(cds.allnn.fil2,file = "cds_all_aligned_fil.rds")