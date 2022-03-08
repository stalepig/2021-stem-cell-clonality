setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

source("camMonocleHelper.R")
library(monocle3)
library(magrittr)
library(dplyr)
library(ggthemes)

# load(file = "cds_epi_filtered2.rds",verbose = T)

## adds Ly6a gate
# cds.epi.fil.fil <- add.gene.sum.column(cds = cds.epi.fil.fil,genes.of.interest = c("Ly6a"),col.name = "idx.Ly6a",use.raw = T)
# colData(cds.epi.fil.fil)$Ly6a.status <- ifelse(colData(cds.epi.fil.fil)$idx.Ly6a > 2,"hi","lo")

## splits into hi and lo populations and does the trajectory analysis separately
## works on the lo population
# barcodes <- colData(cds.epi.fil.fil) %>% data.frame
# cds.epi.lo <- cds.epi.fil.fil[,subset(barcodes,barcodes$Ly6a.status=="lo") %>% row.names]
# cds.epi.lo <- recluster.cds(cds = cds.epi.lo,cluster.method = "louvain")

## assigns clusters
# mtr <- top_markers(cds = cds.epi.lo,genes_to_test_per_group = 40,verbose = T)
# tm <- filter.marker.results(marker_test_res = mtr,num.to.keep = 40,criterion = "marker_score")
# g <- plot_genes_by_group(cds = cds.epi.lo,markers = tm,max.size = 4)
# print(g)
# anno <- read.csv(file = "cluster_assignments_epithelium_Ly6a_lo.csv",header = T)
# cds.epi.lo <- assign.clusters(cds = cds.epi.lo,assignments = anno,column = "Assignment.2",new.col.name = "Assignment")
# cds.epi.lo <- assign.clusters(cds = cds.epi.lo,assignments = anno,column = "Assignment",new.col.name = "Assignment.2")
# colData(cds.epi.lo)$cluster.name <- paste(colData(cds.epi.lo)$Assignment,colData(cds.epi.lo)$Assignment.2,sep = "")
# mtr10 <- top_markers(cds = cds.epi.lo,genes_to_test_per_group = 40,verbose = T,group_cells_by = "cluster.name")
# tm10 <- filter.marker.results(marker_test_res = mtr,num.to.keep = 4,criterion = "pseudo_R2")
# g10 <- plot_genes_by_group(cds = cds.epi.lo,markers = tm10,max.size = 4,group_cells_by = "cluster.name")
# print(g10)

## computes Klf4/Klf5 ratio as a proxy for crypt position -- maybe
# m.klf <- aggregate.expression.by.factor(cds = cds.epi.lo[map.gene.short.names(cds = cds.epi.lo,gene.short.names = c("Klf4","Klf5")),],grouping = "cluster.name")
# m.klf10 <- t(apply(X = m.klf,MARGIN = 1,FUN = function(x) {
#   return((x - min(x))/(max(x) - min(x)))
# }))
# m.klf20 <- apply(X = m.klf10,MARGIN = 2,FUN = function(x) { x[1] - x[2] })
# m.klf30 <- (m.klf20 - min(m.klf20)) / (max(m.klf20) - min(m.klf20))
# colData(cds.epi.lo)$crypt.position <- sapply(X = colData(cds.epi.lo)$cluster.name,FUN = function(x) {
#   return(m.klf30[as.character(x)])
# })

## makes visualizations
# h <- plot_cells(cds = cds.epi.lo,color_cells_by = "crypt.position",show_trajectory_graph = F) +
#   geom_point(mapping=aes(colour=crypt.position))
# print(h)
# colData(cds.epi.lo)$cluster.name <- factor(colData(cds.epi.lo)$cluster.name)
# h10 <- plot_cells(cds = cds.epi.lo,color_cells_by = "cluster.name",show_trajectory_graph = F,group_label_size = 4,cell_size = 1)
# print(h10)

## splits into hi and lo populations and does the trajectory analysis separately
## works on the hi population for DSS
# barcodes <- colData(cds.epi.fil.fil) %>% data.frame
# cds.epi.hi.dss <- cds.epi.fil.fil[,subset(barcodes,barcodes$Ly6a.status=="hi" & barcodes$injury=="DSS") %>% row.names]
# cds.epi.hi.dss <- recluster.cds(cds = cds.epi.hi.dss,cluster.method = "louvain")
# mtr20 <- top_markers(cds = cds.epi.hi.dss,genes_to_test_per_group = 40,verbose = T)
# tm20 <- filter.marker.results(marker_test_res = mtr20,num.to.keep = 4,criterion = "specificity")
# g20 <- plot_genes_by_group(cds = cds.epi.hi.dss,markers = tm20,max.size = 4)
# print(g20)
# anno <- read.csv(file = "cluster_assignments_epithelium_Ly6a_hi_DSS.csv",header = T)
# cds.epi.hi.dss <- assign.clusters(cds = cds.epi.hi.dss,assignments = anno,column = "Assignment.2",new.col.name = "Assignment")
# cds.epi.hi.dss <- assign.clusters(cds = cds.epi.hi.dss,assignments = anno,column = "Assignment",new.col.name = "Assignment.2")
# colData(cds.epi.hi.dss)$cluster.name <- paste(colData(cds.epi.hi.dss)$Assignment,colData(cds.epi.hi.dss)$Assignment.2,sep = "")
# m10.klf <- aggregate.expression.by.factor(cds = cds.epi.hi.dss[map.gene.short.names(cds = cds.epi.hi.dss,gene.short.names = c("Klf4","Klf5")),],grouping = "cluster.name")
# m10.klf10 <- t(apply(X = m10.klf,MARGIN = 1,FUN = function(x) {
#   return((x - min(x))/(max(x) - min(x)))
# }))
# m10.klf20 <- apply(X = m10.klf10,MARGIN = 2,FUN = function(x) { x[1] - x[2] })
# m10.klf30 <- (m10.klf20 - min(m10.klf20)) / (max(m10.klf20) - min(m10.klf20))
# colData(cds.epi.hi.dss)$crypt.position <- sapply(X = colData(cds.epi.hi.dss)$cluster.name,FUN = function(x) {
#   return(m10.klf30[as.character(x)])
# })
# h20 <- plot_cells(cds = cds.epi.hi.dss,color_cells_by = "crypt.position",show_trajectory_graph = F) +
#   geom_point(mapping=aes(colour=crypt.position))
# print(h20)
# colData(cds.epi.hi.dss)$cluster.name <- factor(colData(cds.epi.hi.dss)$cluster.name)
# h30 <- plot_cells(cds = cds.epi.hi.dss,color_cells_by = "cluster.name",show_trajectory_graph = F,group_label_size = 4,cell_size = 1)
# print(h30)

## splits into hi and lo populations and does the trajectory analysis separately
## works on the hi population for Il10-/-
# barcodes <- colData(cds.epi.fil.fil) %>% data.frame
# cds.epi.hi.il10 <- cds.epi.fil.fil[,subset(barcodes,barcodes$Ly6a.status=="hi" & barcodes$injury=="IL10") %>% row.names]
# cds.epi.hi.il10 <- recluster.cds(cds = cds.epi.hi.il10,cluster.method = "louvain")
# mtr30 <- top_markers(cds = cds.epi.hi.il10,genes_to_test_per_group = 40,verbose = T)
# tm30 <- filter.marker.results(marker_test_res = mtr30,num.to.keep = 4,criterion = "marker_score")
# g30 <- plot_genes_by_group(cds = cds.epi.hi.il10,markers = tm30,max.size = 4)
# print(g30)
# anno <- read.csv(file = "cluster_assignments_epithelium_Ly6a_hi_Il10.csv",header = T)
# cds.epi.hi.il10 <- assign.clusters(cds = cds.epi.hi.il10,assignments = anno,column = "Assignment.2",new.col.name = "Assignment")
# cds.epi.hi.il10 <- assign.clusters(cds = cds.epi.hi.il10,assignments = anno,column = "Assignment",new.col.name = "Assignment.2")
# colData(cds.epi.hi.il10)$cluster.name <- paste(colData(cds.epi.hi.il10)$Assignment,colData(cds.epi.hi.il10)$Assignment.2,sep = "")
# m20.klf <- aggregate.expression.by.factor(cds = cds.epi.hi.il10[map.gene.short.names(cds = cds.epi.hi.il10,gene.short.names = c("Klf4","Klf5")),],grouping = "cluster.name")
# m20.klf10 <- t(apply(X = m20.klf,MARGIN = 1,FUN = function(x) {
#   return((x - min(x))/(max(x) - min(x)))
# }))
# m20.klf20 <- apply(X = m20.klf10,MARGIN = 2,FUN = function(x) { x[1] - x[2] })
# m20.klf30 <- (m20.klf20 - min(m20.klf20)) / (max(m20.klf20) - min(m20.klf20))
# colData(cds.epi.hi.il10)$crypt.position <- sapply(X = colData(cds.epi.hi.il10)$cluster.name,FUN = function(x) {
#   return(m20.klf30[as.character(x)])
# })
# h40 <- plot_cells(cds = cds.epi.hi.il10,color_cells_by = "crypt.position",show_trajectory_graph = F) +
#   geom_point(mapping=aes(colour=crypt.position))
# print(h40)
# colData(cds.epi.hi.il10)$cluster.name <- factor(colData(cds.epi.hi.il10)$cluster.name)
# h50 <- plot_cells(cds = cds.epi.hi.il10,color_cells_by = "cluster.name",show_trajectory_graph = F,group_label_size = 4,cell_size = 1)
# print(h50)

## performs graph analysis another way (Moran's I)
## works on Ly6a-lo
# gtr <- graph_test(cds = cds.epi.lo,neighbor_graph = "principal_graph",verbose = T)
# gtr10 <- gtr %>% na.omit %>% filter(morans_I > 0.2)
# gtr20 <- gtr10[order(gtr10$morans_I,decreasing = T),]
# mapping.genes <- c("Dclk1","Aqp8","Muc2","Cd24a","Lgr5","Ascl2","Clca1","Neurod1","Spink4","Ptma","Guca2a","Krt20","Ly6g","Best2","Gsdmc4","Birc5")
# ha <- plot_cells(cds = cds.epi.lo,genes = mapping.genes,show_trajectory_graph = F,norm_method = "log",label_cell_groups = F)

## works on Ly6a-hi DSS since Ly6a-hi cells in Il10-/- are a little bit different
# gtra <- graph_test(cds = cds.epi.hi.dss,neighbor_graph = "principal_graph",verbose = T)
# gtra10 <- gtra %>% na.omit %>% filter(morans_I > 0.2)
# gtra20 <- gtra10[order(gtra10$morans_I,decreasing = T),]
# mapping.genes.dss <- c("Spdef","Saa1","Gpx2","Gsdmc4","Clca4a","Tgm3","Slc9a3","Mgat4c","Krt7","Cldn4")
# hb <- plot_cells(cds = cds.epi.hi.dss,genes = gtra20[201:220,"gene_short_name"],show_trajectory_graph = F,norm_method = "log")
# print(hb)

## makes reciprocal plots
# mapping.genes.merged <- c(mapping.genes,mapping.genes.dss) %>% unique
# hc <- plot_cells(cds = cds.epi.lo,genes = mapping.genes.merged,show_trajectory_graph = F,norm_method = "log",label_cell_groups = F)
# print(hc)
# hd <- plot_cells(cds = cds.epi.hi.dss,genes = mapping.genes.merged,show_trajectory_graph = F,norm_method = "log",label_cell_groups = F)
# print(hd)

## plots base trajectories
# he <- plot_cells(cds = cds.epi.lo,show_trajectory_graph = T,cell_size = 3,alpha=0.7,trajectory_graph_color = "black",trajectory_graph_segment_size = 2,label_leaves = F,label_branch_points = F,label_cell_groups = F)
# print(he)
# hf <- plot_cells(cds = cds.epi.hi.dss,show_trajectory_graph = T,cell_size = 3,alpha=0.7,trajectory_graph_color = "black",trajectory_graph_segment_size = 2,label_leaves = F,label_branch_points = F,label_cell_groups = F)
# print(hf)