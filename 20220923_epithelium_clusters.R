setwd("/Volumes/Macintosh HD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)
library(readxl)
library(magrittr)
library(tibble)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(gplots)
library(dendextend)
library(pheatmap)
source(file = "camMonocleHelper.R")

# load(file = "cds_all_aligned_fil.rds",verbose = T)

## takes epithelium and filters out Il10 samples
# barcodes <- colData(cds.allnn.fil2) %>% data.frame
# cds.epi <- cds.allnn.fil2[,subset(barcodes,barcodes$annotation=="Epithelium" & !barcodes$injury=="IL10") %>% row.names]
# cds.epi <- reduce_dimension(cds = cds.epi,verbose = T)
# cds.epi <- cluster_cells(cds = cds.epi,verbose = T)
# mtr <- top_markers(cds = cds.epi,verbose = T)
# print(plot_genes_by_group(cds = cds.epi,markers = filter.marker.results(mtr,5,"marker_score")))
# assignments <- read_xlsx(path = "assignments_epithelium.xlsx",col_names = T) %>% data.frame
# cds.epi <- assign.clusters(cds = cds.epi,assignments = assignments,column = 2,new.col.name = "annotation")
# cds.epi <- assign.clusters(cds = cds.epi,assignments = assignments,column = 3,new.col.name = "filter")

## filters some more
# barcodes <- colData(cds.epi) %>% data.frame
# cds.epi.fil <- cds.epi[,subset(barcodes,is.na(barcodes$filter)) %>% row.names]
# cds.epi.fil <- reduce_dimension(cds = cds.epi.fil,verbose = T)
# cds.epi.fil <- cluster_cells(cds = cds.epi.fil,verbose = T)
# mtr <- top_markers(cds = cds.epi.fil,verbose = T)
# print(plot_genes_by_group(cds = cds.epi.fil,markers = filter.marker.results(mtr,5,"marker_score")))

## cluster 10 is a mesenchymal doublet
## cluster 8 is outlier
# colData(cds.epi.fil)$cell.group <- clusters(cds.epi.fil)
# barcodes <- colData(cds.epi.fil) %>% data.frame
# cds.epi.fil2 <- cds.epi.fil[,subset(barcodes,!barcodes$cell.group %in% c(8,10)) %>% row.names]
# save(cds.epi.fil2,file = "cds_epi.rds")

# g.epi <- plot_cells(cds = cds.epi.fil2,color_cells_by = "healing",label_cell_groups = F,alpha = 0.5) +
#         scale_colour_manual(values = c("red","orange","yellow","cyan","blue","green")) +
#         theme_classic(base_size = 16)
# print(g.epi)
# ggsave(filename = "plot_umap_epithelium_filtered.png",plot = g.epi)

## computes Ly6a and Birc5 cell abundance
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Ly6a"),col.name = "Expr.Ly6a",use.raw = T,use.norm = F)
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Birc5"),col.name = "Expr.Birc5",use.raw = T,use.norm = F)
# print(hist(log10(colData(cds.epi.fil2)$Expr.Ly6a+0.1)))
# threshold.Ly6a <- summary(colData(cds.epi.fil2)$Expr.Ly6a)[5]
# colData(cds.epi.fil2)$status.Ly6a <- sapply(X = colData(cds.epi.fil2)$Expr.Ly6a,FUN = function(x) {
#   if (x > threshold.Ly6a) "hi"
#   else "lo"
# })
# g.Ly6a.gating <- plot_cells(cds = cds.epi.fil2,color_cells_by = "status.Ly6a",label_cell_groups = F) +
#   theme_classic(base_size = 16)
# ggsave(filename = "plot_epithelium_umap_Ly6a_gating.png",plot = g.Ly6a.gating)
# g.Ly6a.expression <- plot_cells(cds = cds.epi.fil2,genes = c("Ly6a"),label_cell_groups = F) +
#   theme_classic(base_size = 16)
# ggsave(filename = "plot_epithelium_umap_Ly6a_expression.png",plot = g.Ly6a.expression)
# tab <- table(colData(cds.epi.fil2)$sample,colData(cds.epi.fil2)$status.Ly6a)
# tab10 <- sweep(x = tab,MARGIN = 1,STATS = rowSums(tab),FUN = "/") * 100
# df.Ly6a.abund <- data.frame(sample=row.names(tab10),pct.Ly6a.hi=tab10[,"hi"])
# meta <- read_xlsx(path = "metadata_healing.xlsx",col_names = T)
# df.Ly6a.abund10 <- merge(x = df.Ly6a.abund,y = meta,by = "sample",all = F)
# df.Ly6a.abund10$healing <- factor(df.Ly6a.abund10$healing,levels=c("none","early","middle","late","very-late"))
# summ.Ly6a <- df.Ly6a.abund10 %>% group_by(healing) %>% summarize(med.Ly6a = median(pct.Ly6a.hi))
# g.Ly6a.abund <- ggplot(data = df.Ly6a.abund10,mapping = aes(x = healing,y = pct.Ly6a.hi)) +
#   geom_point(size = 5) +
#   geom_line(mapping = aes(y = med.Ly6a,group = 1),data = summ.Ly6a,colour = "red", size=4,alpha=0.6) +
#   theme_classic(base_size = 24)
# print(g.Ly6a.abund)
# ggsave(filename = "plot_Ly6a_abundance.svg",plot = g.Ly6a.abund)

## computes the representation of Ly6a-hi cells in each cluster
# tab <- table(colData(cds.epi.fil2)$annotation,colData(cds.epi.fil2)$status.Ly6a)
# tab10 <- sweep(x = tab,MARGIN = 1,STATS = rowSums(tab),FUN = "/") * 100
# df.Ly6a.ct <- data.frame(cluster=row.names(tab10),Ly6a.hi=tab10[,"hi"],Ly6a.lo=tab10[,"lo"])
# df.Ly6a.ct10 <- melt(data = df.Ly6a.ct,id.vars = c(1),variable.name = "status.Ly6a",value.name = "abundance")
# df.Ly6a.ct10$status.Ly6a <- factor(df.Ly6a.ct10$status.Ly6a,levels=c("Ly6a.lo","Ly6a.hi"))
# df.Ly6a.ct10$cluster <- factor(df.Ly6a.ct10$cluster,levels=c("Prog-Gpx2",
#                                                              "Prog-Gsdmc4",
#                                                              "Diff-Ly6g",
#                                                              "Diff-Dpep1/Guca2a/Mgat4c",
#                                                              "Goblet",
#                                                              "EEC",
#                                                              "Tuft"))
# g.Ly6a.ct <- ggplot(data = df.Ly6a.ct10,mapping = aes(x = cluster,y = abundance)) +
#   geom_col(mapping = aes(fill = status.Ly6a)) +
#   geom_segment(mapping = aes(x=0.5,xend=7.5,y=25,yend=25),colour = "green") +
#   scale_fill_colorblind() +
#   theme_classic(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# ggsave(filename = "plot_bar_Ly6a_by_cluster.svg",plot = g.Ly6a.ct)

## computes cell-type diversity for each healing period
# tab <- table(colData(cds.epi.fil2)$annotation,colData(cds.epi.fil2)$sample)
# tab10 <- sweep(x = tab,MARGIN = 2,STATS = colSums(tab),FUN = "/") * 100
# div.df <- -colSums(apply(X = tab,MARGIN = 2,FUN = function(x) {
#   (x/(sum(x)))*log(0.001+x/(sum(x)))
# })) %>% data.frame
# colnames(div.df) <- c("shannon")
# div.df$sample <- row.names(div.df)
# div.df$simpson <- 1/colSums(apply(X = tab,MARGIN = 2,FUN = function(x) {
#   (x/sum(x))^2
# }))
# meta <- read_xlsx(path = "metadata_healing.xlsx",col_names = T)
# div.df10 <- merge(x = div.df,y = meta,by = "sample",all = F)
# div.df10$healing <- factor(div.df10$healing,levels=c("none","early","middle","late","very-late"))
# g.simpson <- ggplot(data = div.df10,mapping = aes(x = healing,y = simpson)) +
#   geom_point(size = 5) +
#   theme_classic(base_size = 24)
# print(g.simpson)
# ggsave(filename = "plot_diversity_simpson.svg",plot = g.simpson)
# g.shannon <- ggplot(data = div.df10,mapping = aes(x = healing,y = shannon)) +
#   geom_point(size = 5) +
#   theme_classic(base_size = 24)
# print(g.shannon)
# ggsave(filename = "plot_diversity_shannon.svg",plot = g.shannon)

## computes differential gene abundance of Ly6a-hi vs. Ly6a-lo in each of the non-secretory clusters
# barcodes <- colData(cds.epi.fil2) %>% data.frame
# cds.stem <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Prog-Gpx2") %>% row.names]
# agg.stem <- aggregate.expression.by.factor(cds = cds.stem,grouping = "status.Ly6a",do.print = T,use.norm = F)
# write.table(x = agg.stem,file = "table_Ly6a_hi_v_lo_stem.csv",sep = ",",row.names = T)

# cds.prog <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Prog-Gsdmc4") %>% row.names]
# agg.prog <- aggregate.expression.by.factor(cds = cds.prog,grouping = "status.Ly6a",do.print = T,use.norm = F)
# write.table(x = agg.prog,file = "table_Ly6a_hi_v_lo_prog.csv",sep = ",",row.names = T)

# cds.diff1 <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Diff-Ly6g") %>% row.names]
# agg.diff1 <- aggregate.expression.by.factor(cds = cds.diff1,grouping = "status.Ly6a",do.print = T,use.norm = F)
# write.table(x = agg.diff1,file = "table_Ly6a_hi_v_lo_diff1.csv",sep = ",",row.names = T)

# cds.diff2 <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Diff-Dpep1/Guca2a/Mgat4c") %>% row.names]
# agg.diff2 <- aggregate.expression.by.factor(cds = cds.diff2,grouping = "status.Ly6a",do.print = T,use.norm = F)
# write.table(x = agg.diff2,file = "table_Ly6a_hi_v_lo_diff2.csv",sep = ",",row.names = T)

# agg.merge <- cbind(stem.hi = agg.stem[,1], stem.lo = agg.stem[,2],
#                                        prog.hi = agg.prog[,1], prog.lo = agg.prog[,2],
#                                        diff1.hi = agg.diff1[,1],diff1.lo = agg.diff1[,2],
#                                        diff2.hi = agg.diff2[,1],diff2.lo = agg.diff2[,2])
# row.names(agg.merge) <- row.names(agg.stem)
# agg.merge10 <- agg.merge[which(rowSums(agg.merge) > 0.5),]
# agg.ratio <- data.frame(stem = (agg.merge10[,"stem.hi"]+0.01) / (agg.merge10[,"stem.lo"]+0.01),
#                          prog = (agg.merge10[,"prog.hi"]+0.01) / (agg.merge10[,"prog.lo"]+0.01),
#                          diff1 = (agg.merge10[,"diff1.hi"]+0.01) / (agg.merge10[,"diff1.lo"]+0.01),
#                          diff2 = (agg.merge10[,"diff2.hi"]+0.01) / (agg.merge10[,"diff2.lo"]+0.01)) %>%
#                           log2 %>% data.matrix
# agg.spec <- apply(X = agg.ratio,MARGIN = 1,FUN = function(x) {
#   # return(x[which.max(abs(x))] - median(x))
#   # return(max(abs(x)) - median(abs(x)))
#   return(var(x))
# })
# agg.which <- apply(X = agg.ratio,MARGIN = 1,FUN = function(x) {
#   return(which.max(abs(x)))
# })
# agg.mag <- apply(X = agg.ratio,MARGIN = 1,FUN = function(x) {
#   return(x[which.max(abs(x))])
# })
# df.summ <- data.frame(id = names(agg.spec),pop = agg.which,
#                       log2ratio = agg.mag,specificity = agg.spec)
# feats <- rowData(x = cds.epi.fil2) %>% data.frame
# df.summ10 <- merge(x = feats,y = df.summ,by = "id",all = F)
# df.summ10$score <- abs(df.summ10$log2ratio) + abs(df.summ10$specificity)
# df.summ20 <- df.summ10[order(-df.summ10$pop,df.summ10$specificity,decreasing = T),]
# write.table(x = df.summ20,file = "table_Ly6ahi_v_lo_specificity.csv",sep = ",",row.names = F)

## looks just at genes expressed highest in individual clusters
# agg.spec.simp <- apply(X = agg.merge10,MARGIN = 1,FUN = function(x) {
#   x.ord <- x[order(x,decreasing = T)]
#   return(log2(x.ord[1]/x.ord[2]))
# })
# agg.which.simp <- apply(X = agg.merge10,MARGIN = 1,FUN = function(x) {
#   return(which.max(x))
# })
# df.summa <- data.frame(id = names(agg.which.simp),pop = agg.which.simp,
#                        specificity = agg.spec.simp)
# df.summa$pop.name <- colnames(agg.merge10)[df.summa$pop]
# feats <- rowData(x = cds.epi.fil2) %>% data.frame
# df.summa10 <- merge(x = feats,y = df.summa,by = "id",all = F)
# df.summa20 <- df.summa10[order(-df.summa10$pop,df.summa10$specificity,decreasing = T),]
# print(plot_genes_by_group(cds = cds.epi.central,markers = df.summa20[1:40,"gene_short_name"],group_cells_by = "condition",max.size = 4))
# write.table(x = df.summa20,file = "table_Ly6ahi_v_lo_specificity2.csv",sep = ",",row.names = F)
# df.summa30 <- subset(df.summa20,df.summa20$specificity > 1.5)
# mat.summa <- agg.merge10[df.summa30$id,]
# row.names(mat.summa) <- df.summa30$gene_short_name
# mat.summa10 <- mat.summa[!duplicated(row.names(mat.summa)),]
# mat.summa20 <- t(apply(X = mat.summa10,MARGIN = 1,FUN = function(x) {
#   (x - mean(x)) / sd(x)
# }))
# pheatmap(mat = mat.summa20,cluster_cols = F)

## tries a different way - cluster Ly6a-hi epithelium separately
## and then ask how to distinguish the clusters
# barcodes <- colData(cds.epi.fil2) %>% data.frame
# cds.Ly6ahi <- cds.epi.fil2[,subset(barcodes,barcodes$status.Ly6a=="hi") %>% row.names]
# cds.Ly6ahi <- reduce_dimension(cds = cds.Ly6ahi,verbose = T)
# cds.Ly6ahi <- cluster_cells(cds = cds.Ly6ahi)
# print(plot_cells(cds = cds.Ly6ahi))
# mtr <- top_markers(cds = cds.Ly6ahi,genes_to_test_per_group = 100,verbose = T)

## or do some subclustering and then evaluate Ly6a-hi progenitors
# barcodes <- colData(cds.epi.fil2) %>% data.frame
# cds.epi.central <- cds.epi.fil2[,subset(barcodes,barcodes$annotation %in% c("Prog-Gpx2","Prog-Gsdmc4","Diff-Ly6g","Diff-Dpep1/Guca2a/Mgat4c")) %>% row.names]
# cds.epi.central <- reduce_dimension(cds = cds.epi.central,verbose = T)
# cds.epi.central <- cluster_cells(cds = cds.epi.central)
# print(plot_cells(cds = cds.epi.central,group_label_size = 6))
# colData(cds.epi.central)$condition <- paste(colData(cds.epi.central)$annotation,colData(cds.epi.central)$status.Ly6a,sep = "-")
# mtr.central <- top_markers(cds = cds.epi.central,group_cells_by = "condition",genes_to_test_per_group = 100,verbose = T)
# goi <- subset(mtr.central,mtr.central$cell_group=="Prog-Gpx2-hi")$gene_short_name
# print(plot_genes_by_group(cds = cds.epi.central,markers = filter.marker.results(mtr.central,5,"specificity"),group_cells_by = "condition"))
# print(plot_genes_by_group(cds = cds.epi.central,markers = goi,group_cells_by = "condition"))

## computes time evolution of total epithelium
# agg.healing.all <- aggregate.expression.by.factor(cds = cds.epi.fil2,grouping = "healing",do.print = T,use.norm = F)
# feats <- rowData(x = cds.epi.fil2) %>% data.frame
# df.healing <- data.frame(agg.healing.all)
# df.healing$id <- row.names(df.healing)
# df.healing10 <- merge(x = feats,y = df.healing,by = "id",all = F)
# df.healing20 <- df.healing10[-dim(df.healing10)[[2]]]
# write.table(x = df.healing20,file = "table_healing.csv",sep = ",",row.names = F)
# mat.healing <- data.matrix(df.healing20[,3:7])
# row.names(mat.healing) <- df.healing20$gene_short_name
# mat.healing10 <- mat.healing[which(rowSums(mat.healing)>0.5),]
# mat.healing20 <- sweep(x = (mat.healing10+0.01),MARGIN = 1,STATS = mat.healing10[,1]+0.01,FUN = "/")
# mat.healing30 <- log2(mat.healing20)
# mat.healing40 <- mat.healing30[apply(X = mat.healing30,MARGIN = 1,FUN = function(x) {
#   if (max(abs(x)) > 1.5) { T }
#   else { F }
# }),]
# mat.healing50 <- t(apply(X = mat.healing40,MARGIN = 1,FUN = function(x) {
#   return((x - min(x)) / (max(x) - min(x)))
# }))
# mat.healing60 <- mat.healing50[,c("none","early","middle","late","very.late")]
# mat.healing70 <- mat.healing60[!duplicated(row.names(mat.healing60)),]
# hm <- heatmap.2(x = mat.healing60,Colv = F,dendrogram = "row",tracecol = "blue")
# ph <- pheatmap(mat = mat.healing70,cluster_cols = F)
# cutter <- cutree(tree = ph$tree_row,k = 7)
# anno.df <- data.frame(type = cutter)
# row.names(anno.df) <- row.names(mat.healing70)
# anno.df$type <- as.factor(anno.df$type)
# pheatmap(mat = mat.healing70,cluster_cols = F,annotation_row = anno.df)
# anno.df10 <- anno.df
# anno.df10$gene_short_name <- row.names(anno.df10)
# write.table(x = anno.df10,file = "table_epithelium_time_gene_modules.csv",sep = ",",row.names = F)

## looks to see if FPC genes correspond to what type of time evolution
# df.fpc <- merge(x = anno.df10,y = df.summa30,by = "gene_short_name",all.x = T,all.y = F)
# tab <- table(df.fpc$type,df.fpc$pop.name)
# tab.df <- data.frame(tab)
# colnames(tab.df) <- c("type","pop.name","count")
# k <- ggplot(data = tab.df,mapping = aes(x = pop.name,y = count)) +
#   geom_col(mapping = aes(fill = type))
# print(k)
# df.fpc10 <- df.fpc[-which(df.fpc$gene_short_name=="conf"),]
# row.names(df.fpc10) <- df.fpc10$gene_short_name
# df.fpc20 <- df.fpc10[,c("type","pop.name")]
# df.fpc20$type <- as.factor(df.fpc20$type)
# pheatmap(mat = mat.healing70,cluster_cols = F,annotation_row = df.fpc20)
# write.table(x = df.fpc10,file = "table_epithelium_time_gene_modules_FPC.csv",sep = ",",row.names = F)

## plots representative dotted plots
# goi <- c("Ccna2","Icam1","Sulf2","Krt17","Gpx2","Clu","Ass1","Pyy","Scnn1b","Ly6a","Lgr5")
# cds.sub <- cds.epi.fil2[map.gene.short.names(cds = cds.epi.fil2,gene.short.names = goi),]
# colData(cds.sub)$healing <- factor(colData(cds.sub)$healing,levels = c("none","early","middle","late","very-late"))
# h <- plot.expression.with.dots(cds.subset = cds.sub,grouping = "healing",yjitter = 0.2,logticks = T)
# ggsave(filename = "plot_dots_genes_time.svg",plot = h)

## does a different version of analysis using gene regression
# cds.sub <- cds.epi.central[row.names(agg.merge10),]
# gf <- fit_models(cds = cds.sub,model_formula_str = "~condition",verbose = T)
# ct <- coefficient_table(model_tbl = gf)
# ct10 <- filter(ct,q_value < 0.001 & normalized_effect > 3)
# ct20 <- ct10[,c(1,2,8,12,13,15)] %>% data.frame
# write.table(x = ct20,file = "table_Ly6a_diffexp.csv",sep = ",",row.names = F)

# save(cds.epi.central,file = "cds_epi_no_secretory.rds")

## checks cell counts
# tab <- table(colData(cds.epi.fil2)$annotation,colData(cds.epi.fil2)$healing)
# tab10 <- sweep(x = tab,MARGIN = 2,STATS = colSums(tab),FUN = "/") * 100
# tab20 <- tab10[,-6]
# df <- data.frame(tab20)
# colnames(df) <- c("cluster","healing","frequency")
# df$cell.type <- df$cluster
# df$healing <- factor(df$healing,levels = c("none","early","middle","late","very-late","Il10"))
# h.epi <- ggplot(data = df,mapping = aes(x = healing,y = frequency)) +
#   geom_col(mapping = aes(fill = cell.type)) +
#   scale_fill_calc() +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(h.epi)
# ggsave(filename = "plot_abundance_epithelium.svg",plot = h.epi)