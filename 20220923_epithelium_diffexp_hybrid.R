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
library(randomForest)
library(ggrepel)
source(file = "camMonocleHelper.R")

# load(file = "cds_epi.rds")

## makes simple plot of clusters on UMAP
# g <- plot_cells(cds = cds.epi.fil2,color_cells_by = "annotation",show_trajectory_graph = F,label_cell_groups = F)
# print(g)
# ggsave(filename = "plot_epithelim_umap_celltypes.png",plot = g)

## makes plot of markers of each cluster


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
# summ <- div.df10 %>% group_by(healing) %>% summarize(m.simpson = mean(simpson))
# g.simpson <- ggplot(data = div.df10,mapping = aes(x = healing,y = simpson)) +
#   geom_point(size = 5) +
#   geom_line(mapping = aes(y = m.simpson,group = 1),data = summ,size = 3,colour = "red",alpha=0.7) +
#   expand_limits(y = 0) +
#   coord_fixed(ratio = 1) +
#   theme_classic(base_size = 24)
# print(g.simpson)
# ggsave(filename = "plot_diversity_simpson.svg",plot = g.simpson)


## computes Ly6a and Birc5 cell abundance and gates based on raw values
## then performs cluster specific feature identification on normalized values
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Ly6a"),col.name = "Expr.Ly6a",use.raw = T,use.norm = F)
# print(hist(colData(cds.epi.fil2)$Expr.Ly6a))
# threshold.Ly6a <- mean(c(summary(colData(cds.epi.fil2)$Expr.Ly6a)[3],summary(colData(cds.epi.fil2)$Expr.Ly6a)[5]))
# colData(cds.epi.fil2)$status.Ly6a <- sapply(X = colData(cds.epi.fil2)$Expr.Ly6a,FUN = function(x) {
#   if (x > threshold.Ly6a) "hi"
#   else "lo"
# })
# g.Ly6a.gating <- plot_cells(cds = cds.epi.fil2,color_cells_by = "status.Ly6a",label_cell_groups = F) +
#   theme_classic(base_size = 16)
# print(g.Ly6a.gating)
# ggsave(filename = "plot_epithelium_umap_Ly6a_gating_hybrid.png",plot = g.Ly6a.gating)
# tab <- table(colData(cds.epi.fil2)$sample,colData(cds.epi.fil2)$status.Ly6a)
# tab10 <- sweep(x = tab,MARGIN = 1,STATS = rowSums(tab),FUN = "/") * 100
# df.Ly6a.abund <- data.frame(sample=row.names(tab10),pct.Ly6a.hi=tab10[,"hi"])
# meta <- read_xlsx(path = "metadata_healing.xlsx",col_names = T)
# df.Ly6a.abund10 <- merge(x = df.Ly6a.abund,y = meta,by = "sample",all = F)
# df.Ly6a.abund10$healing <- factor(df.Ly6a.abund10$healing,levels=c("none","early","middle","late","very-late"))
# summ <- df.Ly6a.abund10 %>% group_by(healing) %>% summarize(med.pct.Ly6a.hi = median(pct.Ly6a.hi)) %>% data.frame
# g.Ly6a.abund <- ggplot(data = df.Ly6a.abund10,mapping = aes(x = healing)) +
#   geom_point(size = 5,alpha=0.8,mapping = aes(y = pct.Ly6a.hi)) +
#   geom_line(data = summ,mapping = aes(y = med.pct.Ly6a.hi,group=1),colour = "red",size=4,alpha=0.5) +
#   theme_classic(base_size = 24)
# print(g.Ly6a.abund)
# ggsave(filename = "plot_Ly6a_abundance_hybrid.svg",plot = g.Ly6a.abund)

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
#   geom_segment(mapping = aes(x=0.5,xend=7.5,y=31.8,yend=31.8),colour = "green") +
#   scale_fill_colorblind() +
#   theme_classic(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(g.Ly6a.ct)
# ggsave(filename = "plot_bar_Ly6a_by_cluster_hybrid.svg",plot = g.Ly6a.ct)

## computes the representation of Ly6a-hi cells in each cluster - RAW
# tab <- table(colData(cds.epi.fil2)$annotation,colData(cds.epi.fil2)$status.Ly6a)
# tab10 <- tab
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
# g.Ly6a.ct.raw <- ggplot(data = df.Ly6a.ct10,mapping = aes(x = cluster,y = abundance)) +
#   geom_col(mapping = aes(fill = status.Ly6a)) +
#   geom_segment(mapping = aes(x=0.5,xend=7.5,y=400,yend=400),colour = "green") +
#   scale_fill_colorblind() +
#   theme_classic(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(g.Ly6a.ct.raw)
# ggsave(filename = "plot_bar_Ly6a_by_cluster_raw_hybrid.svg",plot = g.Ly6a.ct.raw)


## computes differential gene abundance of Ly6a-hi vs. Ly6a-lo in each of the clusters that are big enough
# barcodes <- colData(cds.epi.fil2) %>% data.frame
# 
# cds.stem <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Prog-Gpx2") %>% row.names]
# agg.stem <- aggregate.expression.by.factor(cds = cds.stem,grouping = "status.Ly6a",do.print = T)
# write.table(x = agg.stem,file = "table_Ly6a_hi_v_lo_stem_hybrid.csv",sep = ",",row.names = T)
# 
# cds.prog <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Prog-Gsdmc4") %>% row.names]
# agg.prog <- aggregate.expression.by.factor(cds = cds.prog,grouping = "status.Ly6a",do.print = T)
# write.table(x = agg.prog,file = "table_Ly6a_hi_v_lo_prog_hybrid.csv",sep = ",",row.names = T)
# 
# cds.diff1 <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Diff-Ly6g") %>% row.names]
# agg.diff1 <- aggregate.expression.by.factor(cds = cds.diff1,grouping = "status.Ly6a",do.print = T)
# write.table(x = agg.diff1,file = "table_Ly6a_hi_v_lo_diff1_hybrid.csv",sep = ",",row.names = T)
# 
# cds.diff2 <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Diff-Dpep1/Guca2a/Mgat4c") %>% row.names]
# agg.diff2 <- aggregate.expression.by.factor(cds = cds.diff2,grouping = "status.Ly6a",do.print = T)
# write.table(x = agg.diff2,file = "table_Ly6a_hi_v_lo_diff2_hybrid.csv",sep = ",",row.names = T)
# 
# cds.gob <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Goblet") %>% row.names]
# agg.gob <- aggregate.expression.by.factor(cds = cds.gob,grouping = "status.Ly6a",do.print = T)
# write.table(x = agg.gob,file = "table_Ly6a_hi_v_lo_gob_hybrid.csv",sep = ",",row.names = T)
# 
# agg.merge <- cbind(stem.hi = agg.stem[,1], stem.lo = agg.stem[,2],
#                                        prog.hi = agg.prog[,1], prog.lo = agg.prog[,2],
#                                        diff1.hi = agg.diff1[,1],diff1.lo = agg.diff1[,2],
#                                        diff2.hi = agg.diff2[,1],diff2.lo = agg.diff2[,2],
#                                        gob.hi = agg.gob[,1],gob.lo = agg.gob[,2])
# row.names(agg.merge) <- row.names(agg.stem)
# agg.merge10 <- agg.merge[which(rowSums(agg.merge) > 0.05),]
# agg.ratio <- apply(X = agg.merge10,MARGIN = 1,FUN = function(x) {
#   vec <- c((x[1]-x[2])/(x[2]+0.001),
#            (x[3]-x[4])/(x[4]+0.001),
#            (x[5]-x[6])/(x[6]+0.001),
#            (x[7]-x[8])/(x[8]+0.001),
#            (x[9]-x[10])/(x[10]+0.001))
#   return(vec)
# }) %>% t
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
#                       fold_change = agg.mag,specificity = agg.spec)
# feats <- rowData(x = cds.epi.fil2) %>% data.frame
# df.summ10 <- merge(x = feats,y = df.summ,by = "id",all = F)
# df.summ10$score <- abs(df.summ10$fold_change) + abs(df.summ10$specificity)
# df.summ20 <- df.summ10[order(-df.summ10$pop,df.summ10$specificity,decreasing = T),]
# write.table(x = df.summ20,file = "table_Ly6ahi_v_lo_specificity_hybrid.csv",sep = ",",row.names = F)

## looks just at genes expressed highest in individual clusters
# agg.spec.simp <- apply(X = agg.merge10,MARGIN = 1,FUN = function(x) {
#   x.ord <- x[order(x,decreasing = T)]
#   return((x.ord[1]-x.ord[2])/(median(x)+0.001))
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
# write.table(x = df.summa20,file = "table_Ly6ahi_v_lo_specificity2_hybrid.csv",sep = ",",row.names = F)
# df.summa30 <- subset(df.summa20,df.summa20$specificity > 1.5)
# mat.summa <- agg.merge10[df.summa30$id,]
# row.names(mat.summa) <- df.summa30$gene_short_name
# mat.summa10 <- mat.summa[!duplicated(row.names(mat.summa)),]
# mat.summa20 <- t(apply(X = mat.summa10,MARGIN = 1,FUN = function(x) {
#   (x - mean(x)) / sd(x)
# }))
# pheatmap(mat = mat.summa20,cluster_cols = F)

## loads pre-computed time evolution matrix
# df.healing20 <- read.csv(file = "table_healing.csv",header = T)
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

## looks to see if FPC genes correspond to what type of healing time evolution
## remove type 6 as these are likely artifacts of squamous cell abundance at very-late healing
# anno.df10 <- read.csv(file = "table_epithelium_time_gene_modules.csv",header = T)
# anno.df20 <- subset(anno.df10,!anno.df10$type %in% c(6))
# anno.df20$new.type <- sapply(X = anno.df20$type,FUN = function(x) {
#   if (x == 7) { 6 }
#   else { x }
# })
# df.fpc <- merge(x = anno.df20,y = df.summa30,by = "gene_short_name",all.x = T,all.y = F)
# tab <- table(df.fpc$new.type,df.fpc$pop.name)
# tab.df <- data.frame(tab)
# colnames(tab.df) <- c("new.type","pop.name","count")
# k <- ggplot(data = tab.df,mapping = aes(x = pop.name,y = count)) +
#   geom_col(mapping = aes(fill = new.type))
# print(k)
# ggsave(filename = "plot_bar_epithelium_gene_type_by_cluster_hybrid.svg",plot = k)
# df.fpc10 <- df.fpc[-which(df.fpc$gene_short_name=="conf"),]
# row.names(df.fpc10) <- df.fpc10$gene_short_name
# df.fpc20 <- df.fpc10[,c("new.type","pop.name")]
# df.fpc20$new.type <- as.factor(df.fpc20$new.type)
# mat.healing80 <- mat.healing70[row.names(df.fpc20),]
# df.fpc30 <- df.fpc20
# df.fpc30$pop.name <- sapply(X = df.fpc30$pop.name,FUN = function(txt) {
#   if (txt %in% c("diff1.hi","stem.hi","stem.lo")) { txt }
#   else { NA }
# })
# annoCol <- list(pop.name=c(diff1.hi = "blue",stem.hi = "red",stem.lo = "green"),
#                 new.type=c(`1` = "red",`2` = "orange",`3` = "yellow",`4` = "green",`5` = "blue",`6` = "purple"))
# pheatmap(mat = mat.healing80,cluster_cols = F,annotation_row = df.fpc30,annotation_colors = annoCol)
# write.table(x = df.fpc10,file = "table_epithelium_time_gene_modules_FPC_hybrid.csv",sep = ",",row.names = F)

## analyzes proliferative cells in Ly6a-hi vs. lo groups
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Birc5"),col.name = "Expr.Birc5",use.raw = T,use.norm = F)
# print(hist(colData(cds.epi.fil2)$Expr.Birc5))
# bcvec <- colData(cds.epi.fil2)$Expr.Birc5[order(colData(cds.epi.fil2)$Expr.Birc5,decreasing = F)]
# threshold.Birc5 <- bcvec[floor(length(bcvec)*0.95),drop = T]
# colData(cds.epi.fil2)$status.Birc5 <- sapply(X = colData(cds.epi.fil2)$Expr.Birc5,FUN = function(x) {
#   if (x > threshold.Birc5) "prolif"
#   else "non-prolif"
# })
# colData(cds.epi.fil2)$healing <- factor(colData(cds.epi.fil2)$healing)
# tab.birc.df <- table(colData(cds.epi.fil2)$status.Birc5,colData(cds.epi.fil2)$status.Ly6a,colData(cds.epi.fil2)$healing) %>% data.frame
# colnames(tab.birc.df) <- c("Birc5","Ly6a","Healing","Freq")
# g.birc <- ggplot(data = tab.birc.df,mapping = aes(x = Ly6a,y = Freq)) +
#   facet_wrap(facets = .~Healing) +
#   geom_col(mapping = aes(fill = Birc5))
# print(g.birc)
# summ.birc <- split(x = tab.birc.df,f = list(tab.birc.df$Healing),drop = T) %>%
#   lapply(FUN = function(df) {
#     subdf <- df[which(df$Birc5 == "prolif"),]
#     cell.sum <- sum(subdf$Freq)
#     retdf <- subdf
#     retdf$Prop <- retdf$Freq / cell.sum * 100
#     norm.factor <- cell.sum / sum(df$Freq)
#     retdf$Norm.Prop <- retdf$Prop * norm.factor
#     retdf
#   }) %>% do.call(what = "rbind",args = .)
# summ.birc$Healing <- factor(summ.birc$Healing,levels=c("none","early","middle","late","very-late"))
# summ.birc$Ly6a <- factor(summ.birc$Ly6a,levels = c("lo","hi"))
# g10.birc <- ggplot(data = summ.birc,mapping = aes(x = Healing,y = Norm.Prop)) +
#   geom_col(mapping = aes(fill = Ly6a)) +
#   scale_fill_colorblind() +
#   ylab("% of Birc5+ cells") +
#   theme_classic(base_size = 24)
# print(g10.birc)
# ggsave(filename = "plot_bar_Birc5_by_Ly6a.svg",plot = g10.birc)

## analyzes predictors of Ly6a+ epithelium
# bcrf <- colData(cds.epi.fil2) %>% data.frame
# agg.hivlo <- aggregate.expression.by.factor(cds = cds.epi.fil2,grouping = "status.Ly6a",do.print = T,use.norm = T)
# feats <- rowData(x = cds.epi.fil2) %>% data.frame
# aggdf.hivlo <- data.frame(agg.hivlo)
# aggdf.hivlo$id <- row.names(aggdf.hivlo)
# aggdf10.hivlo <- merge(x = feats,y = aggdf.hivlo, by = "id",all = F)
# aggdf10.hivlo$sum <- aggdf10.hivlo$hi + aggdf10.hivlo$lo
# write.table(x = aggdf10.hivlo,file = "table_Ly6a_hi_v_lo.csv",sep = ",",row.names = F)
# aggdf20.hivlo <- subset(aggdf10.hivlo,aggdf10.hivlo$sum > 0.3)
# aggdf20.hivlo$FC <- (aggdf20.hivlo$hi - aggdf20.hivlo$lo) / aggdf20.hivlo$sum
# aggdf20.hivlo$product <- aggdf20.hivlo$FC*aggdf20.hivlo$sum
# aggdf20.hivlo$highlight <- ifelse(aggdf20.hivlo$FC*aggdf20.hivlo$sum > 0.18,"signature","background")
# aggdf20.hivlo$label <- mapply(FUN = function(x,y,txt) {
#   if (x*y > 0.192) { txt } 
#   else { NA }
# },aggdf20.hivlo$FC,aggdf20.hivlo$sum,aggdf20.hivlo$gene_short_name)
# write.table(x = aggdf20.hivlo,file = "table_Ly6a_hi_v_lo_thresholded.csv",sep = ",",row.names = F)
# g.repro <- ggplot(data = aggdf20.hivlo,mapping = aes(x = FC,y = sum)) +
#   geom_point(mapping = aes(colour = highlight))
# print(g.repro)

## applies a stringent criterion to determine FPC genes
## have to be temporally elevated in "group 4" pattern as this group has no expression at baseline
# df.fpc10 <- read.csv(file = "table_epithelium_time_gene_modules_FPC_hybrid.csv",header = T)
# df.fpc.a <- subset(df.fpc10,df.fpc10$pop.name == "stem.hi" & df.fpc10$new.type == 4)
# write.table(x = df.fpc.a,file = "table_FPC_signature_conservative.csv",sep = ",",row.names = F)

## plots representative dotted plots
# goi <- c("Birc5","Ifitm1","Sulf2","Krt17","Gpx2","Clu","Ass1","Ido1","Ly6a","Lgr5")
# cds.sub <- cds.epi.fil2[map.gene.short.names(cds = cds.epi.fil2,gene.short.names = goi),]
# colData(cds.sub)$healing <- factor(colData(cds.sub)$healing,levels = c("none","early","middle","late","very-late"))
# h <- plot.expression.with.dots(cds.subset = cds.sub,grouping = "healing",yjitter = 0.2,logticks = F)
# print(h)
# ggsave(filename = "plot_dots_genes_time.svg",plot = h)
