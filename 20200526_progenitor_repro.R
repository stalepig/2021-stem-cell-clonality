setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
library(dplyr)
library(scales)
library(umap)
source("camMonocleHelper.R")

# load(file = "cds_epi_filtered2.rds",verbose = T)

## pulls out cluster 2 - progenitor cell cluster
# barcodes <- colData(cds.epi.fil.fil) %>% data.frame
# cds.progenitor <- cds.epi.fil.fil[,subset(barcodes,barcodes$cell.group==2) %>% row.names]
# cds.progenitor <- recluster.cds(cds = cds.progenitor)
# cds.progenitor <- choose_cells(cds = cds.progenitor)
# cds.progenitor <- recluster.cds(cds = cds.progenitor)

## performs differential expression testing
# gene_fits <- fit_models(cds = cds.progenitor,model_formula_str = "~sample.name",verbose = T)
# gene_fits10 <- subset(gene_fits,status == "OK")
# fit_coeffs <- coefficient_table(gene_fits10)
# fit_coeffs10 <- subset(fit_coeffs,!fit_coeffs$term %in% c("(Intercept)"))
# fit_coeffs20 <- fit_coeffs10[,-c(3:4)]
# fit_coeffs30 <- subset(fit_coeffs20,fit_coeffs20$p_value<0.001)
# fit_coeffs40 <- fit_coeffs30[order(fit_coeffs30$term,fit_coeffs30$normalized_effect),]

## does a more sensitive version where all the dss samples are pooled together
# gene_fits.ms <- fit_models(cds = cds.progenitor,model_formula_str = "~injury",verbose = T)
# fit_coeffs.ms <- subset(gene_fits.ms,gene_fits.ms$status=="OK") %>% coefficient_table
# fit_coeffs.ms10 <- fit_coeffs.ms[,-c(3:4)]
# fit_coeffs.ms20 <- subset(fit_coeffs.ms10,!fit_coeffs.ms10$term %in% c("(Intercept)"))
# fit_coeffs.ms30 <- subset(fit_coeffs.ms20,fit_coeffs.ms20$p_value < 0.001)
# fit_coeffs.ms40 <- fit_coeffs.ms30[order(fit_coeffs.ms30$normalized_effect),]

## performs Moran's I
# graph_test_res <- graph_test(cds = cds.progenitor,neighbor_graph = "principal_graph",verbose = T)
# graph_test_res10 <- filter(graph_test_res,status=="OK") %>%
#   filter(log10(q_value) < -50) %>%
#   filter(abs(morans_I) > 0.1)
# graph_test_res20 <- graph_test_res10[order(graph_test_res10$morans_I,decreasing = T),]
# gene_modules <- find_gene_modules(cds = cds.progenitor[graph_test_res20$id,],max_components = 5,verbose = T)
# g <- plot_cells(cds = cds.progenitor,genes=graph_test_res20[81:100,"gene_short_name"],show_trajectory_graph = F)
# print(g)

# goi <- c("Birc5","Ccnb1","Mki67")
# goi <- c("Chga","Syp","Gcg")
# goi <- c("Lgr5","Ascl2","Hopx","Ephb3","Lrig1","Clu","Ephb2")
# goi <- c("Spink4","Muc2","Reg4","Kit","Atoh1")
# goi.label <- "idx.dcs"
# cds.progenitor <- add.gene.sum.column(cds = cds.progenitor,genes.of.interest = goi,col.name=goi.label)
# h <- plot_cells(cds = cds.progenitor,color_cells_by = goi.label,show_trajectory_graph = T,trajectory_graph_color = "green",trajectory_graph_segment_size = 3) +
#   geom_point(size=2,mapping=aes_string(colour=goi.label),alpha=0.7) +
#   scale_colour_gradient(low="blue",high="red")
# print(h)

## really zooms in on the "stem" cells
# k <- plot_cells(cds = cds.progenitor,color_cells_by = "idx.stem",show_trajectory_graph = T,trajectory_graph_color = "green") +
#   facet_wrap(facets = .~injury) +
#   geom_point(mapping=aes(colour=idx.stem))
# print(k)
# cds.stem <- choose_cells(cds = cds.progenitor)
# cds.stem <- recluster.cds(cds = cds.stem)

## looks for reprogramming markers by graph test
# gtr <- graph_test(cds = cds.stem,neighbor_graph = "principal_graph")
# gtr10 <- filter(gtr,status=="OK") %>%
#   filter(q_value < 0.01) %>%
#   filter(morans_I > 0.1)
# gtr20 <- gtr10[order(gtr10$morans_I,decreasing = T),]

## tried various ways to break down data - nothing satisfactory here
# exprs.subset <- data.matrix(exprs(cds.stem[gtr20$id,]))
# um <- umap(d = exprs.subset)
# pc <- prcomp(x = exprs.subset,retx = T,center = T,scale. = T)
# gm <- find_gene_modules(cds = cds.stem[gtr20$id,],max_components = 5,verbose = T,umap.fast_sgd = F,)
# print(qplot(x = gm$dim_1,y = gm$dim_2,geom = "density2d"))
# m <- data.matrix(cbind(gm$dim_1,gm$dim_2))
# centers <- list(c(-2.5,-0.1),c(0.0,-0.5),c(0.4,0.0),c(-1.0,0.3),c(4.7,-0.1),c(-3.0,2.1))
# center.rows <- lapply(X = centers,FUN = function(cent) {
#   return(find.closest.row(input.matrix = m,center.point = cent))
# }) %>% do.call("rbind",.)
# kclust <- kmeans(x = m,centers = m[center.rows,])

## uses Ly6a to "gate" on reprogramming-associated markers, and then will use these to narrow down the results from graph_test
## this linear model is not quite satisfactory, as it does seem not take into account cells that do not express the transcript
# cds.stem <- add.gene.sum.column(cds = cds.stem,genes.of.interest = c("Ly6a"),col.name = "idx.Ly6a",use.raw=T)
# colData(cds.stem)$Ly6a.status <- ifelse(log10(colData(cds.stem)$idx.Ly6a+0.1) > 0.783,"hi","lo")
# colData(cds.stem)$Ly6a.status <- factor(colData(cds.stem)$Ly6a.status,levels=c("lo","hi"))
# gene_fits.Ly6a <- fit_models(cds = cds.stem,model_formula_str = "~Ly6a.status",verbose = T)
# fit_coeffs.Ly6a <- filter(gene_fits.Ly6a,status=="OK") %>% coefficient_table()
# fit_coeffs.Ly6a10 <- fit_coeffs.Ly6a[,-c(3:4)] %>%
#   filter(p_value < 0.001) %>%
#   filter(!term %in% c("(Intercept)"))
# fit_coeffs.Ly6a20 <- fit_coeffs.Ly6a10[order(fit_coeffs.Ly6a10$normalized_effect),]
# kk <- plot_cells(cds = cds.stem,genes = c("Ly6a"),show_trajectory_graph = F,cell_size = 2)
# print(kk)

## does the Ly6a gating differential expression, but the "simple way"
# exprs.stem <- exprs(cds.stem)
# Ly6a.status <- colData(x = cds.stem)$Ly6a.status
# names(Ly6a.status) <- colData(cds.stem) %>% data.frame %>% row.names
# aggdf <- data.frame(id=vector(mode = "character",length = 0),
#                     lo=vector(mode = "numeric",length = 0),
#                     hi=vector(mode = "numeric",length = 0))
# for (i in c(1:dim(exprs.stem)[1])) {
#   if (i %% 1000 == 0) { print(i) }
#   theVec <- as.vector(exprs.stem[i,])
#   names(theVec) <- colnames(exprs.stem)
#   id <- row.names(exprs.stem)[i]
#   agg <- tapply(X = theVec,INDEX = Ly6a.status,FUN = mean)
#   aggrow <- data.frame(id=id,lo=agg[1],hi=agg[2])
#   aggdf <- rbind(aggdf,aggrow)
# }
# aggdf10 <- filter(aggdf,lo+hi>0.7)
# aggdf20 <- merge(x = aggdf10,y = data.frame(rowData(cds.stem)),by = "id",all = F)
# aggdf20$log2ratio <- log2(aggdf20$hi/aggdf20$lo)
# aggdf30 <- filter(aggdf20,abs(log2ratio)>1.75)
# aggdf30 <- aggdf30[order(aggdf30$log2ratio),]

## calculates specificity for reprogramming-associated genes
# cds.epi.reclust <- cluster_cells(cds = cds.epi.fil.fil,cluster_method = "louvain",verbose = T)
# exprs.epi <- exprs(cds.epi.reclust[aggdf30$id %>% as.character,])
# colData(cds.epi.reclust)$cell_group <- as.factor(clusters(cds.epi.reclust))
# cell_group <- colData(cds.epi.reclust)$cell_group
# names(cell_group) <- row.names(colData(cds.epi.reclust))
# aggm <- matrix(nrow = dim(exprs.epi)[1],ncol = length(levels(cell_group)))
# for (i in c(1:dim(exprs.epi)[1])) {
#   theVec <- as.vector(exprs.epi[i,])
#   names(theVec) <- colnames(exprs.epi)
# 
#   aggvec <- tapply(X = theVec,INDEX = cell_group,FUN = mean)
#   aggm[i,] <- aggvec
# }
# row.names(aggm) <- aggdf30$id
# colnames(aggm) <- as.character(c(1:length(levels(cell_group))))
# aggdf30$max.cluster <- apply(X = aggm,MARGIN = 1,FUN = which.max)
# aggdf30$min.cluster <- apply(X = aggm,MARGIN = 1,FUN = which.min)
# aggdf30$median.expression <- apply(X = aggm,MARGIN = 1,FUN = median)
# aggdf30$quar75.expression <- apply(X = aggm,MARGIN = 1,FUN = function(x) {
#   return(x[order(x)[round(x = length(x)*0.75,digits = 0)]])
# })
# aggdf30$quar90.expression <- apply(X = aggm,MARGIN = 1,FUN = function(x) {
#   return(x[order(x)[round(x = length(x)*0.9,digits = 0)]])
# })
# aggdf30$max.expression <- apply(X = aggm,MARGIN = 1,FUN = max)

## makes specificity plots
# genes.to.plot <- c("Ly6a","Clu")
# gene.ids <- aggdf20[which(aggdf20$gene_short_name %in% genes.to.plot),"id"]
# m.plot <- aggm[as.character(gene.ids),]
# df.plot <- data.frame(id = row.names(m.plot),m.plot)
# df.plot10 <- melt(data = df.plot,id.vars = c("id"))
# colnames(df.plot10) <- c("id","cell.group","abundance")
# df.plot20 <- merge(x = df.plot10,y = rowData((cds.epi.reclust)) %>% data.frame,by="id",all=F)
# oo <- ggplot(data = df.plot20,mapping = aes(x=cell.group,y=abundance)) +
#   facet_wrap(facets = .~gene_short_name,scales = "free_y") +
#   geom_bar(stat="identity") +
#   xlab("Cell cluster") +
#   ylab("Abundance (mean tpm)") +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_blank())
# print(oo)

## now restricts analysis to genes highest only in clusters 15 or 18
# aggdf40 <- filter(aggdf30,max.cluster %in% c(15,18))
# aggdf40$specificity <- log2(aggdf40$hi / aggdf40$quar90.expression)
# aggdf40 <- aggdf40[order(aggdf40$specificity,decreasing = T),]
# ll <- plot_genes_violin(cds_subset = cds.epi.reclust["ENSMUSG00000075602",],group_cells_by = "cell_group")
# ll <- plot_genes_violin(cds_subset = cds.epi.reclust["ENSMUSG00000006800",],group_cells_by = "cell_group")
# print(ll)

## now takes the specific genes and calculates the dynamics of emergence
# aggdf50 <- filter(aggdf40,specificity > 1.5)
# m.dyn <- aggregate.expression.by.factor(cds = cds.epi.reclust[aggdf50$id %>% as.character,],grouping = "sample.name")
# m.dyn.norm <- t(apply(X = m.dyn,MARGIN = 1,FUN = function(x) { (x-min(x))/(max(x)-min(x)) }))
# m.dyn.norm <- m.dyn.norm[order(m.dyn.norm[,2],decreasing=T),]
# df.dyn <- data.frame(id=row.names(m.dyn.norm),m.dyn.norm) %>% melt(id.vars="id")
# colnames(df.dyn) <- c("id","sample.name","abundance")
# df.dyn.agg <- aggregate(formula=abundance~sample.name,data=df.dyn,FUN=mean)
# df.dyn.agg$stdev <- tapply(X = df.dyn$abundance,INDEX = df.dyn$sample.name,FUN = sd)
# df.dyn.agg$stderr <- df.dyn.agg$stdev / sqrt(dim(m.dyn.norm)[1])
# df.dyn10 <- merge(x = df.dyn,y = rowData(cds.epi.reclust) %>% data.frame,by = "id",all=F)
# df.dyn20 <- dcast(data = df.dyn10,formula = id+gene_short_name~sample.name,value.var = "abundance")
# df.dyn20 <- df.dyn20[order(df.dyn20$D06,decreasing = T),]
# genes.to.highlight <- c("Ass1","Hs3st1","Clu","Cd47")
# df.dyn20$highlight <- ifelse(df.dyn20$gene_short_name %in% genes.to.highlight,"yes","no")
# df.dyn30 <- melt(data = df.dyn20,id.vars = c("id","gene_short_name","highlight"))
# colnames(df.dyn30) <- c("id","gene_short_name","highlight","sample.name","abundance")
# df.dyn30$highlight.color <- ifelse(df.dyn30$highlight=="yes",df.dyn30$gene_short_name,"no")
# df.dyn40 <- subset(df.dyn30,!sample.name %in% c("IL10"))
# df.dyn40$day <- as.numeric(substr(x = df.dyn40$sample.name,start = 2,stop = 3))
# l <- ggplot(data = df.dyn40,mapping = aes(x=day,y=abundance,group=gene_short_name)) +
#   geom_line(mapping = aes(size=highlight,colour=highlight.color,alpha=highlight)) +
#   scale_size_manual(values=c(0.3,2)) +
#   scale_alpha_manual(values=c(0.2,1)) +
#   scale_color_pander() +
#   ylab("Relative abundance") +
#   xlab("Exp d") +
#   theme_classic(base_size = 24)
# print(l)

## compares Il10-/- levels to wildtype levels in a w-plot
# df.il <- data.frame(id=row.names(m.dyn),wt=m.dyn[,1],ko=m.dyn[,8])
# df.il10 <- melt(data = df.il,id.vars = c("id"))
# colnames(df.il10) <- c("id","sample.name","abundance")
# df.il10$abundance.log10 <- log10(df.il10$abundance)
# df.il20 <- dcast(data = df.il10,formula = id~sample.name,value.var = "abundance.log10")
# n <- ggplot(data = df.il10,mapping = aes(x=sample.name,y = abundance.log10)) +
#   geom_point() +
#   geom_segment(data=df.il20, mapping=aes(x=1,xend=2,y=wt,yend=ko),alpha=0.5) +
#   xlab("Sample") +
#   ylab("Mean abundance (log10)") +
#   theme_classic(base_size = 24)
# print(n)
# tt <- t.test(x = df.il20$wt,y = df.il20$ko)$p.value

## calculates relative abundance of reprogrammed stem cells over time
# barcodes <- colData(cds.stem) %>% data.frame
# m.stemcount <- table(barcodes$sample.name,barcodes$Ly6a.status)
# df.stemcount <- data.frame(sample.name=row.names(m.stemcount),lo=m.stemcount[,1],hi=m.stemcount[,2])
# df.stemcount10 <- melt(data = df.stemcount,id.vars = c("sample.name"))
# colnames(df.stemcount10) <- c("sample.name","Ly6a.level","num.cells")
# m.epicount <- table(colData(cds.epi.fil.fil)$sample)
# df.epicount <- data.frame(m.epicount)
# colnames(df.epicount) <- c("sample.num","epi")
# df.counts <- cbind(df.stemcount,y = df.epicount)
# df.counts$lo.frac <- df.counts$lo / df.counts$y.epi
# df.counts$hi.frac <- df.counts$hi / df.counts$y.epi
# df.counts$total.frac <- (df.counts$hi + df.counts$lo) / df.counts$y.epi
# df.counts10 <- melt(data = df.counts,id.vars = c("sample.name"),measure.vars = c("lo.frac","hi.frac","total.frac"))
# df.counts20 <- subset(df.counts10,substr(df.counts10$sample.name,1,1) == "D")
# df.counts20$day <- substr(df.counts20$sample.name,2,3) %>% as.numeric
# o <- ggplot(data = df.stemcount10,mapping = aes(x=sample.name,y=num.cells)) +
#   geom_bar(stat="identity",mapping=aes(fill=Ly6a.level),position="stack") +
#   ylab("# progenitor cells") +
#   xlab("Sample") +
#   scale_fill_pander() +
#   theme_classic(base_size = 24) +
#   theme(axis.text.x = element_text(angle=90,vjust = 0.5))
# print(o)
# o10 <- ggplot(data = df.counts20,mapping = aes(x = day,y = value)) +
#   geom_line(mapping = aes(group = variable,colour = variable),size=3) +
#   scale_color_colorblind() +
#   theme_classic(base_size = 24)
# print(o10)