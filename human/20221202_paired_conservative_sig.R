setwd("/Volumes/Macintosh HD/Users/cambrian/Dropbox/notacliu_shared/data/2022/07_01_22_Regev_UC_reanalysis/")

library(monocle3)
library(dplyr)
library(readxl)
library(ggthemes)
source("/Volumes/Macintosh HD/Users/cambrian/Dropbox/scrna-seq/camMonocleHelper.R")

# load(file = "cds_filtered_paired.rds",verbose = T)
# goi <- toupper(x = read.csv(file = "table_FPC_signature_conservative.csv",header = T)$gene_short_name)
# goi.human <- subset(goi,goi %in% rowData(cds.paired)$gene_short_name)

## does simple aggregation
# cds.paired.sub <- cds.paired[goi.human,]
# colData(cds.paired.sub)$Condition <- paste(colData(cds.paired.sub)$Subject,colData(cds.paired.sub)$Health,sep = ";")
# agg <- aggregate.expression.by.factor(cds = cds.paired.sub,grouping = "Condition",do.print = T,use.norm = T)
# aggdf <- data.frame(t(agg))
# aggdf$Condition <- row.names(aggdf)
# aggdf10 <- melt(data = aggdf,id.vars = c("Condition"),variable.name = "Gene",value.name = "Expression")
# aggdf10$Subject <- sapply(X = aggdf10$Condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = ";")[[1]][1]
# })
# aggdf10$Health <- sapply(X = aggdf10$Condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = ";")[[1]][2]
# })

## plots simple aggregation
# h <- ggplot(data = aggdf10,mapping = aes(x = Health,y = Expression)) +
#   facet_wrap(facets = .~Gene,scales = "free_y") +
#   geom_point() +
#   geom_line(mapping = aes(group = Subject))
# print(h)

# h10 <- plot_genes_by_group(cds = cds.paired.sub,markers = goi.human,group_cells_by = "Cell.Type")
# print(h10)

## aggregates by human cell type to prove that these are, on the whole, progenitor cell markers
# agg <- aggregate.expression.by.factor(cds = cds.paired.sub,grouping = "Cell.Type",do.print = T,use.norm = T)
# agg10 <- scale(x = t(agg),center = F,scale = T)
# aggdf <- data.frame(agg10)
# aggdf$Cell.Type <- row.names(aggdf)
# aggdf10 <- melt(data = aggdf,id.vars = "Cell.Type")
# summ <- aggdf10 %>% group_by(Cell.Type) %>% summarize(Med.Expr = median(value))
# h20 <- ggplot(data = aggdf10,mapping = aes(x = Cell.Type,y = value)) +
#   geom_point() +
#   geom_col(mapping = aes(y = Med.Expr),data = summ,alpha=0.3) +
#   ylab("FPC signature\nexpression (norm)") +
#   xlab("Cell cluster") +
#   coord_fixed(ratio = 2) +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(h20)
# ggsave(filename = "plot_signature_expression_by_cluster.svg",plot = h20)

## now checks whether the expression is increased in inflamed conditions of each of the patients
# barcodes <- colData(cds.paired.sub) %>% data.frame
# cds.pp.sub <- cds.paired.sub[,subset(barcodes,barcodes$Cell.Type %in% c("Sec-Ent-DUOX2","Prog","Prog-Mitotic")) %>% row.names]
# colData(cds.pp.sub)$Condition <- with(colData(cds.pp.sub),paste(Subject,Health,Cell.Type,sep = ".."))
# agg <- aggregate.expression.by.factor(cds = cds.pp.sub,grouping = "Condition",do.print = T,use.norm = T)
# aggdf <- data.frame(agg)
# aggdf$gene_short_name <- row.names(aggdf)
# aggdf10 <- melt(data = aggdf,id.vars = "gene_short_name",variable.name = "condition",value.name = "expr")
# aggdf10$patient <- sapply(X = aggdf10$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = "\\.\\.")[[1]][1]
# })
# aggdf10$health <- sapply(X = aggdf10$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = "\\.\\.")[[1]][2]
# })
# aggdf10$cell.type <- sapply(X = aggdf10$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = "\\.\\.")[[1]][3]
# })
# spl <- split(x = aggdf10,f = list(aggdf10$gene_short_name,aggdf10$patient,aggdf10$cell.type),drop = T)
# comb <- lapply(X = spl,FUN = function(df) {
#   diff <- df[2,"expr"] - df[1,"expr"]
#   data.frame(gene_short_name = df[1,"gene_short_name"],
#              diffexp = diff,
#              cell.type = df[1,"cell.type"],
#              patient = df[1,"patient"])
# }) %>% do.call("rbind",.) %>% na.omit
# summ <- comb %>% group_by(patient,cell.type) %>% summarize(med.diffexp = mean(diffexp,na.rm = T))
# summa <- comb %>% group_by(patient) %>% summarize(idx = mean(diffexp,na.rm = T))
# comb10 <- merge(x = comb,y = summa,by = "patient",all = F)
# h30 <- ggplot(data = comb10,mapping = aes(x = patient,y = diffexp)) +
#   facet_grid(facets = cell.type~.,rows = 3) +
#   geom_point(mapping = aes(colour = idx)) +
#   geom_col(mapping = aes(y = med.diffexp),data = summ,alpha = 0.3) +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(h30)
# ggsave(filename = "plot_signature_progenitors_by_patient.png",plot = h30)

## performs statistical test for each patient
# summ.t <- comb %>% group_by(patient) %>% summarize(pval = t.test(diffexp)$p.value)
# summ.t$fdr <- p.adjust(summ.t$pval,method = "bonferroni")

## looks to see if the overall trend is significant
# h40 <- ggplot(data = summ,mapping = aes(x = cell.type,y = med.diffexp)) +
#   geom_boxplot(size=2) +
#   geom_segment(mapping = aes(x = 0,xend = 3.5,y = 0.0,yend = 0.0),colour = "red") +
#   theme_classic(base_size = 16) +
#   coord_fixed(ratio = 10)
# print(h40)
# ggsave(filename = "plot_signature_progenitors_all_patients.png",plot = h40)
# print(t.test(x = subset(summ,summ$cell.type == "Prog")$med.diffexp))
# print(t.test(x = subset(summ,summ$cell.type == "Prog.Mitotic")$med.diffexp))
# print(t.test(x = subset(summ,summ$cell.type == "Sec.Ent.DUOX2")$med.diffexp))

## makes a few umap plots
# h50 <- plot_cells(cds = cds.paired.sub,genes = goi.human,show_trajectory_graph = F)
# ggsave(filename = "plot_signature_UMAP.png",plot = h50)

## gets a marker plot
# mtr <- top_markers(cds = cds.paired,group_cells_by = "Cell.Type",verbose = T)
# tm <- filter.marker.results(marker_test_res = mtr,num.to.keep = 2,criterion = "marker_score")
# h60 <- plot_genes_by_group(cds = cds.paired,markers = tm,group_cells_by = "Cell.Type")
# ggsave(filename = "plot_marker_genes.svg",plot = h60)