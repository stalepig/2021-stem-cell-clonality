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

## computes Ly6a and Birc5 cell abundance and gates based on raw values
## then performs cluster specific feature identification on normalized values
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Ly6a"),col.name = "Expr.Ly6a",use.raw = T,use.norm = F)
# print(hist(colData(cds.epi.fil2)$Expr.Ly6a))
# threshold.Ly6a <- mean(c(summary(colData(cds.epi.fil2)$Expr.Ly6a)[3],summary(colData(cds.epi.fil2)$Expr.Ly6a)[5]))
# colData(cds.epi.fil2)$status.Ly6a <- sapply(X = colData(cds.epi.fil2)$Expr.Ly6a,FUN = function(x) {
#   if (x > threshold.Ly6a) "hi"
#   else "lo"
# })

## gets just stem cell cluster
# barcodes <- colData(cds.epi.fil2) %>% data.frame
# cds.epi.stem <- cds.epi.fil2[,subset(barcodes,barcodes$annotation == "Prog-Gpx2") %>% row.names]
# goi <- c("Sulf2","Clu","Ass1","Lgr5","Hopx","Gpx2")
# cds.stem.sub <- cds.epi.stem[map.gene.short.names(cds = cds.epi.stem,gene.short.names = goi),]
# colData(cds.stem.sub)$condition <- paste(colData(cds.stem.sub)$healing,colData(cds.stem.sub)$status.Ly6a,sep=";")
# agg <- aggregate.expression.by.factor(cds = cds.stem.sub,grouping = "condition",do.print = T,use.norm = F)
# aggdf <- data.frame(t(agg))
# aggdf$condition <- row.names(aggdf)
# aggdf10 <- melt(data = aggdf,id.vars = "condition",variable.name = "id",value.name = "expression")
# aggdf20 <- merge(x = data.frame(rowData(cds.stem.sub)),y = aggdf10,by = "id",all = F)
# aggdf20$healing <- sapply(X = aggdf20$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = ";")[[1]][1]
# })
# aggdf20$Ly6a <- sapply(X = aggdf20$condition,FUN = function(txt) {
#   strsplit(x = as.character(txt),split = ";")[[1]][2]
# })
# aggdf20$healing <- factor(aggdf20$healing,levels=c("none","early","middle","late","very-late"))
# g <- ggplot(data = aggdf20,mapping = aes(x = healing,y = expression)) +
#   facet_wrap(facets = .~gene_short_name,scales = "free_y") +
#   geom_col(mapping = aes(fill = Ly6a)) +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(g)
# ggsave(filename = "plot_bar_FPC_markers.svg",plot = g)
# print(table(colData(cds.stem.sub)$condition))

# goi <- c("Sulf2","Clu","Ass1","Ly6a","Lgr5","Hopx")
# colData(cds.epi.fil2)$healing <- factor(colData(cds.epi.fil2)$healing,levels = c("none","early","middle","late","very-late"))
# cds.sub <- cds.epi.fil2[map.gene.short.names(cds = cds.epi.fil2,gene.short.names = goi),]
# h <- plot.expression.with.dots(cds.subset = cds.sub,grouping = "healing",yjitter = 0.2,logticks = T,use.norm = F)
# h <- h + theme_classic(base_size = 20) +
#   theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
# print(h)
# ggsave(filename = "plot_dots_genes_time.svg",plot = h)

