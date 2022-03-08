setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggplot2)
library(magrittr)
library(reshape2)
library(ggthemes)
library(dplyr)
library(scales)
library(umap)
library(readxl)
source("camMonocleHelper.R")

# load(file = "cds_epi_curated.rds",verbose = T)

# agg <- aggregate.expression.by.factor(cds = cds.epi.anno,grouping = "Ly6a.status",do.print = T)
# aggdf <- data.frame(agg)
# aggdf$id <- row.names(aggdf)
# feats <- rowData(cds.epi.anno) %>% data.frame
# d <- merge(x = feats,y = aggdf,by = "id")
# d$total <- d$hi + d$lo
# d$log2total <- log2(d$total)
# d$log2ratio <- log2(d$hi/d$lo)

## filters out low-expressing genes
# d10 <- subset(d, d$total > 0.02)

## assigns "significance"
# d10$mark <- sapply(X = d10$log2ratio,FUN = function(x) {
#   if (x > 3) { "up" }
#   else if (x < -2) { "down" }
#   else { "none" }
# })

## plots relationship between log expression and log2ratio
# g <- ggplot(data = d10,mapping = aes(x = log2total, y = log2ratio)) +
#   geom_point(mapping = aes(colour = mark)) +
#   ylab("Ly6a-hi vs. Ly6a-lo, log2") +
#   xlab("Total expression") +
#   theme_classic(base_size = 36) +
#   scale_color_colorblind()
# print(g)

## outputs data for pathway analysis
# d10.up <- subset(d10,d10$mark == "up")
# d10.down <- subset(d10,d10$mark == "down")
# write.table(x = d10.up,file = "table_Ly6ahi_up.csv",sep = ",",row.names = F)
# write.table(x = d10.down,file = "table_Ly6ahi_down.csv",sep = ",",row.names = F)

## plots UMAP split by Ly6a expression
# g20 <- plot_cells(cds = cds.epi.anno,color_cells_by = "Ly6a.status",label_cell_groups = F,show_trajectory_graph = F,cell_size = 1) +
#   scale_color_colorblind()
# print(g20)

## Processes the Ly6a-hi vs. Ly6a-lo genes and makes a mouse version of the gene set (all for Ly6a!)
# d20 <- d10[,c("gene_short_name","total","hi","lo")]
# colnames(d20) <- c("Name","Description","Ly6ahi","Ly6alo")
# write.table(x = d20,file = "gsea_Ly6a_all_mouse.txt",sep = "\t",row.names = F,quote = F)

## Makes a gene set for human GSEA
# hgene <- read.csv(file = "homologene.csv",header = F)
# colnames(hgene) <- c("ID","Organism","V3","Gene_symbol","V5","V6")
# hgene$Gene_symbol <- as.character(hgene$Gene_symbol)
# hgene10 <- dcast(data = hgene,formula = ID~Organism,value.var = "Gene_symbol",fun.aggregate = function(x) {x[1]})
# hgene20 <- hgene10[,c("ID","9606","10090")]
# colnames(hgene20) <- c("ID","Human","Mouse")
# d30 <- merge(x = hgene20,y = d20,by.x = "Mouse",by.y = "Name",all.x = F,all.y = T) %>%
#   na.omit
# d40 <- d30[,-c(1)]
# d50 <- d40[,c(2,1,4:5)]
# colnames(d50) <- c("Name","Description","Ly6ahi","Ly6alo")
# write.table(x = d50,file = "gsea_Ly6a_all_human.txt",sep = "\t",row.names = F,quote = F)

## Loads enrichr results
# da <- read_xlsx(path = "Ly6a_Pathway_Scores.xlsx",col_names = T) %>% data.frame
# da$logqval <- log10(da$qval)
# da$Score <- ifelse(da$Direction == "up",-1 * da$logqval,da$logqval)
# da$Annotation <- factor(da$Annotation,levels = as.character(da[order(da$Score,decreasing = F),"Annotation"]))
# ga <- ggplot(data = da,mapping = aes(y = Annotation,x = Score,fill = Category)) +
#   geom_col() +
#   scale_fill_wsj() +
#   theme_classic(base_size = 16)
# print(ga)

## Loads GSEA results
# db <- read_xlsx(path = "gsea/AllSymbols.Gsea.1616602432291/gsea_report_for_Ly6ahi_1616602432291.xlsx",col_names = T)
# db10 <- db[1:50,]
# gb <- ggplot(data = db10,mapping = aes(x = NES,y = NAME)) +
#   geom_col() +
#   theme_classic(base_size = 12)
print(gb)
