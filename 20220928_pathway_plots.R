setwd("/Volumes/Macintosh HD/Users/cambrian/Dropbox/scrna-seq/")

library(enrichR)
library(reshape2)
library(pheatmap)
library(dplyr)
source(file = "camMonocleHelper.R")

dbopt <- listEnrichrDbs()

## looks at epithelium over time

# d <- read.csv(file = "table_epithelium_time_gene_modules_FPC_hybrid.csv",header = T)
# dbsoi <- c("MSigDB_Hallmark_2020")
# res <- query.enrichr(genes = d$gene_short_name,grouping = d$type,dbs = dbsoi,limit.result = 5,q.value = 0.05)

## makes time plot of epithelium pathways
# df.healing <- read.csv(file = "table_healing.csv",header = T)
# df.goi <- merge(x = df.healing, y = d,by= "gene_short_name",all = F)
# df.goi10 <- melt(data = df.goi,id.vars = c("type","id.x"),
#                  measure.vars = c("none","very.late","late","middle","early"),
#                  variable.name = "healing",value.name = "expression")
# summ <- df.goi10 %>% group_by(healing,type) %>% summarize(expr = median(expression)) %>% data.frame
# d10 <- merge(x = summ,y = res,by = "type",all = F)
# d10$heatscore <- d10$expr * -log10(d10$P.value)
# d20 <- dcast(data = d10,formula = Term~healing,value.var = "heatscore",fun.aggregate = median)
# m <- data.matrix(d20[,-1])
# row.names(m) <- d20$Term
# m10 <- m[,c(1,5,4,3,2)]
# m20 <- sweep(x = m10,MARGIN = 1,STATS = rowMins(m10),FUN = "-")
# m30 <- sweep(x = m20,MARGIN = 1,STATS = rowMaxs(m20),FUN = "/")
# ph <- pheatmap(mat = m30,cluster_cols = F)

## looks at differential display of FPC vs homeostatic stem cell genes
# da <- read.csv(file = "table_Ly6ahi_v_lo_specificity2_hybrid.csv",header = T)
# da10 <- subset(da,da$specificity > 1)
# da20 <- subset(da10,da10$pop.name %in% c("stem.hi","stem.lo"))
# dbsoi <- c("MSigDB_Hallmark_2020","ChEA_2022","KEGG_2021_Human","GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021","Reactome_2022")
# dbsoi <- c("ChEA_2022","KEGG_2021_Human","Reactome_2022")
# resa <- query.enrichr(genes = da20$gene_short_name,grouping = da20$pop.name,dbs = dbsoi,limit.result = 5,q.value = 0.01)

## looks at FPC specific genes
# d10 <- d[!is.na(d$pop.name),]
# dbsoi <- c("MSigDB_Hallmark_2020","WikiPathway_2021_Human")
# enr <- query.enrichr(genes = d10$gene_short_name,grouping = d10$pop.name,dbs = dbsoi,limit.result = 10,q.value = 0.15)

## looks at mesenchyme 1 genes
# d.mesen <- read.csv(file = "table_mesenchyme_1_time_gene_modules.csv",header = T)
# dbsoi <- c("MSigDB_Hallmark_2020")
# query.enrichr(genes = d.mesen$gene_short_name,grouping = d.mesen$type,dbs = dbsoi,limit.result = 8,q.value = 0.05)

## looks at mesenchyme 2 genes
# d.mesen <- read.csv(file = "table_mesenchyme_2_time_gene_modules.csv",header = T)
# dbsoi <- c("MSigDB_Hallmark_2020")
# query.enrichr(genes = d.mesen$gene_short_name,grouping = d.mesen$type,dbs = dbsoi,limit.result = 8,q.value = 0.05)