setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/scrna-seq/")

library(enrichR)
library(reshape2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggthemes)
source(file = "camMonocleHelper.R")

dbopt <- listEnrichrDbs()

dbsoi <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021","KEGG_2019_Mouse","Mouse_Gene_Atlas",
           "MSigDB_Hallmark_2020","MSigDB_Oncogenic_Signatures","Tabula_Muris","ChEA_2022","WikiPathway_2021_Human","Reactome_2022")

# d <- read.csv(file = "table_FPC_signature_conservative.csv",header = T)

# res <- query.enrichr(genes = d$gene_short_name,grouping = d$pop.name,dbs = dbsoi,limit.result = 3,q.value = 0.01)

# res$Term <- factor(res$Term,levels=res[order(-log10(res$Adjusted.P.value),decreasing = F),]$Term)
res10 <- res[order(-log10(res$Adjusted.P.value),decreasing = T),]
res20 <- res10[1:10,]
g <- ggplot(data = res20,mapping = aes(x = -log10(Adjusted.P.value),y = Term)) +
  geom_col(mapping = aes(fill = database),alpha=1) +
  scale_fill_colorblind() +
  theme_classic(base_size = 12)
print(g)
ggsave(filename = "plot_bar_FPC_specific_pathways.svg",plot = g)