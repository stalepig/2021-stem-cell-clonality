setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)
source("camMonocleHelper.R")

load(file = "cds_all_aligned_fil.rds")

cds.il10 <- subset.cds.cells(cds = cds.allnn.fil2,column = "injury",elements = c("None","IL10"))
cds.il10.epi <- subset.cds.cells(cds = cds.il10,column = "annotation",elements = c("Epithelium"))

cds.il10.epi <- recluster.cds(cds = cds.il10.epi)
mtr <- top_markers(cds = cds.il10.epi)
colData(cds.il10.epi)$cell.group <- clusters(cds.il10.epi)
cds.il10.epi.stem <- subset.cds.cells(cds = cds.il10.epi,column = "cell.group",elements = c(2))

goi <- read.csv(file = "table_FPC_signature_conservative.csv",header = T)
cds.sub <- subset.cds.features(cds = cds.il10.epi.stem,gene.short.names = goi$gene_short_name)

metrics <- table(colData(cds.il10.epi.stem)$injury)
g <- plot_genes_by_group(cds = cds.sub,markers = goi$gene_short_name,group_cells_by = "injury")
print(g)

agg <- aggregate.expression.by.factor(cds = cds.sub,grouping = "injury",do.print = T,use.norm = T) * 100
aggdf <- data.frame(agg)
aggdf$id <- row.names(aggdf)
aggdf10 <- melt(data = aggdf,id.vars = "id",variable.name = "injury",value.name = "expr")
aggdf10$logexpr <- log2(aggdf10$expr)
aggdf10$injury <- factor(aggdf10$injury,c("None","IL10"))

sigvec <- agg[,"IL10"] - agg[,"None"]
t.test(x = sigvec) %>% print

# g <- ggplot(data = aggdf10,mapping = aes(x = injury,y = logexpr)) +
#   geom_point(size = 3) +
#   geom_line(mapping = aes(group=id)) +
#   ylab("Expression of FPC gene in progenitor cells (log)") +
#   theme_classic(base_size = 24) +
#   coord_fixed(ratio = 0.3)
# print(g)
# ggsave(filename = "plot_il10_comparison.svg",plot = g)