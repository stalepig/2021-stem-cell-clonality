setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)
library(ggplot2)
source("camMonocleHelper.R")

# load(file = "cds_all_aligned_fil.rds",verbose = T)

# feats <- rowData(cds.allnn.fil2) %>% data.frame
# goi <- feats[grep(pattern = "^Ifn",x = feats$gene_short_name),"gene_short_name"]

# g <- plot_genes_by_group(cds = cds.allnn.fil2,markers = goi,group_cells_by = "annotation.full",ordering_type = "maximal_on_diag")
# print(g)

# colData(cds.allnn.fil2)$time.cluster <- paste(colData(cds.allnn.fil2)$annotation,colData(cds.allnn.fil2)$healing,sep = "-")
# h <- plot_genes_by_group(cds = cds.allnn.fil2,markers = goi,group_cells_by = "time.cluster",ordering_type = "maximal_on_diag")
# print(h)

## makes the dotted plot
# cds.sub <- subset.cds.cells(cds = cds.allnn.fil2,column = "injury",elements = c("None","DSS"))
# g10 <- plot_genes_by_group(cds = cds.sub,markers = goi,group_cells_by = "annotation.full",ordering_type = "maximal_on_diag")
# print(g10)
# ggsave(filename = "plot_dots_Ifn_expression.svg",plot = g10)

## makes the UMAP
# cds.allnn.fil2 <- add.gene.sum.column(cds = cds.allnn.fil2,genes.of.interest = c("Ifng"),col.name = "idx.Ifng")
# barcodes <- colData(cds.allnn.fil2) %>% data.frame
# umap.df <- reducedDims(cds.allnn.fil2)[["UMAP"]] %>% data.frame
# barcodes10 <- cbind(barcodes,umap.df)
# g <- ggplot(data = barcodes10,mapping = aes(x = X1,y = X2)) +
#   geom_point(mapping = aes(size = idx.Ifng,alpha = idx.Ifng,colour = idx.Ifng)) +
#   scale_size_continuous(range = c(0.1,6)) +
#   scale_alpha_continuous(range = c(0.1,1)) +
#   ylab("UMAP-2") +
#   xlab("UMAP-1") +
#   labs(colour = "Ifng\nexpression") +
#   scale_colour_gradient2(low = "black",mid = "red",high = "blue",midpoint = 0.5) +
#   theme_classic(base_size = 24)
# print(g)
# ggsave(filename = "plot_umap_Ifng.png",plot = g)

## plots epithelial cluster expression of Ifngr1 over time of DSS healing
nonzero.mean <- function(vec) {
  vecfil <- subset(vec,vec > 0)
  if (length(vecfil) > 0) {
    return(mean(vecfil,na.rm=T))
  }
  else {
    return(0)
  }
}
fraction.expressing <- function(vec) {
  vecfil <- subset(vec,vec>0)
  return(length(vecfil) / length(vec))
}
# load(file = "cds_epi.rds")
# cds.epi.fil2 <- add.gene.sum.column(cds = cds.epi.fil2,genes.of.interest = c("Ifngr1"),col.name = "idx.Ifngr1",use.raw = T)
# barcodes<- colData(cds.epi.fil2) %>% data.frame
# summ <- barcodes %>% group_by(annotation,healing) %>% summarize(expr.ifngr1 = nonzero.mean(idx.Ifngr1),frac.ifngr1 = fraction.expressing(idx.Ifngr1),dumb.avg = mean(idx.Ifngr1))
k <- ggplot(data = summ,mapping = aes(x = healing,y = annotation)) +
  geom_point(mapping = aes(size = frac.ifngr1,colour = expr.ifngr1))
print(k)