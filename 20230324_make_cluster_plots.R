setwd("/Users/cambrian/Dropbox/scrna-seq")

library(monocle3)
source("camMonocleHelper.R")
library(ggplot2)
library(ggthemes)

## does this for all cells
# load(file = "cds_all_aligned_fil.rds",verbose = T)
# tm <- top_markers(cds = cds.allnn.fil2,group_cells_by = "annotation.full",verbose = T)
# tm10 <- tm[-grep(pattern = "^mt-",x = tm$gene_short_name),]
# goi <- filter.marker.results(marker_test_res = tm10,num.to.keep = 3,criterion = "marker_score")
# g <- plot_genes_by_group(cds = cds.allnn.fil2,markers = goi,group_cells_by = "annotation.full",max.size = 4) +
#   coord_fixed(ratio = 0.5)
# print(g)
# ggsave(filename = "plot_dots_global_markers_pub.svg",plot = g)
# write.table(x = tm,file = "table_global_cluster_markers.csv",sep = ",",row.names = F)



## does this for epithelium only
# load(file = "cds_epi.rds")
# tm <- top_markers(cds = cds.epi.fil2,group_cells_by = "annotation",verbose = T)
# tm10 <- tm[!grepl(pattern = "^mt-",x = tm$gene_short_name),]
# goi <- c(filter.marker.results(marker_test_res = tm10,num.to.keep = 5,criterion = "marker_score"),map.gene.short.names(cds = cds.epi.fil2,gene.short.names = c("Mgat4c")))
# g <- plot_genes_by_group(cds = cds.epi.fil2,markers = goi,group_cells_by = "annotation",max.size = 4) +
#   coord_fixed(ratio = 0.5)
# print(g)
# ggsave(filename = "plot_dots_epithelial_markers_pub.svg",plot = g)
# write.table(x = tm,file = "table_epithelial_cluster_markers.csv",sep = ",",row.names = F)

## plots reprogramming for stromal cells
# load(file = "cds_all_aligned_fil.rds",verbose = T)
# barcodes <- cbind(data.frame(colData(cds.allnn.fil2)),reducedDims(cds.allnn.fil2)[["UMAP"]])
# colnames(barcodes)[(dim(barcodes)[[2]]-1):(dim(barcodes)[[2]])] <- c("UMAP.1","UMAP.2")
# barcodes$injury <- factor(barcodes$injury,levels=c("None","DSS","IL10"))

## subsamples for plotting
# tab <- table(barcodes$healing,barcodes$annotation)
# barcodes10 <- subset(barcodes,barcodes$annotation %in% c("Endothelium","Epithelium","Immune","Lymphatic","Mesenchyme"))
# barcodes20 <- split(barcodes10,list(barcodes10$annotation)) %>%
#   lapply(FUN = function(df) {
#     sample.size <- min(table(df$injury))
#     retdf <- split(df,list(df$injury)) %>%
#       lapply(FUN = function(indf) {
#         idx <- sample(x = dim(indf)[[1]],size = sample.size,replace = F)
#         return(indf[idx,])
#       }) %>% do.call(what = "rbind",args = .)
#     retdf$sample.size <- sample.size
#     return(retdf)
#   }) %>% do.call(what = "rbind",args = .)

g <- ggplot(data = barcodes20,mapping = aes(x = UMAP.1,y = UMAP.2)) +
  facet_wrap(facets = .~annotation,scales = "free") +
  geom_point(mapping = aes(colour = injury,alpha = 1/sample.size),size=0.25) +
  scale_colour_colorblind() +
  guides(alpha = F) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  scale_alpha_continuous(range = c(0.2,1.0)) +
  theme_classic(base_size = 22)
print(g)
ggsave(filename = "plot_umap_stromal_reprogramming.png",plot = g)