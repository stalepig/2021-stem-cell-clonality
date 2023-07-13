setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/notacliu_shared/data/2022/07_01_22_Regev_UC_reanalysis/")

library(monocle3)
library(dplyr)
library(readxl)
library(ggthemes)
source("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/camMonocleHelper.R")

load(file = "cds_filtered.rds")

barcodes <- colData(cds.fil) %>% data.frame
umap <- reducedDims(cds.fil)$UMAP %>% data.frame
colnames(umap) <- c("umap.1","umap.2")
barcodes10 <- cbind(barcodes,umap)

## filters out the patients that don't have inflamed/non-inflamed pairs
tab <- table(barcodes$Subject,barcodes$Health)
to.keep <- apply(X = tab,MARGIN = 1,FUN = function(x) {
  if (x[2] > 0 & x[3] > 0) {
    T
  } else {
    F
  }
})
patients.to.keep <- row.names(tab)[to.keep]
barcodes20 <- subset(barcodes10,barcodes10$Subject %in% patients.to.keep)
cds.paired <- cds.fil[,subset(barcodes10,barcodes10$Subject %in% patients.to.keep) %>% row.names]

save(cds.paired,file = "cds_filtered_paired.rds")

## plots umap for patients to keep
# g <- ggplot(data = barcodes20,mapping = aes(x = umap.1,y = umap.2)) +
#   facet_wrap(facets = .~Subject) +
#   geom_point(mapping = aes(colour = Health),size=0.25,alpha=0.5)
# print(g)
# ggsave(filename = "plot_umap_by_patient.png",plot = g,dpi = 300)

## do some simple aggregation
goi <- c("ANXA2","S100A11","S100A16","CAPG","APRT","ASS1","GPX2","SULF2")
cds.paired.sub <- cds.paired[goi,]
colData(cds.paired.sub)$Condition <- paste(colData(cds.paired.sub)$Subject,colData(cds.paired.sub)$Health,sep = ";")
agg <- aggregate.expression.by.factor(cds = cds.paired.sub,grouping = "Condition",do.print = T)
aggdf <- data.frame(t(agg))
aggdf$Condition <- row.names(aggdf)
aggdf10 <- melt(data = aggdf,id.vars = c("Condition"),variable.name = "Gene",value.name = "Expression")
aggdf10$Subject <- sapply(X = aggdf10$Condition,FUN = function(txt) {
  strsplit(x = as.character(txt),split = ";")[[1]][1]
})
aggdf10$Health <- sapply(X = aggdf10$Condition,FUN = function(txt) {
  strsplit(x = as.character(txt),split = ";")[[1]][2]
})

## plots simple aggregation
h <- ggplot(data = aggdf10,mapping = aes(x = Health,y = Expression)) +
  facet_wrap(facets = .~Gene,scales = "free_y") +
  geom_point() +
  geom_line(mapping = aes(group = Subject))
print(h)
