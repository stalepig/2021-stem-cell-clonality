setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

source("camMonocleHelper.R")
library(randomForest)
library(MASS)
library(ROCR)
library(ggplot2)
library(ggthemes)

# load(file = "cds_epi_stem.rds")
# markers.df <- read.csv(file = "table_stem_repro_dynamics.csv",header = T)
# 
# cds.stem <- add.gene.sum.column(cds = cds.stem,genes.of.interest = c("Ly6a"),col.name = "idx.Ly6a",use.raw=T)
# colData(cds.stem)$Ly6a.status <- ifelse(log10(colData(cds.stem)$idx.Ly6a+0.1) > 0.783,"hi","lo")
# colData(cds.stem)$Ly6a.status <- factor(colData(cds.stem)$Ly6a.status,levels=c("lo","hi"))
# 
# cds.stem.subset <- cds.stem[as.character(markers.df$id),]
# expr.mat <- t(as.matrix(exprs(cds.stem.subset)))
# d <- data.frame(expr.mat)
# d$Ly6a.status <- colData(cds.stem.subset)$Ly6a.status
# 
# rf <- randomForest(formula = Ly6a.status ~ .,data = d)
# summ <- data.frame(importance(rf))
# summ$id <- row.names(summ)
# feats <- rowData(cds.stem.subset)
# summ10 <- merge(x = data.frame(feats),y = summ,by = "id",all = F)
# summ20 <- summ10[order(summ10$MeanDecreaseGini,decreasing = T),]
# 
# d10 <- d[,which(colnames(d) %in% summ20[1:50,"id"]),drop = F]
# d10$Ly6a.status <- d$Ly6a.status
# ld <- lda(formula = Ly6a.status ~ .,data = d10)
# ld.predict <- predict(object = ld,newdata = d10)
# 
# 
# summa <- data.frame(pred = ld.predict$class,actual = d10$Ly6a.status,cell_id = row.names(d10))
# summa$correct <- mapply(FUN = function(pred,act) {
#   if (pred == act) {
#     1
#   } else {
#     0
#   }
# },summa$pred,summa$actual)
# print(sum(summa$correct) / dim(summa)[[1]] * 100)

# summ30 <- summ20[1:20,]
# summ30$gene_short_name <- factor(summ30$gene_short_name,levels=rev(summ30$gene_short_name))
# g <- ggplot(data = summ30,mapping = aes(x = MeanDecreaseGini,y = gene_short_name)) +
#   geom_col() +
#   theme_classic(base_size = 24) +
#   coord_fixed(ratio = 2)
# print(g)
# ggsave(filename = "plot_ASSCA_random_forest.svg",plot = g)


## loads the global CDS to evaluate ASSCA expression in other clusters
# load(file = "cds_full.rds")
goi <- c("S100a11","Anxa2","Capg","S100a16","Aprt")
goi.ids <- map.gene.short.names(cds = cds.full,gene.short.names = goi)
cds.full.subset <- cds.full[goi.ids,]
h <- plot_genes_violin(cds_subset = cds.full.subset,group_cells_by = "annotation",log_scale = F,ncol = 2) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5))
print(h)
ggsave(filename = "plot_ASSCA_violin.svg",plot = h)

## in case you need ROC curve
# pred <- prediction(predictions = data.frame(ld.predict$posterior[,2]),labels = d10$Ly6a.status)
# roc.perf <- performance(prediction.obj = pred,measure = "spec",x.measure = "sens")
# plot(roc.perf)