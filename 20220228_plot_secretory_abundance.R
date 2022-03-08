setwd("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/")

library(monocle3)
library(ggthemes)
source("camMonocleHelper.R")

# load(file = "cds_epi_curated.rds",verbose = T)
# barcodes <- colData(cds.epi.anno) %>% data.frame
# 
# goi <- c("Spink4","Neurod1","Dclk1","Clca1","Aqp8","Muc2","Ly6a","Krt20")
# 
# # cds.bleh <- cds.epi.anno[,subset(barcodes,!barcodes$sample.name == "IL10") %>% row.names]
# cds.bleh <- cds.epi.anno
# 
# cds.subset <- cds.bleh[map.gene.short.names(cds = cds.bleh,gene.short.names = goi),]
# colData(cds.subset)$Ly6a.status <- factor(colData(cds.subset)$Ly6a.status,levels = c("lo","hi"))
# 
# g <- plot_genes_violin(cds_subset = cds.subset,group_cells_by = "Ly6a.status",ncol = 2) +
#   scale_fill_excel_new() +
#   theme_clean(base_size = 16)
# print(g)

# ggsave(filename = "plot_secretory_Ly6a_expression.svg",plot = g)

## does significance testing
# cds.bleh <- add.gene.sum.column(cds = cds.bleh,genes.of.interest = c("Spink4"),col.name = "goi",use.raw = T)
# d <- data.frame(Ly6a.status = colData(cds.bleh)$Ly6a.status,goi = colData(cds.bleh)$goi)
# d10 <- subset(d,d$goi > -1)
# d10$log.goi <- log2(d10$goi)
# sig <- wilcox.test(formula = log.goi ~ Ly6a.status,data = d10)
# print(sig)
# summ <- d10 %>% group_by(Ly6a.status) %>% summarize(med.UMI = median(log.goi))
# print(summ)

## does significance testing the monocle3 way
# gene.fits <- fit_models(cds = cds.bleh,model_formula_str = "~Ly6a.status",verbose = T)
# fit.coefs <- coefficient_table(gene.fits)
# fit.coefs10 <- fit.coefs %>% filter(term == "Ly6a.statuslo")
# fit.coefs20 <- fit.coefs10 %>% filter(status == "OK")
# fit.coefs30 <- fit.coefs20[,c(1:4,7:15)] %>% data.frame

# write.table(x = fit.coefs30,file = "table_Ly6ahi_v_low_regress.csv",sep = ",",row.names = F)

# print(plot.expression.with.dots(cds.subset = cds.subset,grouping = "Ly6a.status"))