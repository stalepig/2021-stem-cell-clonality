setwd("~/Dropbox/scrna-seq/")

library(monocle3)
library(dplyr)
library(ggplot2)
library(reshape2)

# load(file = "cds_epi.rds")
# feats <- rowData(cds.epi.fil2) %>% data.frame
# refmarkers <- read.csv(file = "table_Ly6ahi_v_lo_specificity2_hybrid.csv",header = T)
# fpc.lib <- subset(refmarkers,refmarkers$pop.name=="stem.hi" & refmarkers$specificity>1 & refmarkers$gene_short_name != "conf")

# d.diff1 <- read.csv(file = "table_Ly6a_hi_v_lo_diff1_hybrid.csv",header = T)
# d.diff1$id <- row.names(d.diff1)
# d.diff1.clean <- merge(x = feats,y = d.diff1,by = "id",all = F)
# d.diff1.clean10 <- subset(d.diff1.clean,d.diff1.clean$gene_short_name %in% fpc.lib$gene_short_name)
# d.diff1.clean10$delta <- with(d.diff1.clean10,hi-lo)
print(t.test(d.diff1.clean10$delta-0.000,alternative = "g"))

# d.diff2 <- read.csv(file = "table_Ly6a_hi_v_lo_diff2_hybrid.csv",header = T)
# d.diff2$id <- row.names(d.diff2)
# d.diff2.clean <- merge(x = feats,y = d.diff2,by = "id",all = F)
# d.diff2.clean10 <- subset(d.diff2.clean,d.diff2.clean$gene_short_name %in% fpc.lib$gene_short_name)
# d.diff2.clean10$delta <- with(d.diff2.clean10,hi-lo)
print(t.test(d.diff2.clean10$delta-0.000,alternative = "g"))
# 
# d.goblet <- read.csv(file = "table_Ly6a_hi_v_lo_gob_hybrid.csv",header = T)
# d.goblet$id <- row.names(d.goblet)
# d.goblet.clean <- merge(x = feats,y = d.goblet,by = "id",all = F)
# d.goblet.clean10 <- subset(d.goblet.clean,d.goblet.clean$gene_short_name %in% fpc.lib$gene_short_name)
# d.goblet.clean10$delta <- with(d.goblet.clean10,hi-lo)
print(t.test(d.goblet.clean10$delta-0.000,alternative = "g"))

# d.prog <- read.csv(file = "table_Ly6a_hi_v_lo_prog_hybrid.csv",header = T)
# d.prog$id <- row.names(d.prog)
# d.prog.clean <- merge(x = feats,y = d.prog,by = "id",all = F)
# d.prog.clean10 <- subset(d.prog.clean,d.prog.clean$gene_short_name %in% fpc.lib$gene_short_name)
# d.prog.clean10$delta <- with(d.prog.clean10,hi-lo)
print(t.test(d.prog.clean10$delta-0.000,alternative = "g"))

# d <- data.frame(gene_short_name = d.prog.clean10$gene_short_name,
#                 prog = d.prog.clean10$delta,
#                 diff1 = d.diff1.clean10$delta,
#                 diff2 = d.diff2.clean10$delta,
#                 goblet = d.goblet.clean10$delta)
# d10 <- melt(data = d,id.vars = "gene_short_name",variable.name = "cell.group",value.name = "delta")

# g <- ggplot(data = d10,mapping = aes(x = cell.group,y = delta)) +
#   geom_violin(size=1) +
#   geom_segment(mapping = aes(x = 0,xend = 5,y = 0,yend = 0),colour = "red") +
#   theme_classic(base_size = 24) +
#   coord_fixed(ratio = 50)
# print(g)
# ggsave(filename = "plot_FPC_expression_diff_cells.svg",plot = g)