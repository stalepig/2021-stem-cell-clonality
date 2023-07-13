setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/notacliu_shared/data/2022/07_01_22_Regev_UC_reanalysis/")

library(Matrix)
library(monocle3)

# exprs <- readMM(file = "gene_sorted-Epi.matrix.mtx")

# barcodes <- read.table(file = "Epi.barcodes2.tsv",sep = "\t",header = F)
# colnames(barcodes) <- c("barcode")
# row.names(barcodes) <- barcodes$barcode

# feats <- read.table(file = "Epi.genes.tsv",sep = "\t",header = F)
# colnames(feats) <- c("gene_short_name")
# row.names(feats) <- feats$gene_short_name

# colnames(exprs) <- barcodes$barcode
# row.names(exprs) <- feats$gene_short_name

## makes monocle3 cds object
# cds.full <- new_cell_data_set(expression_data = exprs,cell_metadata = barcodes,gene_metadata = feats)
# cds.full <- preprocess_cds(cds = cds.full,num_dim = 100,verbose = T)
# cds.full <- reduce_dimension(cds = cds.full,verbose = T)
# cds.full <- cluster_cells(cds = cds.full,verbose = T)
# cds.full <- learn_graph(cds = cds.full,verbose = T)

# meta <- read.table(file = "all.meta2.txt",header = T,sep = "\t")
# colData(cds.full) <- merge(x = colData(cds.full),y = meta,by.x = "barcode",by.y = "NAME",all = F)
# row.names(colData(cds.full)) <- colData(cds.full)$barcode