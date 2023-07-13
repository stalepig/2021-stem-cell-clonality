setwd("/Volumes/Macintosh HD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)
library(readxl)
library(magrittr)
library(tibble)

# load(file = "cds_all.rds",verbose = T)
# barcodes <- colData(cds.all) %>% data.frame
# meta <- read_xlsx(path = "metadata_batch.xlsx",col_names = T)
# colData(cds.all) <- merge(x = colData(cds.all),y = meta,by = "sample",all = T,)
# row.names(colData(cds.all)) <- paste(colData(cds.all)$barcode,colData(cds.all)$sample,sep = "_")

# cds.allnn <- align_cds(cds = cds.all,alignment_group = "batch",verbose = T)
# cds.allnn <- reduce_dimension(cds = cds.allnn)
# cds.allnn <- cluster_cells(cds = cds.allnn,verbose = T)

# save(cds.allnn,file = "cds_all_aligned.rds")