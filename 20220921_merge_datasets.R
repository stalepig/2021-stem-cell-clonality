setwd("/Volumes/MacintoshHD/Users/cambrian/Dropbox/scrna-seq/")

library(monocle3)

# run1.d00 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D0/")
# run1.d06 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D3/")
# run1.d09 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D9/")
# run1.d12 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D12/")
# run1.d16 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D16/")
# run1.d19 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D19/")
# run1.d23 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/D23/")
# run1.il10 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/novaseq/out/IL10/")
# 
# run2.d00 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2020/02_14_20_metaplasia_single_cell/D00_out/")
# run2.d07 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2020/02_14_20_metaplasia_single_cell/D07_out/")
# run2.d42 <- load_cellranger_data(pipestance_path = "/Volumes/SHAQ/2020/02_14_20_metaplasia_single_cell/D42_out/")
# 
# run3.d07 <- load_cellranger_data(pipestance_path = "/Volumes/LEBRON/2022/08_08_22_DSS_scRNASeq_Chicago/D07/")
# run3.d11 <- load_cellranger_data(pipestance_path = "/Volumes/LEBRON/2022/08_08_22_DSS_scRNASeq_Chicago/D11/")
# run3.d18 <- load_cellranger_data(pipestance_path = "/Volumes/LEBRON/2022/08_08_22_DSS_scRNASeq_Chicago/D18/")

# cds.all <- combine_cds(cds_list = list(run1.d00,run1.d06,run1.d09,run1.d12,run1.d16,run1.d19,run1.d23,run1.il10,
#                                        run2.d00,run2.d07,run2.d42,
#                                        run3.d07,run3.d11,run3.d18),keep_all_genes = F)

# save(cds.all,file = "cds_all.rds")