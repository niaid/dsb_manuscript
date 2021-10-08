# R4 Seurat 4 (see session info) 
suppressMessages(library(dsb))
suppressMessages(library(ggridges))
suppressMessages(library(magrittr))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# project info 
project_title = "joint_model"
figpath = paste0(here("V2/joint_clustering/"), project_title, "/figures/"); dir.create(figpath, recursive = TRUE)
datapath = paste0(here("V2/joint_clustering/"), project_title, "/generated_data/"); dir.create(datapath, recursive = TRUE)


# Read seurat version 2 object and extract data 
s = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds'))
rna_raw = s@raw.data
adt_raw = s@assay$CITE@raw.data
adt_dsb = s@assay$CITE@data
md = s@meta.data
rm(s); gc()

# joint lcustering pipeline  
# create v4 object starting with rna 
h1 = CreateSeuratObject(counts = rna_raw,meta.data = md)

# normalize RNA and sc transcorm 
DefaultAssay(h1) <- 'RNA'
h1 = NormalizeData(h1) %>% FindVariableFeatures(selection.method = 'vst', verbose = TRUE) %>% ScaleData() %>% RunPCA()

# define non staining proteins and define raw protein data 
non_staining = c("CX3CR1-PROT", "CD357-PROT", "CD275-PROT", "CD152-PROT", "CD294-PROT",
                 "CD70-PROT", "CD90-PROT", "Mouse IgG2bkIsotype-PROT", "MouseIgG1kappaisotype-PROT",  
                 "RatIgG2bkIsotype-PROT", "CD206-PROT", "MouseIgG2akappaisotype-PROT", "CD223-PROT", 
                 "CD138-PROT", "CD274-PROT", "CD137-PROT", "CD273-PROT","CD80-PROT")

adt_raw = as.matrix(adt_raw[setdiff(rownames(adt_raw), non_staining), ])
rownames(adt_raw) = str_replace_all(rownames(adt_raw),pattern = "_", replacement = "-")

# set clr assay 
adt_assay = CreateAssayObject(counts = adt_raw)
h1[["CLR"]] <- adt_assay

# run clr Joint clustering workflow 
DefaultAssay(h1) <- 'CLR'
h1 = NormalizeData(h1, normalization.method = 'CLR', margin = 2) %>% ScaleData() 

# hack seurat to use normalized protein values as a dimensionality reduction object.
VariableFeatures(h1) <- rownames(adt_raw)
# true pca
h1 = RunPCA(h1, reduction.name = 'pCLR', features = VariableFeatures(h1))
# make matrix of norm values to add as dr embeddings
pseudo = t(h1@assays$CLR@data)
p_colnames = paste('pseudo', 1:69, sep = "_")
colnames(pseudo) = p_colnames
# add to object 
h1@reductions$pCLR@cell.embeddings = pseudo

# run WNN 
h1 = FindMultiModalNeighbors(
  object = h1, 
  reduction.list = list("pca", "pCLR"),
  weighted.nn.name = "clr_pseudo_wnn", 
  knn.graph.name = "clr_pseudo_knn",
  modality.weight.name = "clr_pseudo_weight",
  snn.graph.name = "clr_pseudo_snn",
  dims.list = list(1:30, 1:69)
)
h1 = FindClusters(h1, graph.name = "clr_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3, resolution = 2, random.seed = 1990,  verbose = FALSE)
h1 = FindClusters(h1, graph.name = "clr_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3, resolution = 3, random.seed = 1990,  verbose = FALSE)
h1 = RunUMAP(h1, nn.name = "clr_pseudo_wnn", reduction.name = "CLR_wnn_umap", reduction.key = "CLR_wnnUMAP_", seed.use = 1990)

####################################
# dsb workflow 
# adt_assay = CreateAssayObject(counts = adt_raw)
h1[["dsb"]] <- adt_assay

# run dsb workflow without pca 
DefaultAssay(h1) <- 'dsb'
# add dsb normalized matrix to data slot 
rownames(adt_dsb) = str_replace_all(rownames(adt_dsb), pattern = "_", replacement = "-")
h1@assays$dsb@data = as.matrix(adt_dsb[setdiff(rownames(adt_dsb), non_staining), ])
h1 = ScaleData(h1, assay = 'dsb')

# hack seurat to use normalized protein values as a dimensionality reduction object.
VariableFeatures(h1) <- rownames(adt_raw)
# true pca
h1 = RunPCA(h1, reduction.name = 'pdsb', features = VariableFeatures(h1))
# make matrix of norm values to add as dr embeddings
pseudo = t(h1@assays$dsb@data)
fake_colnames = paste('pseudo', 1:69, sep = "_")
colnames(pseudo) = fake_colnames
# add to object 
h1@reductions$pdsb@cell.embeddings = pseudo

# run WNN 
h1 = FindMultiModalNeighbors(
  object = h1,
  reduction.list = list("pca", "pdsb"),
  weighted.nn.name = "dsb_pseudo_wnn", 
  knn.graph.name = "dsb_pseudo_knn",
  modality.weight.name = "dsb_pseudo_weight",
  snn.graph.name = "dsb_pseudo_snn",
  dims.list = list(1:30, 1:69)
)

h1 = FindClusters(h1, graph.name = "dsb_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3, resolution = 2, random.seed = 1990,  verbose = FALSE)
h1 = FindClusters(h1, graph.name = "dsb_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3, resolution = 3, random.seed = 1990,  verbose = FALSE)
h1 = RunUMAP(h1, nn.name = "dsb_pseudo_wnn", reduction.name = "dsb_wnn_umap", reduction.key = "dsb_wnnUMAP_", seed.use = 1990)

# find dsb model RNA markers
# enable parallelization 
library(future)
plan("multiprocess", workers = 4)

# roc test for high res WNN clusters 
Idents(h1) = 'dsb_pseudo_knn_res.3'
DefaultAssay(h1) = 'RNA'
vf = VariableFeatures(s, assay = "RNA")
dsb_rna_m = FindAllMarkers(object = h1, assay = 'RNA', test.use = 'roc', features = vf)
write_delim(dsb_rna_m,file = paste0(datapath, 'dsb_rna_m.txt'))

# save object 
saveRDS(h1,file = paste0(datapath, 'h1_WNN.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
# [8] tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       cowplot_1.1.1      magrittr_2.0.1    
# [15] ggridges_0.5.3     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2        splines_4.0.5        
# [7] listenv_0.8.0         scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1     viridis_0.5.1         fansi_0.4.2          
# [13] tensor_1.5            cluster_2.1.1         ROCR_1.0-11           openxlsx_4.2.3        limma_3.46.0          globals_0.14.0       
# [19] modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0 colorspace_2.0-0      rvest_0.3.6           ggrepel_0.9.1        
# [25] haven_2.3.1           xfun_0.21             crayon_1.4.1          jsonlite_1.7.2        spatstat.data_2.1-0   survival_3.2-10      
# [31] zoo_1.8-8             glue_1.4.2            polyclip_1.10-0       gtable_0.3.0          leiden_0.3.7          car_3.0-10           
# [37] future.apply_1.7.0    abind_1.4-5           scales_1.1.1          pheatmap_1.0.12       DBI_1.1.1             rstatix_0.7.0        
# [43] miniUI_0.1.1.1        Rcpp_1.0.6            viridisLite_0.3.0     xtable_1.8-4          reticulate_1.18       spatstat.core_2.0-0  
# [49] foreign_0.8-81        bit_4.0.4             mclust_5.4.7          htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2   
# [55] ellipsis_0.3.1        ica_1.0-2             pkgconfig_2.0.3       farver_2.0.3          uwot_0.1.10           dbplyr_2.1.0         
# [61] deldir_0.2-10         utf8_1.1.4            tidyselect_1.1.0      labeling_0.4.2        rlang_0.4.10          reshape2_1.4.4       
# [67] later_1.1.0.1         munsell_0.5.0         cellranger_1.1.0      tools_4.0.5           cli_2.3.0             generics_0.1.0       
# [73] broom_0.7.5           evaluate_0.14         fastmap_1.1.0         yaml_2.2.1            goftest_1.2-2         knitr_1.31           
# [79] bit64_4.0.5           fs_1.5.0              fitdistrplus_1.1-3    zip_2.1.1             RANN_2.6.1            pbapply_1.4-3        
# [85] future_1.21.0         nlme_3.1-152          mime_0.10             xml2_1.3.2            hdf5r_1.3.3           compiler_4.0.5       
# [91] rstudioapi_0.13       curl_4.3              plotly_4.9.3          png_0.1-7             ggsignif_0.6.0        spatstat.utils_2.1-0 
# [97] reprex_1.0.0          stringi_1.5.3         RSpectra_0.16-0       lattice_0.20-41       Matrix_1.3-2          ggsci_2.9            
# [103] vctrs_0.3.6           pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18     
# [109] data.table_1.14.0     irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       R6_2.5.0              promises_1.2.0.1     
# [115] rio_0.5.16            KernSmooth_2.23-18    gridExtra_2.3         parallelly_1.23.0     codetools_0.2-18      MASS_7.3-53.1        
# [121] assertthat_0.2.1      rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           parallel_4.0.5       
# [127] hms_1.0.0             grid_4.0.5            rpart_4.1-15          rmarkdown_2.7         carData_3.0-4         Rtsne_0.15           
# [133] ggpubr_0.4.0          shiny_1.6.0           lubridate_1.7.9.2   