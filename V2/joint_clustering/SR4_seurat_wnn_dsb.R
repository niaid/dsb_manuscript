# R4 Seurat4 
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

# joint clustering pipeline  
# create v4 object starting with rna 
h1 = CreateSeuratObject(counts = rna_raw,meta.data = md)

# normalize RNA and vst residual var gene pca 
DefaultAssay(h1) <- 'RNA'
h1 = NormalizeData(h1) %>% 
  FindVariableFeatures(selection.method = 'vst', verbose = TRUE) %>% 
  ScaleData() %>%
  RunPCA()

# define isotype controls and select non staining proteins and define raw protein data 
rownames(adt_raw) = str_replace_all(rownames(adt_raw),pattern = "_", replacement = "-")
non_staining = c("CX3CR1-PROT", "CD357-PROT", "CD275-PROT", "CD152-PROT", "CD294-PROT",
                 "CD70-PROT", "CD90-PROT", "Mouse IgG2bkIsotype-PROT", "MouseIgG1kappaisotype-PROT",  
                 "RatIgG2bkIsotype-PROT", "CD206-PROT", "MouseIgG2akappaisotype-PROT", "CD223-PROT", 
                 "CD138-PROT", "CD274-PROT", "CD137-PROT", "CD273-PROT","CD80-PROT")
adt_raw = as.matrix(adt_raw[setdiff(rownames(adt_raw), non_staining), ])

# set clr assay as RAW adt then norlaized in seurat
adt_assay = CreateAssayObject(counts = adt_raw)
h1[["CLR"]] <- adt_assay

# run clr Joint clustering workflow 
DefaultAssay(h1) <- 'CLR'
h1 = NormalizeData(h1, normalization.method = 'CLR', margin = 2) %>% ScaleData() 

# hack seurat to use normalized protein values as a dimensionality reduction object.
# these variable veatures and PCs are not used in clustering below
# they are added to prevent the function stop() call when these slots are empty
VariableFeatures(h1) <- rownames(adt_raw)
h1 = RunPCA(h1, reduction.name = 'pCLR', features = VariableFeatures(h1))

# make matrix of norm values to add as dr embeddings
pseudo = t(h1@assays$CLR@data)
p_colnames = paste('pseudo', 1:69, sep = "_")
colnames(pseudo) = p_colnames
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
h1 = FindClusters(h1, graph.name = "clr_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3,
                  resolution = 3, random.seed = 1990,  verbose = FALSE)
h1 = RunUMAP(h1, nn.name = "clr_pseudo_wnn", reduction.name = "CLR_wnn_umap",
             reduction.key = "CLR_wnnUMAP_", seed.use = 1990)

####################################
## Code below repeats analysis above with dsb normalized values. 
# the same PCs from RNA are combined with dsb normalized protein values 
# dsb workflow 

# add raw counts to a new dsb assay 
# the new assay is NOT used 
# below the dsb normalized values are added directly as a dimensionality reduction object 
h1[["dsb"]] <- adt_assay
DefaultAssay(h1) <- 'dsb'

# add dsb normalized matrix to data slot 
rownames(adt_dsb) = str_replace_all(rownames(adt_dsb), pattern = "_", replacement = "-")
stopifnot(isTRUE(all.equal(colnames(adt_dsb), rownames(h1@meta.data))))
adt_dsb_add = as.matrix(adt_dsb[setdiff(rownames(adt_dsb), non_staining), ])
h1@assays$dsb@data = adt_dsb_add

# need to scale the raw data in the dsb slot to prevent the function stop() 
# below, the scaled values are not used, the dsb values are used directly 
h1 = ScaleData(h1, assay = 'dsb')

# add the dsb values directly 
# hack seurat to use normalized protein values as a dimensionality reduction object.
# run pca to prevent the function stop()
VariableFeatures(h1) <- rownames(adt_raw)
h1 = RunPCA(h1, reduction.name = 'pdsb', features = VariableFeatures(h1))

# make matrix of norm values to add as dr embeddings
pseudo = t(h1@assays$dsb@data)
fake_colnames = paste('pseudo', 1:69, sep = "_")
colnames(pseudo) = fake_colnames
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
h1 = FindClusters(h1, graph.name = "dsb_pseudo_knn", n.start = 5, n.iter = 5, algorithm = 3,
                  resolution = 3, random.seed = 1990,  verbose = FALSE)
h1 = RunUMAP(h1, nn.name = "dsb_pseudo_wnn", reduction.name = "dsb_wnn_umap", 
             reduction.key = "dsb_wnnUMAP_", seed.use = 1990)

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
#   [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6      
# [9] ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       cowplot_1.1.1      magrittr_2.0.1     ggridges_0.5.3     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        rprojroot_2.0.2       mclust_5.4.7         
# [7] fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0   leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1        
# [13] RSpectra_0.16-0       lubridate_1.7.9.2     xml2_1.3.2            codetools_0.2-18      splines_4.0.5         polyclip_1.10-0      
# [19] jsonlite_1.7.2        broom_0.7.5           ica_1.0-2             cluster_2.1.2         dbplyr_2.1.0          png_0.1-7            
# [25] uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5        httr_1.4.2           
# [31] backports_1.2.1       assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2        cli_2.5.0            
# [37] limma_3.46.0          later_1.1.0.1         htmltools_0.5.1.1     tools_4.0.5           igraph_1.2.6          gtable_0.3.0         
# [43] glue_1.4.2            RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7       cellranger_1.1.0     
# [49] vctrs_0.3.6           nlme_3.1-152          lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10            
# [55] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1        
# [61] zoo_1.8-8             scales_1.1.1          spatstat.core_2.0-0   hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0 
# [67] parallel_4.0.5        RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15         
# [73] stringi_1.5.3         rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41       ROCR_1.0-11          
# [79] tensor_1.5            patchwork_1.1.1       htmlwidgets_1.5.3     tidyselect_1.1.0      parallelly_1.23.0     RcppAnnoy_0.0.18     
# [85] plyr_1.8.6            R6_2.5.0              generics_0.1.0        DBI_1.1.1             withr_2.4.1           pillar_1.4.7         
# [91] haven_2.3.1           mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10       abind_1.4-5           future.apply_1.7.0   
# [97] modelr_0.1.8          crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3          grid_4.0.5           
# [103] readxl_1.3.1          data.table_1.14.0     reprex_1.0.0          digest_0.6.27         xtable_1.8-4          httpuv_1.5.5         
# [109] munsell_0.5.0         viridisLite_0.3.0    