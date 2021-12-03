# R4 Seurat 4
suppressMessages(library(Seurat))
suppressMessages(library(here))

# Read seurat version 2 object and extract data 
s = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds'))
rna_raw = s@raw.data
adt_raw = s@assay$CITE@raw.data
md = s@meta.data

# create S4 object and normalize with CLR across cells  
s = CreateSeuratObject(counts = rna_raw, meta.data = md)
adt_assay = CreateAssayObject(counts = adt_raw)
s[["CLR"]] = adt_assay
DefaultAssay(s) = 'CLR'
s = NormalizeData(s, normalization.method = 'CLR', margin = 2) 

# save
saveRDS(as.matrix(s@assays$CLR@data), file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_CLR_acroass_cells_matrix.rds"))

##########################################
### note this version not used in publication: 
# repeat for unstained cells (load Seurat V2 object)
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))
u_adt = CreateAssayObject(counts = un@assay$CITE@raw.data)
un = CreateSeuratObject(counts = u_adt, assay = "CITE")
un = NormalizeData(object = un, normalization.method = "CLR", margin = 2)
saveRDS(un@assays$CITE@data, file = here("V2/dsb_normalize_cluster_pipeline/generated_data/unstainedcells_CLR_acroass_cells_matrix.rds"))

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
#   [1] here_1.0.1         SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-152          matrixStats_0.58.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.18      RColorBrewer_1.1-2    httr_1.4.2           
# [7] rprojroot_2.0.2       sctransform_0.3.2     tools_4.0.5           R6_2.5.0              irlba_2.3.3           rpart_4.1-15         
# [13] KernSmooth_2.23-18    uwot_0.1.10           mgcv_1.8-34           DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-0     
# [19] tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.5        plotly_4.9.3          scales_1.1.1          lmtest_0.9-38        
# [25] spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.4-3         goftest_1.2-2         stringr_1.4.0         digest_0.6.27        
# [31] spatstat.utils_2.1-0  pkgconfig_2.0.3       htmltools_0.5.1.1     parallelly_1.23.0     fastmap_1.1.0         htmlwidgets_1.5.3    
# [37] rlang_0.4.10          shiny_1.6.0           generics_0.1.0        zoo_1.8-8             jsonlite_1.7.2        ica_1.0-2            
# [43] dplyr_1.0.4           magrittr_2.0.1        patchwork_1.1.1       Matrix_1.3-2          Rcpp_1.0.6            munsell_0.5.0        
# [49] abind_1.4-5           reticulate_1.18       lifecycle_1.0.0       stringi_1.5.3         MASS_7.3-53.1         Rtsne_0.15           
# [55] plyr_1.8.6            grid_4.0.5            parallel_4.0.5        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
# [61] crayon_1.4.1          miniUI_0.1.1.1        deldir_0.2-10         lattice_0.20-41       cowplot_1.1.1         splines_4.0.5        
# [67] tensor_1.5            pillar_1.4.7          igraph_1.2.6          spatstat.geom_2.0-1   future.apply_1.7.0    reshape2_1.4.4       
# [73] codetools_0.2-18      leiden_0.3.7          glue_1.4.2            data.table_1.14.0     vctrs_0.3.6           png_0.1-7            
# [79] httpuv_1.5.5          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.0-0   polyclip_1.10-0      
# [85] tidyr_1.1.2           scattermore_0.7       future_1.21.0         assertthat_0.2.1      ggplot2_3.3.3         mime_0.10            
# [91] xtable_1.8-4          later_1.1.0.1         survival_3.2-10       viridisLite_0.3.0     tibble_3.0.6          cluster_2.1.2        
# [97] globals_0.14.0        fitdistrplus_1.1-3    ellipsis_0.3.1        ROCR_1.0-11    