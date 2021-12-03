# run umap on DSB normalized values 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# load python virtual environment for python umap learn via reticulate 
suppressMessages(library(reticulate))
use_virtualenv("r-reticulate")
library(umap)

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# DSB protein data matrix without isotype controls 
h1 = readRDS("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds")
proteins = rownames(h1@assay$CITE@data)
proteins = proteins[-c(67:70)]
adt = h1@assay$CITE@data[proteins, ]
adt = t(adt)

# Set parameters and run umap 
config = umap.defaults
config$n_neighbors = 35
config$min_dist = 0.6
ump = umap(adt,config = config)

# extract cell embeddings 
embeddings = ump$layout %>% as.data.frame()
colnames(embeddings) = c("UMAP1", "UMAP2")

# save combined data 
df = cbind(embeddings, h1@meta.data, as.data.frame(t(h1@assay$CITE@data)))
saveRDS(df, file = paste0(datapath, "h1_sng_metadata_umapdim_prots_dataframe.rds"))

sesessionInfo()
# R version 3.5.3 Patched (2019-03-11 r77192)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] umap_0.2.3.1    reticulate_1.12 here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1   
# [12] tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5        htmlTable_1.13.1    base64enc_0.1-3    
# [10] rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16   
# [19] splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6        packrat_0.5.0       broom_0.5.2        
# [28] ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1   
# [37] lazyeval_0.2.2      cli_1.1.0           lars_1.2            acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1         
# [46] RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10   
# [55] fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8     
# [64] MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  pbapply_1.4-0       gridExtra_2.3      
# [73] rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0      
# [82] SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          htmlwidgets_1.3    
# [91] bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0        
# [100] pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4       
# [109] crayon_1.3.4        hdf5r_1.2.0         KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7     
# [118] R.utils_2.8.0       openssl_1.4         stats4_3.5.3        munsell_0.5.0       askpass_1.1        