# S4 R4 gap statistic 
suppressMessages(library(dsb))
suppressMessages(library(tidyverse))
suppressMessages(library(factoextra))
suppressMessages(library(cluster))
suppressMessages(library(here))

# save paths 
figpath = here("V2/si/figures/"); dir.create(figpath)
datapath = here("V2/si/generated_data/"); dir.create(datapath)

# load S2 object and extract data 
s = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds'))
adt_raw = s@assay$CITE@raw.data
dsb = s@assay$CITE@data
md = s@meta.data
rm(s); gc()

# get CLR normalized ADT matrix (across cells)
clr = 
  Seurat::CreateSeuratObject(counts = adt_raw) %>%  
  Seurat::NormalizeData(normalization.method = 'CLR', margin = 2) %>% 
  Seurat::GetAssayData(slot = 'data')
rownames(clr) = str_replace_all(string = rownames(clr), pattern = "-", replacement = "_")

# remove isotypes from data 
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",
             "RatIgG2bkIsotype_PROT", "MouseIgG2akappaisotype_PROT")
dsb = dsb[setdiff(rownames(dsb), isotypes), ]
clr = clr[setdiff(rownames(clr), isotypes), ]

### Compare Gap statistic or k medoids clustering with dsb and clr normalized values 
# calculate gap statistic 
set.seed(1990)
gap_dsb <- clusGap(t(dsb), FUN = cluster::clara, K.max = 20, B = 20, d.power = 2)
gap_clr <- clusGap(t(clr), FUN = cluster::clara, K.max = 20, B = 20, d.power = 2)

# save output 
saveRDS(object = gap_dsb,file = paste0(datapath, 'gap_dsb.rds'))
saveRDS(object = gap_clr,file = paste0(datapath, 'gap_clr.rds'))

gap_dsb = readRDS(here('V2/si/generated_data/gap_dsb.rds'))
gap_clr = readRDS(here('V2/si/generated_data/gap_clr.rds'))

# compare clr and dsb gap statistic 
p3 = fviz_gap_stat(gap_dsb)
p4 = fviz_gap_stat(gap_clr)
p3$data$method = 'dsb'
p4$data$method = 'clr'
d = rbind(p3$data, p4$data)
cs = seq(from =2, to = 20,by = 2)
p= 
  ggplot(d, aes(x = clusters, y = gap, ymin = ymin, ymax = ymax, color = method, group = method)) + 
  theme_bw() + 
  geom_errorbar() +
  geom_point() + 
  ylim(c(1, 3)) + 
  geom_line() +
  theme(axis.text.x = element_text(size = 6)) + 
  xlab('number of clusters (k)') + ylab('Gap Statistic')+ 
  scale_color_manual(values = c('grey60', 'deepskyblue3')) +
  theme(legend.position = c(0.8, 0.3)) 
ggsave(p,filename = paste0(figpath,'clust_gap_full.pdf'), width = 3, height = 3)

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
#   [1] SeuratObject_4.0.0 Seurat_4.0.1       here_1.0.1         cluster_2.1.2      factoextra_1.0.7   forcats_0.5.1      stringr_1.4.0     
# [8] dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0   
# [15] dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        ggridges_0.5.3        mclust_5.4.7         
# [7] rprojroot_2.0.2       fs_1.5.0              spatstat.data_2.1-0   rstudioapi_0.13       leiden_0.3.7          listenv_0.8.0        
# [13] ggrepel_0.9.1         lubridate_1.7.9.2     xml2_1.3.2            codetools_0.2-18      splines_4.0.5         polyclip_1.10-0      
# [19] jsonlite_1.7.2        broom_0.7.5           ica_1.0-2             dbplyr_2.1.0          png_0.1-7             uwot_0.1.10          
# [25] spatstat.sparse_2.0-0 shiny_1.6.0           sctransform_0.3.2     compiler_4.0.5        httr_1.4.2            backports_1.2.1      
# [31] assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2        limma_3.46.0          cli_2.5.0            
# [37] later_1.1.0.1         htmltools_0.5.1.1     tools_4.0.5           igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
# [43] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6          
# [49] nlme_3.1-152          lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10             miniUI_0.1.1.1       
# [55] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1         zoo_1.8-8            
# [61] scales_1.1.1          spatstat.core_2.0-0   spatstat.utils_2.1-0  hms_1.0.0             promises_1.2.0.1      parallel_4.0.5       
# [67] RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15          stringi_1.5.3        
# [73] rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41       tensor_1.5            ROCR_1.0-11          
# [79] patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1         tidyselect_1.1.0      parallelly_1.23.0     RcppAnnoy_0.0.18     
# [85] plyr_1.8.6            magrittr_2.0.1        R6_2.5.0              generics_0.1.0        DBI_1.1.1             mgcv_1.8-34          
# [91] pillar_1.4.7          haven_2.3.1           withr_2.4.1           fitdistrplus_1.1-3    abind_1.4-5           survival_3.2-10      
# [97] future.apply_1.7.0    modelr_0.1.8          crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3         
# [103] grid_4.0.5            readxl_1.3.1          data.table_1.14.0     reprex_1.0.0          digest_0.6.27         xtable_1.8-4         
# [109] httpuv_1.5.5          munsell_0.5.0         viridisLite_0.3.0   
