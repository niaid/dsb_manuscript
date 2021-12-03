# manual gates 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
'%ni%' = Negate('%in%')

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# source vis script 
source(here("V2/functions/analysis_functions.R"))
# source dsb gate functions 
source(here("V2/dsb_normalize_cluster_pipeline/manual_gates/dsb_gates.r"))

# data set-up. 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
h1 = SubsetData(h1, subset.name = "CD11c_PROT", accept.low = -10)
h1 = SubsetData(h1, subset.name = "CD14_PROT", accept.low = -4, accept.high = 12)

# lineage 
df = h1@assay$CITE@data %>% t %>% as.data.frame()
p = ggplot(df, aes(x = CD19_PROT, y = CD3_PROT )) + 
  geom_bin2d(bins = 300)+
  scale_fill_viridis_c(option = "B") + 
  geom_density_2d(color = "white")+
  geom_hline(yintercept = 5) +
  geom_vline(xintercept =8)
ggsave(p, filename = paste0(figpath, "dsb_lineage.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath, "dsb_lineage.pdf"), width = 6, height = 5)

# lineage negative and outlier cell removal for plots
linneg = GateLinneg(h1, return.seurat = T)
linneg = SubsetData(linneg, subset.name = "CD14_PROT",accept.high = 14)


# plot 
p = GenePlot4(linneg, gene1 = "CD14_PROT", gene2 = "CD16_PROT", pt.size = 0.5) +
  geom_hline(yintercept = 4) + geom_vline(xintercept = 2)
ggsave(p, filename = paste0(figpath, "dsb_linneg_cd14_cd16.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath, "dsb_linneg_cd14_cd16.pdf"), width = 6, height = 5)

# Classical monocytes 
cd14 = GateMono14(linneg, return.seurat = T)
subset1 = SubsetData(linneg, cells.use = linneg@cell.names[!linneg@cell.names %in% cd14@cell.names])

# plot 
p = GenePlot4(subset1, gene1 = "CD56_PROT", gene2 = "CD16_PROT", pt.size = 0.5)
ggsave(p, filename = paste0(figpath,"dsb_subset1_cd56_cd16.png"),  width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_subset1_cd56_cd16.pdf"),  width = 6, height = 5)

# dump gate 
neg16neg56 = Gateneg16neg56(subset1,return.seurat = T)
subset2 = SubsetData(subset1, cells.use = subset1@cell.names[!(subset1@cell.names %in% neg16neg56@cell.names)])

# non classical and intermediate
p = GenePlot4(subset2, gene1 = "CD11c_PROT", gene2 = "CD14_PROT", pt.size = 0.5) +
  geom_hline(yintercept = 2.5) + geom_vline(xintercept = 9)
ggsave(p, filename = paste0(figpath,"dsb_subset2_cd11c_cd14.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_subset2_cd11c_cd14.pdf"), width = 6, height = 5)


#### T cell gating 
linneg = GateLinnegT(h1, return.seurat = T)
df = linneg@assay$CITE@data %>% t %>% as.data.frame()
p = ggplot(df, aes(x = CD4_PROT, y = CD8_PROT )) + 
  geom_bin2d(bins = 300)+
  scale_fill_viridis_c(option = "B") + 
  geom_density_2d(color = "white")+
  geom_hline(yintercept = 5) +
  geom_vline(xintercept =6.5)
ggsave(p, filename = paste0(figpath, "dsb_Tcells.pdf"), width = 6, height = 5)

# cd4 and CD8 T cells 
cd4 = GateT4(linneg, return.seurat = T)
cd8 = GateT8(linneg, return.seurat = T)
p =GenePlot4(cd4, gene1 = "CD45RO_PROT", gene2 = "CD62L_PROT", pt.size = 0.5) +
  geom_hline(yintercept = 6) + geom_vline(xintercept = 6)
ggsave(p, filename = paste0(figpath,"dsb_cd4.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_cd4.pdf"), width = 6, height = 5)

p =GenePlot4(cd8, gene1 = "CD45RO_PROT", gene2 = "CD62L_PROT", pt.size = 0.5) +
  geom_hline(yintercept = 6) + geom_vline(xintercept = 6)
ggsave(p, filename = paste0(figpath,"dsb_cd8.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_cd8.pdf"), width = 6, height = 5)

# B cell 
bc = GateBC(h1, return.seurat = T)
p =GenePlot4(bc, gene1 = "IgD_PROT", gene2 = "CD27_PROT", pt.size = 0.5) +
  geom_hline(yintercept = 4) + geom_vline(xintercept = 10)
ggsave(p, filename = paste0(figpath,"Bcell.pdf"), width = 6, height = 5)

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
#   [1] viridis_0.5.1     viridisLite_0.3.0 here_0.1          magrittr_2.0.1    forcats_0.4.0     stringr_1.4.0     dplyr_0.8.5       purrr_0.3.3      
# [9] readr_1.3.1       tidyr_1.0.2       tibble_2.1.1      tidyverse_1.2.1   Seurat_2.3.4      Matrix_1.2-15     cowplot_0.9.4     ggplot2_3.1.1    
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5       
# [8] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7        
# [15] lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0         
# [22] robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6        packrat_0.5.0       broom_0.5.2         ica_1.0-2          
# [29] cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3      httr_1.4.0          backports_1.1.4    
# [36] assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0           lars_1.2            acepack_1.4.1       htmltools_0.3.6     tools_3.5.3        
# [43] igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0   
# [50] vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11        
# [57] lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8     
# [64] MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2 
# [71] reticulate_1.12     pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3      
# [78] foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5        
# [85] pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          labeling_0.3       
# [92] htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          R6_2.4.0            generics_0.0.2      snow_0.4-3         
# [99] gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14
# [106] mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0        
# [113] KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7     
# [120] R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0      

