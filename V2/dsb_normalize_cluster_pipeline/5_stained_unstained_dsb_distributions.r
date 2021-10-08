# makes visualizations of dsb, umap results used in Figure 2 and  norm comparison in in Fig 1
set.seed(1)
'%ni%' = Negate('%in%')
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(dsb))

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# unstained control dsb norm with empty drops 
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))
b1cells = un@meta.data[un@meta.data$batch == 1, ]$barcode_check
b2cells = un@meta.data[un@meta.data$batch == 2, ]$barcode_check
pos_adt1 = un@assay$CITE@raw.data[ ,b1cells]
pos_adt2 = un@assay$CITE@raw.data[ ,b2cells]
pos_adt = c(pos_adt1, pos_adt2)

# apply denoised scaled by background protein normalization per batch 
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")
un_dsb = list()
for (i in 1:length(neg_adt)) {
  un_dsb[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                       empty_drop_matrix = neg_adt[[i]], 
                                       denoise.counts = TRUE, 
                                       use.isotype.control = TRUE, 
                                       isotype.control.name.vec = isotypes)
}
un_dsb = do.call(cbind, un_dsb)
un = SetAssayData(un, assay.type = "CITE",slot = "data",new.data = un_dsb)

###### nondenoised 
un_dsb_nd = list()
for (i in 1:length(neg_adt)) {
  un_dsb_nd[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                     empty_drop_matrix = neg_adt[[i]], 
                                     denoise.counts = FALSE)
}
un_dsb_nd = do.call(cbind, un_dsb_nd)
un = SetAssayData(un, assay.type = "nond",slot = "data",new.data = un_dsb_nd)

# load raw and dsb normalized h1 data 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
raw = h1@assay$CITE@raw.data
dsb = h1@assay$CITE@data
md = h1@meta.data

# CLR normalize data across genes (original implementation)
h1 = SetAssayData(h1,assay.type = "CLR", slot = "raw.data",new.data = h1@assay$CITE@raw.data)
h1 = NormalizeData(h1, normalization.method = "genesCLR", assay.type = "CLR")

# repeat for unstained 
un = SetAssayData(un,assay.type = "CLR", slot = "raw.data",new.data = un@assay$CITE@raw.data)
un = NormalizeData(un, normalization.method = "genesCLR", assay.type = "CLR")

# add other comparison transformation / normalization methods 
# normalize / transform both unstained spike in and h1 stained cells 
# log10 transform 
log10 = log10(as.matrix(h1@assay$CITE@raw.data + 1))
log10u = log10(as.matrix(un@assay$CITE@raw.data + 1))

# library size scaling factors with log transform 
lognorm = NormalizeData(h1, assay.type = "CITE",normalization.method = "LogNormalize", scale.factor = 1e4)
lognormu = NormalizeData(un, assay.type = "CITE",normalization.method = "LogNormalize", scale.factor = 1e4)

# arcsin transform
arcsin = asinh(as.matrix(h1@assay$CITE@raw.data))
arcsinu = asinh(as.matrix(un@assay$CITE@raw.data))

# square root transformation 
sqr = sqrt(as.matrix(h1@assay$CITE@raw.data))
sqru = sqrt(as.matrix(un@assay$CITE@raw.data))

# CLR acrlss genes 
clr = h1@assay$CLR@data
clru = un@assay$CLR@data

# dsb without step 2 cell-cell variation removal (denoise = false)
dsb_nodenoise = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/NonDenoised_dsb_Mtx.rds"))
dsb_nodenoiseu = un_dsb_nd

# CLR across cells as used in updated versions of Seurat from script 1.1_calc_clr_across_cells_s4.r 
cells_clr = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_CLR_acroass_cells_matrix.rds'))
cells_clru = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/unstainedcells_CLR_acroass_cells_matrix.rds'))

# rm outlier cells 
h1 = SubsetData(h1, subset.name = "CD11c_PROT", accept.low = -10)
h1 = SubsetData(h1, subset.name = "CD14_PROT", accept.low = -4, accept.high = 12)
cells_plot = rownames(h1@meta.data)

# outlier cells removed 
dim(md)[1] - dim(h1@meta.data)[1]
#  58 cells 

# manual gate plot layers 
mg_layer = list(theme_bw(),  
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                geom_point(size = 0.1, shape = 16, alpha = 0.4) ,  
                geom_vline(xintercept = 0, size = 0.9, color = "black") ,
                geom_hline(yintercept = 0, size = 0.9, color = "black") , 
                theme(plot.title = element_text(size = 18) )
)

# arcsin
p1 = ggplot(as.data.frame(t(arcsin))[cells_plot, ], aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +
  ggtitle("Arcsin transformation \nasinh(x)") + xlim(-0.7, 8) + ylim(-0.7, 6) + 
  theme(plot.title = element_text(size = 12)) + 
  geom_density_2d(data = as.data.frame(t(arcsinu)), color = 'red', alpha = 0.9, size = 1)

# log norm global library size scaling (as in RNA) for comparison
p2 = ggplot(as.data.frame(t(as.matrix(lognorm@assay$CITE@data)))[cells_plot, ], aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +
  ggtitle("library size scaling\nlog(1 + x/sum(x) ) * 1e4") + 
  theme(plot.title = element_text(size = 12)) + 
  geom_density_2d(data = as.data.frame(t(as.matrix(lognormu@assay$CITE@data))),color = 'red', alpha = 0.9, size = 1)

# log 10 
p3 = ggplot(as.data.frame(t(log10))[cells_plot, ],  aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +
  ggtitle("Log transformation\nlog10(1 + x)") + xlim(-0.5, 3) + ylim(-0.5, 2.5) + 
  theme(plot.title = element_text(size = 12)) + 
  geom_density_2d(data = as.data.frame(t(log10u)),color = 'red', alpha = 0.9, size = 1)

# CLR 
p4 = ggplot(as.data.frame(t(h1@assay$CLR@data))[cells_plot, ], aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + 
  ggtitle("CLR \nTransformation (across proteins)") + xlim(-0.5, 3) + ylim(-0.5, 3) + 
  theme(plot.title = element_text(size = 12)) + 
  geom_density_2d(data = as.data.frame(t(un@assay$CLR@data)),color = 'red', alpha = 0.9, size = 1)

# dsb
p5 = ggplot(as.data.frame(t(h1@assay$CITE@data))[cells_plot, ], aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer + 
  ggtitle("dsb: both protein-specific ambient noise & \ncell-intrinsic technical component removed") +
  theme(plot.title = element_text(size = 11)) + 
  geom_density_2d(data = as.data.frame(t(un@assay$CITE@data)), color = 'red', alpha = 0.9, size = 1)

# dsb non-denoised
p6 = ggplot(as.data.frame(t(dsb_nodenoise))[cells_plot, ], aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer + 
  ggtitle("dsb: \nonly protein-specific ambient noise removed") + xlim(c(-3, 13)) +  ylim(c(-3,13)) + 
  theme(plot.title = element_text(size = 10)) + 
  geom_density_2d(data = as.data.frame(t(un_dsb_nd)), color = 'red', alpha = 0.9, size = 1)

# clr across cells
p7 =  ggplot(as.data.frame(t(cells_clr))[cells_plot, ], aes(x = `CD4-PROT`, y = `CD14-PROT`)) +
  mg_layer + 
  ggtitle("CLR Transformation\n (across cells)") + 
  theme(plot.title = element_text(size = 12)) + 
  geom_density_2d(data = as.data.frame(t(cells_clru)), color = 'red', alpha = 0.9, size = 1)

# combine 
p = cowplot::plot_grid(p5, p6, p4, p2, p3, p7,nrow = 2)

# save 
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained2.png"), width = 4.75, height = 3.25, dpi = "retina", scale = 2)
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained.png"), width = 9.5, height = 6.5, dpi = "retina")
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained.pdf"), width = 12.5, height = 8.5)
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained.png"), width = 13, height = 9)

sessionInfo()
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
#   [1] dsb_0.1.0       here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1   
# [11] tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5        htmlTable_1.13.1   
# [9] base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10     
# [17] xml2_1.2.0          codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3      
# [25] jsonlite_1.6        packrat_0.5.0       broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0        
# [33] compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      limma_3.38.3        cli_1.1.0           lars_1.2           
# [41] acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3     
# [49] Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10    fpc_2.2-1          
# [57] gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8     
# [65] MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12    
# [73] pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3    
# [81] caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1     
# [89] bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4         
# [97] R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71     
# [105] withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4       
# [113] hdf5r_1.2.0         KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7     
# [121] R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0 