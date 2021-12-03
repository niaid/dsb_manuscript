# distribution of data with different normalization methods 
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

# unstained control seurat (v2) object
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))
# unnormalized seurat (v2) object 
h1 = readRDS(file = here('data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds'))
# empty drop adt matrices (hashing negatives)
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# define adt matrices 
raw_adt_unstained = as.matrix(un@assay$CITE@raw.data)
raw_adt_background_drops = as.matrix(do.call(cbind, neg_adt))
raw_adt_cells = as.matrix(h1@assay$CITE@raw.data)

# use minimal matching metadata in object creation
cols = intersect(colnames(h1@meta.data), colnames(un@meta.data))
cellmd = h1@meta.data[ ,cols]
unstainedmd = un@meta.data[ ,cols]

# merge cells and unstained metadata 
md_full = rbind(cellmd, unstainedmd)

# define the droplet class as either a stained cell or unstained control cell
md_full$drop_class = if_else(md_full$barcode_check %in% rownames(un@meta.data),
                             'unstained_control', 'stained_cell')

# combine rna (to initialize object) and prot for unstained cells and stained cells, not neg.  
md_full = md_full %>%
  remove_rownames() %>% 
  column_to_rownames('barcode_check')
rna_full = cbind(h1@raw.data, un@raw.data)
adt_full = cbind(raw_adt_cells, raw_adt_unstained)

# make combined (Seurat 4) object min cells is only a gene filter (RNA) for reducing object size
s = CreateSeuratObject(counts = rna_full, min.cells = 100, meta.data = md_full)

# create an assay object for adt 
adt_assay_raw = CreateAssayObject(counts = adt_full)

# normalize
# clr across proteins 
s[['CLRprot']] = adt_assay_raw
s = NormalizeData(s, assay = 'CLRprot',normalization.method = 'CLR',margin = 1)

# clr across cells 
s[['CLRcells']] = adt_assay_raw
s = NormalizeData(s, assay = 'CLRcells',normalization.method = 'CLR',margin = 2)

# library size scaling factors 
s[['libsize']] = adt_assay_raw
s = NormalizeData(s, assay = 'libsize', normalization.method = 'LogNormalize', scale.factor = 1e4)

# log transformation 
s[['log']] = adt_assay_raw
log_transform = log10(1 + adt_full)
# format the externally transformed matrix to match protein _PROT -> -PROT string coersion in Seurat 
rownames(log_transform) = str_replace(string = rownames(log_transform),
                                      pattern = "_", replacement = "-")
s = SetAssayData(object = s, slot = 'data',assay = 'log',new.data = log_transform)

# dsb
s[['dsb']] = adt_assay_raw
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")
# run dsb norm with both batches of cells and unstained cells merged 
dsb = DSBNormalizeProtein(cell_protein_matrix = adt_full, 
                          empty_drop_matrix = raw_adt_background_drops,
                          denoise.counts = TRUE,
                          use.isotype.control = TRUE,
                          isotype.control.name.vec = isotypes)
rownames(dsb) = str_replace(string = rownames(dsb),pattern = "_",replacement = "-")
s = SetAssayData(object = s, slot = 'data',assay = 'dsb',new.data = dsb)

# dsb both batches merged without denoising 
s[['dsbnondenoised']] = adt_assay_raw
dsb_nodenoise = DSBNormalizeProtein(cell_protein_matrix = adt_full, 
                                    empty_drop_matrix = raw_adt_background_drops,
                                    denoise.counts = FALSE)
rownames(dsb_nodenoise) = str_replace(string = rownames(dsb_nodenoise),
                                      pattern = "_",replacement = "-")
s = SetAssayData(s,slot = 'data', assay = 'dsbnondenoised',new.data = dsb_nodenoise)

# dsb with multi batch processing 
# batch split matrices using cels across stained and unstained controls in each batch 
adt_b1 = adt_full[ ,rownames(md_full[md_full$batch == '1', ])]
adt_b2 = adt_full[ ,rownames(md_full[md_full$batch == '2', ])]

# combine the 2 batches into a list indexed by batch 
batchlist = list(adt_b1, adt_b2)

# normalize by batch with the corresponding empty drops from each batch adt_neg
dsb_multi = list()
for (i in 1:length(neg_adt)) {
  dsb_multi[[i]] =  DSBNormalizeProtein(cell_protein_matrix = batchlist[[i]], 
                                        empty_drop_matrix = neg_adt[[i]], 
                                        denoise.counts = TRUE, 
                                        use.isotype.control = TRUE, 
                                        isotype.control.name.vec = isotypes)
}
# merge and add to object 
dsb_multi_ = do.call(cbind, dsb_multi)
dsb_multi_ = as.matrix(dsb_multi_)
rownames(dsb_multi_) = str_replace(string = rownames(dsb_multi_), pattern = "_",replacement = "-")
dsb_multi_2 = dsb_multi_[ ,match(x = rownames(md_full), table = colnames(dsb_multi_))]
s[['dsbdenoisedmultibatch']] = adt_assay_raw
s = SetAssayData(object = s,assay = 'dsbdenoisedmultibatch',slot = 'data',new.data = dsb_multi_2)

# repeat dsb multibatch norm without step II 
dsb_multi_nd = list()
for (i in 1:length(neg_adt)) {
  dsb_multi_nd[[i]] =  DSBNormalizeProtein(cell_protein_matrix = batchlist[[i]], 
                                           empty_drop_matrix = neg_adt[[i]], 
                                           denoise.counts = FALSE)
}
# merge and add to object 
dsb_multi_nd = do.call(cbind, dsb_multi_nd)
dsb_multi_nd = as.matrix(dsb_multi_nd)
rownames(dsb_multi_nd) = str_replace(string = rownames(dsb_multi_nd), pattern = "_",replacement = "-")
dsb_multi_nd = dsb_multi_nd[ ,match(x = rownames(md_full), table = colnames(dsb_multi_nd))]
s[['dsbNONdenoisedmultibatch']] = adt_assay_raw
s = SetAssayData(object = s,assay = 'dsbNONdenoisedmultibatch',slot = 'data',new.data = dsb_multi_nd)

##################
# Visualization of protein distributions 
mg_layer = list(theme_bw(),  
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                geom_vline(xintercept = 0, size = 0.7, color = "black") ,
                geom_hline(yintercept = 0, size = 0.7, color = "black") , 
                geom_point(size = 0.1, shape = 16, alpha = 0.4) ,  
                theme(plot.title = element_text(size = 18) )
)

norms = names(s@assays)[-1]
plist = list()
for (i in 1:length(norms)) {
  
  dat = 
    GetAssayData(object = s,slot = 'data',assay = norms[i]) %>% 
    as.matrix() %>% 
    t() %>% 
    as.data.frame()
  dat = cbind(s@meta.data, dat)
  
  ### Main plot showin in figure 
  plist[[i]] = ggplot(dat, aes(x= `CD4-PROT`, y = `CD14-PROT`)) + 
    mg_layer + 
    ggtitle(norms[i]) + 
    # label the unstained controls with 2d density 
    geom_density_2d(data = dat %>% filter(drop_class == 'unstained_control'),
                    color = 'red',
                    alpha = 0.9,
                    size = 0.95)

}
# adjust aspect ratio and center on main distributions 
# outlier subsets remove less than 20 cells out of 53974
plist[[7]] = plist[[7]] + ylim(c(-5,12)) + theme(aspect.ratio = 0.85)
plist[[8]] = plist[[8]] + ylim(c(-5,12)) + theme(aspect.ratio = 0.85)
plist[[2]] = plist[[2]] + ylim(c(-0.5,3)) + xlim(c(-0.6,5.1)) + theme(aspect.ratio = 0.85)
plist[[3]] = plist[[3]] + ylim(c(-0.5,6)) + xlim(c(-0.6,8)) + theme(aspect.ratio = 0.85)
p = cowplot::plot_grid(plist[[7]] ,  plist[[8]] , plist[[2]] , plist[[3]], nrow = 2)
ggsave(p, filename = paste0(figpath, 'MAIN_dsbmerged_normalization_unstained.png'), width = 6, height = 6)

# save plot list object 
saveRDS(object = plist, file = paste0(datapath,'plist.rds'))


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
#   [1] dsb_0.1.0          here_1.0.1         magrittr_2.0.1     forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4       
# [7] purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0   
# [13] SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        ggridges_0.5.3       
# [6] mclust_5.4.7          rprojroot_2.0.2       fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0  
# [11] farver_2.0.3          leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1         lubridate_1.7.9.2    
# [16] xml2_1.3.2            codetools_0.2-18      splines_4.0.5         polyclip_1.10-0       jsonlite_1.7.2       
# [21] broom_0.7.5           ica_1.0-2             cluster_2.1.2         dbplyr_2.1.0          png_0.1-7            
# [26] uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5       
# [31] httr_1.4.2            backports_1.2.1       assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0        
# [36] lazyeval_0.2.2        limma_3.46.0          cli_2.5.0             later_1.1.0.1         htmltools_0.5.1.1    
# [41] tools_4.0.5           igraph_1.2.6          gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
# [46] reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6          
# [51] nlme_3.1-152          lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10            
# [56] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0        
# [61] MASS_7.3-53.1         zoo_1.8-8             scales_1.1.1          spatstat.core_2.0-0   hms_1.0.0            
# [66] promises_1.2.0.1      spatstat.utils_2.1-0  parallel_4.0.5        RColorBrewer_1.1-2    reticulate_1.18      
# [71] pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15          stringi_1.5.3         rlang_0.4.10         
# [76] pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41       ROCR_1.0-11           tensor_1.5           
# [81] labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1         tidyselect_1.1.0     
# [86] parallelly_1.23.0     RcppAnnoy_0.0.18      plyr_1.8.6            R6_2.5.0              generics_0.1.0       
# [91] DBI_1.1.1             withr_2.4.1           pillar_1.4.7          haven_2.3.1           mgcv_1.8-34          
# [96] fitdistrplus_1.1-3    survival_3.2-10       abind_1.4-5           future.apply_1.7.0    modelr_0.1.8         
# [101] crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3          isoband_0.2.3        
# [106] grid_4.0.5            readxl_1.3.1          data.table_1.14.0     reprex_1.0.0          digest_0.6.27        
# [111] xtable_1.8-4          httpuv_1.5.5          munsell_0.5.0         viridisLite_0.3.0     
