suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(tidyverse))

# figpath 
figpath = here('V2/si/figures/de/'); dir.create(figpath)
datapath = here('V2/si/generated_data/')

# Read seurat version 2 object and extract data 
s = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds'))
sadt = s@assay$CITE@raw.data
sdsb = s@assay$CITE@data
srna = s@raw.data
smd = s@meta.data

# init new object Version 4 
s = CreateSeuratObject(counts = srna, meta.data = smd)
adt = CreateAssayObject(counts = sadt)
s[['CITE']] = adt
s[['CLR']] = adt
s = SetAssayData(s,slot = 'data', assay = 'CITE',new.data = sdsb)
DefaultAssay(s) = 'CLR'
s = NormalizeData(object = s, normalization.method = 'CLR', margin = 2)

# define non staining / isotypes
non_staining = c("CX3CR1_PROT", "CD357_PROT", "CD275_PROT", "CD152_PROT", "CD294_PROT",
                 "CD70_PROT", "CD90_PROT", "Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",  
                 "RatIgG2bkIsotype_PROT", "CD206_PROT", "MouseIgG2akappaisotype_PROT", "CD223_PROT", 
                 "CD138_PROT", "CD274_PROT", "CD137_PROT", "CD273_PROT","CD80_PROT")
# define markers to test 
non_staining = str_replace_all(non_staining, pattern = "_", replacement = "-")
prot_test = setdiff(rownames(s@assays$CITE@counts), non_staining)

# define cell types to test as the coarse celltypes from baseline paper (p3dist_1 clusters C0-C9)
Idents(s) = "celltype_label_1"
ct = unique(s@meta.data$celltype_label_1); print(ct) 

# differential expression with log fold change for proteins celltype i vs all 
de = list()
for (i in 1:length(ct)) {
  
  print(ct[i]) 
  # dsb 
  dsb_ = FindMarkers(object = s, ident.1 = ct[i], features = prot_test,
                     logfc.threshold = 0.3, only.pos = TRUE,
                     assay = "CITE", slot = "data") %>%   
    mutate(method = 'dsb') %>% 
    mutate(celltype = ct[i]) %>% 
    rownames_to_column('protein')
  # CLR 
  clr_ = FindMarkers(object = s, ident.1 = ct[i], features = prot_test,  
                     logfc.threshold = 0.3, only.pos = TRUE,
                     assay = "CLR", slot = "data") %>% 
    mutate(method = 'CLR') %>% 
    mutate(celltype = ct[i]) %>% 
    rownames_to_column('protein')
  de[[i]] = rbind(dsb_ , clr_)
}
de_test = do.call(rbind, de)
# remove prot strings from protein names 
de_test$protein = str_replace_all(string = de_test$protein, pattern = "-PROT",replacement = "")
de_test$protein = str_replace_all(string = de_test$protein, pattern = " ",replacement = "")
  
# visualize log fold change differences 
# global plot 
point = list(
  theme_bw(), 
  xlab("Log2 Fold Change"), 
  theme(axis.text.y = element_text(size = 7, color = 'black')), 
  geom_point(alpha = 0.8, shape = 16, size = 2, show.legend = F), 
  ylab(""),
  scale_color_manual(values = c('grey60', 'deepskyblue3'))
)

# shorten name of C2 p3dist 1 for label 
de_test$celltype[de_test$celltype == "Classical Monocytes and mDC"] = "Classical Mono and mDC"
ct[ct == "Classical Monocytes and mDC"] = "Classical Mono and mDC"

# remove CD103 to avoid confusion as these cells are within mem population in p3dist_1
# but they cluster separately in p3_dist_3
de_test = de_test %>% filter(!protein == "CD103")

# save plot 
for (i in 1:length(ct)) {
  d = de_test %>% filter(celltype == ct[i]) %>% filter(p_val_adj < 0.01)
  p = ggplot(d, aes(x = avg_log2FC, y = reorder(protein, avg_log2FC), color = method)) + 
    point + 
    ggtitle(ct[i]) + 
    theme(title = element_text(size = 7)) + theme(aspect.ratio = 3)
  p
  ggsave(p, filename = paste0(figpath,ct[i],'_de.pdf'),width = 1.8, height = 3)
}

# save results 
data.table::fwrite(x = de_test, file = paste0(datapath,'de_test.txt'))
de_test = data.table::fread(file = here('V2/si/generated_data/de_test.txt'))

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
#   [1] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3     
# [9] tidyverse_1.3.0    here_1.0.1         SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        ggridges_0.5.3        rprojroot_2.0.2       fs_1.5.0             
# [8] rstudioapi_0.13       spatstat.data_2.1-0   farver_2.0.3          leiden_0.3.7          listenv_0.8.0         ggrepel_0.9.1         lubridate_1.7.9.2    
# [15] xml2_1.3.2            codetools_0.2-18      splines_4.0.5         polyclip_1.10-0       jsonlite_1.7.2        broom_0.7.5           ica_1.0-2            
# [22] cluster_2.1.2         dbplyr_2.1.0          png_0.1-7             uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0
# [29] compiler_4.0.5        httr_1.4.2            backports_1.2.1       assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2       
# [36] limma_3.46.0          cli_2.5.0             later_1.1.0.1         htmltools_0.5.1.1     tools_4.0.5           igraph_1.2.6          gtable_0.3.0         
# [43] glue_1.4.2            RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6          
# [50] nlme_3.1-152          lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10             miniUI_0.1.1.1        lifecycle_1.0.0      
# [57] irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1         zoo_1.8-8             scales_1.1.1          spatstat.core_2.0-0  
# [64] hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0  parallel_4.0.5        RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3        
# [71] gridExtra_2.3         rpart_4.1-15          stringi_1.5.3         rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41      
# [78] ROCR_1.0-11           tensor_1.5            labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.3     cowplot_1.1.1         tidyselect_1.1.0     
# [85] parallelly_1.23.0     RcppAnnoy_0.0.18      plyr_1.8.6            magrittr_2.0.1        R6_2.5.0              generics_0.1.0        DBI_1.1.1            
# [92] withr_2.4.1           pillar_1.4.7          haven_2.3.1           mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10       abind_1.4-5          
# [99] future.apply_1.7.0    modelr_0.1.8          crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3          grid_4.0.5           
# [106] readxl_1.3.1          data.table_1.14.0     reprex_1.0.0          digest_0.6.27         xtable_1.8-4          httpuv_1.5.5          munsell_0.5.0        
# [113] viridisLite_0.3.0    