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
  d = de_test %>% filter(celltype == ct[i])
  p = ggplot(d, aes(x = avg_log2FC, y = reorder(protein, avg_log2FC), color = method)) + 
    point + 
    ggtitle(ct[i]) + 
    theme(title = element_text(size = 7))
  ggsave(p, filename = paste0(figpath,ct[i],'_de.pdf'),width = 2, height = 3)
}

# save results 
data.table::fwrite(x = de_test, file = paste0(datapath,'de_test.txt'))

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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cluster_2.1.2            factoextra_1.0.7         dsb_0.1.0                variancePartition_1.20.0 Biobase_2.50.0          
# [6] BiocGenerics_0.36.1      scales_1.1.1             BiocParallel_1.24.1      limma_3.46.0             magrittr_2.0.1          
# [11] mclust_5.4.7             forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4             
# [16] readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3            tidyverse_1.3.0         
# [21] here_1.0.1               SeuratObject_4.0.0       Seurat_4.0.1            
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2       
# [6] splines_4.0.5         listenv_0.8.0         scattermore_0.7       digest_0.6.27         foreach_1.5.1        
# [11] htmltools_0.5.1.1     tensor_1.5            doParallel_1.0.16     ROCR_1.0-11           openxlsx_4.2.3       
# [16] globals_0.14.0        modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0 prettyunits_1.1.1    
# [21] colorspace_2.0-0      rvest_0.3.6           ggrepel_0.9.1         haven_2.3.1           crayon_1.4.1         
# [26] jsonlite_1.7.2        lme4_1.1-26           spatstat.data_2.1-0   survival_3.2-10       zoo_1.8-8            
# [31] iterators_1.0.13      glue_1.4.2            polyclip_1.10-0       gtable_0.3.0          leiden_0.3.7         
# [36] car_3.0-10            future.apply_1.7.0    abind_1.4-5           DBI_1.1.1             rstatix_0.7.0        
# [41] miniUI_0.1.1.1        Rcpp_1.0.6            progress_1.2.2        viridisLite_0.3.0     xtable_1.8-4         
# [46] reticulate_1.18       spatstat.core_2.0-0   foreign_0.8-81        htmlwidgets_1.5.3     httr_1.4.2           
# [51] gplots_3.1.1          RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2             farver_2.0.3         
# [56] pkgconfig_2.0.3       uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10         labeling_0.4.2       
# [61] tidyselect_1.1.0      rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0        
# [66] cellranger_1.1.0      tools_4.0.5           cli_2.5.0             generics_0.1.0        broom_0.7.5          
# [71] ggridges_0.5.3        fastmap_1.1.0         goftest_1.2-2         fs_1.5.0              fitdistrplus_1.1-3   
# [76] zip_2.1.1             caTools_1.18.1        RANN_2.6.1            pbapply_1.4-3         future_1.21.0        
# [81] nlme_3.1-152          mime_0.10             xml2_1.3.2            pbkrtest_0.5-0.1      compiler_4.0.5       
# [86] rstudioapi_0.13       curl_4.3              plotly_4.9.3          png_0.1-7             ggsignif_0.6.0       
# [91] spatstat.utils_2.1-0  reprex_1.0.0          statmod_1.4.35        stringi_1.5.3         lattice_0.20-41      
# [96] Matrix_1.3-2          nloptr_1.2.2.2        ggsci_2.9             vctrs_0.3.6           pillar_1.4.7         
# [101] lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18      data.table_1.14.0    
# [106] cowplot_1.1.1         bitops_1.0-6          irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1      
# [111] colorRamps_2.3        R6_2.5.0              promises_1.2.0.1      rio_0.5.16            KernSmooth_2.23-18   
# [116] gridExtra_2.3         parallelly_1.23.0     codetools_0.2-18      boot_1.3-27           MASS_7.3-53.1        
# [121] gtools_3.8.2          assertthat_0.2.1      rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2    
# [126] mgcv_1.8-34           hms_1.0.0             grid_4.0.5            rpart_4.1-15          minqa_1.2.4          
# [131] carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0          shiny_1.6.0           lubridate_1.7.9.2 