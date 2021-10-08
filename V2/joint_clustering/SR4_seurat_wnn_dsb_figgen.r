# R4 Seurat 4 (see session info) 
suppressMessages(library(dsb))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# project info 
project_title = "joint_model"
figpath = paste0(here("V2/joint_clustering/"), project_title, "/figures/"); dir.create(figpath, recursive = TRUE)
datapath = paste0(here("V2/joint_clustering/"), project_title, "/generated_data/"); dir.create(datapath, recursive = TRUE)

# read joint analysis script 
s = readRDS(here('V2/joint_clustering/joint_model/generated_data/h1_WNN.rds'))

# set theme for umap 
boxbox = list(theme_bw(),
              theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                    axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),axis.text.y = element_blank()))

# dsb map 
cu = c(BuenColors::jdb_palette(name = 'lawhoops', 20), rep('pink', 19))
p1 = AugmentPlot(DimPlot(s, reduction = 'dsb_wnn_umap', cols = cu, 
                         group.by = 'dsb_pseudo_knn_res.3',
                         label = TRUE, repel = TRUE, label.size = 8)) + 
  boxbox + 
  NoLegend() + 
  xlab('dsb protein RNA multimodal UMAP 1') + 
  ylab('dsb protein RNA multimodal UMAP 2') + 
  ggtitle('dsb WNN')
ggsave(p1, filename = paste0(figpath, "p3_wnncluster_dsb.pdf"), width = 3.5, height = 3.5)


# set up dsb and clr comparison data
prots = rownames(s@assays$dsb@data)
dsb_df = cbind(as.data.frame(t(s@assays$dsb@data)), s@meta.data) %>% 
  group_by(dsb_pseudo_knn_res.3) %>% 
  gather(prot, count, `AnnexinV-PROT`:`CD20-PROT`) %>% 
  group_by(prot, dsb_pseudo_knn_res.3) %>% 
  summarize(median_dsb = median(count), 
  )

clr_df = cbind(as.data.frame(t(s@assays$CLR@data)), s@meta.data) %>% 
  group_by(dsb_pseudo_knn_res.3) %>% 
  gather(prot, count, `AnnexinV-PROT`:`CD20-PROT`) %>% 
  group_by(prot, dsb_pseudo_knn_res.3) %>% 
  summarize(median_clr = median(count), 
  )

# add neg median for each protein  to dsb_df
neg_drop = readRDS(file = here('data/V2_Data/background_data/adt_neg_dmx.rds'))
rownames(neg_drop) = str_replace_all(string = rownames(neg_drop), pattern = "_",replacement = "-")
neg_drop = neg_drop[prots, ]

# calculate medain negative drops 
neg_med = data.frame(median_neg = apply(log10(neg_drop + 1), 1, median) ) %>% rownames_to_column("prot")

# calculate negative drop 98th percentile
q98 = function(x){ quantile(x,probs = 0.98)}
neg_q98 = data.frame(q98_neg = apply(log10(neg_drop + 1), 1, q98) ) %>% rownames_to_column("prot")

# add negative drop log10 median 
d = full_join(dsb_df, neg_med, by = "prot")
d = full_join(dsb_df, neg_q98, by = "prot")

# add clr median 
d$median_clr = clr_df$median_clr

# highlight cluster 3
d2 = d %>% filter(dsb_pseudo_knn_res.3 == '3')
d2$prot = str_replace_all(string = d2$prot,pattern = "-PROT", replacement = "")

# dsb vs empty drop median 
p = ggplot(d2, aes(x = q98_neg, y = median_dsb, label = prot)) + 
  xlab("background drop 98th percentile") +  
  ylab("dsb normalized median") + 
  ggpubr::stat_cor(method = 'pearson',label.x.npc = 0.60, label.y.npc = 0.9) + 
  geom_point() + 
  geom_point(data = d2 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = 'deepskyblue3') + 
  ggrepel::geom_text_repel(data = d2 %>% filter(median_dsb > 5),
                           segment.size = 0.5, size = 3, force = 2, max.overlaps = 20) + 
  geom_hline(yintercept = 3.5, color = 'red') + 
  theme_bw() 
ggsave(p, filename = paste0(figpath, "q98_background_vs_dsb_c3_joint.pdf"), width = 3.5, height = 3.5)

# CLR vs median negative
p = ggplot(d2, aes(x = q98_neg, y = median_clr, label = prot)) + 
  xlab("background drop 98th percentile") +  
  ylab("CLR normalized median") + 
  geom_point() + 
  geom_point(data = d2 %>% filter(median_dsb > 5), shape = 21, size = 2.5, fill = 'grey60') + 
  geom_point(data = d2 %>% filter(median_dsb < 3.5 & median_clr > 1),
             shape = 21, size = 2.5, fill = 'red') + 
  ggpubr::stat_cor(method = 'pearson',label.x.npc = 0.05, label.y.npc = 0.9) + 
  ggrepel::geom_text_repel(data = d2 %>% filter(median_dsb > 5),
                           segment.size = 0.5, size = 3, force = 2) + 
  ggrepel::geom_text_repel(data = d2 %>% filter(median_dsb < 3.5 & median_clr > 1),
                           segment.size = 0.5, size = 2, max.overlaps = 30,nudge_x = 0.3,
                           force = 2) + 
  theme_bw() 
ggsave(p, filename = paste0(figpath, "q98_background_vs_clr_c3_joint.pdf"), width = 3.5, height = 3.5)

# clr vs dsb 
p = ggplot(d2, aes(x = median_clr, y = median_dsb, label = prot)) + 
  ylab("dsb normalized median") +  
  xlab("CLR normalized median") + 
  geom_point() + 
  geom_point(data = d2 %>% filter(median_dsb > 5), shape = 21, size = 2.5, fill = 'grey60') + 
  geom_point(data = d2 %>% filter(median_dsb < 3.5 & median_clr > 1),
             shape = 21, size = 2.5, fill = 'red') + 
  ggrepel::geom_text_repel(data = d2 %>% filter(median_dsb > 5),
                           segment.size = 0.5, size = 3, force = 2) + 
  ggrepel::geom_text_repel(data = d2 %>% filter(median_dsb < 3.5 & median_clr > 1),
                           segment.size = 0.5, size = 2, max.overlaps = 30,nudge_x = 0.3,
                           force = 2) + 
  theme_bw() 
ggsave(p, filename = paste0(figpath, "dsb_vs_clr_c3_joint.pdf"), width = 3.5, height = 3.5)

# de gene cluster 3 
dsb_rna_m = data.table::fread(file = here('git_ignore/seurat_v4_workflow/generated_data/dsb_rna_m.txt'))
p = ggplot(dsb_rna_m %>% filter(cluster == '3' & avg_log2FC > 0), aes(x = avg_log2FC, y = myAUC, label = gene)) +
  theme_bw() + 
  ggrepel::geom_text_repel(data = dsb_rna_m %>% filter( cluster == '3' & myAUC > 0.9),
                           segment.size = 0.5, size = 3, force = 2, max.overlaps = 20) + 
  geom_point() + 
  xlab('average mRNA log2 fold change') + 
  ylab('AUC Cluster 3 classifier')
ggsave(p, filename = paste0(figpath, 'c3_degenes.pdf'), width = 3.5, height = 3.5)

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
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2        tibble_3.0.6      
# [9] ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       magrittr_2.0.1     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1          backports_1.2.1       BuenColors_0.5.6      plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2        splines_4.0.5        
# [8] listenv_0.8.0         scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1     fansi_0.4.2           tensor_1.5            cluster_2.1.2        
# [15] ROCR_1.0-11           openxlsx_4.2.3        limma_3.46.0          globals_0.14.0        modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0
# [22] colorspace_2.0-0      rvest_0.3.6           ggrepel_0.9.1         haven_2.3.1           crayon_1.4.1          jsonlite_1.7.2        spatstat.data_2.1-0  
# [29] survival_3.2-10       zoo_1.8-8             glue_1.4.2            polyclip_1.10-0       gtable_0.3.0          leiden_0.3.7          car_3.0-10           
# [36] future.apply_1.7.0    abind_1.4-5           scales_1.1.1          pheatmap_1.0.12       DBI_1.1.1             rstatix_0.7.0         miniUI_0.1.1.1       
# [43] Rcpp_1.0.6            viridisLite_0.3.0     xtable_1.8-4          reticulate_1.18       spatstat.core_2.0-0   foreign_0.8-81        mclust_5.4.7         
# [50] htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2             pkgconfig_2.0.3       farver_2.0.3         
# [57] uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10         utf8_1.1.4            tidyselect_1.1.0      labeling_0.4.2        rlang_0.4.10         
# [64] reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0         cellranger_1.1.0      tools_4.0.5           cli_2.5.0             generics_0.1.0       
# [71] broom_0.7.5           ggridges_0.5.3        fastmap_1.1.0         goftest_1.2-2         fs_1.5.0              fitdistrplus_1.1-3    zip_2.1.1            
# [78] RANN_2.6.1            pbapply_1.4-3         future_1.21.0         nlme_3.1-152          mime_0.10             xml2_1.3.2            compiler_4.0.5       
# [85] rstudioapi_0.13       plotly_4.9.3          curl_4.3              png_0.1-7             ggsignif_0.6.0        spatstat.utils_2.1-0  reprex_1.0.0         
# [92] stringi_1.5.3         lattice_0.20-41       Matrix_1.3-2          vctrs_0.3.6           pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1  
# [99] lmtest_0.9-38         RcppAnnoy_0.0.18      data.table_1.14.0     cowplot_1.1.1         irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1      
# [106] R6_2.5.0              promises_1.2.0.1      KernSmooth_2.23-18    gridExtra_2.3         rio_0.5.16            parallelly_1.23.0     codetools_0.2-18     
# [113] MASS_7.3-53.1         assertthat_0.2.1      rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           parallel_4.0.5       
# [120] hms_1.0.0             grid_4.0.5            rpart_4.1-15          FSA_0.9.0             carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0         
# [127] shiny_1.6.0           lubridate_1.7.9.2  