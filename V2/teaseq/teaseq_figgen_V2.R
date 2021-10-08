# R4 Seurat 4 (see session info)  
suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# project info 
project_title = "TEA-seq"
figpath = paste0(here("V2/teaseq/figures_v2/"), project_title, "/")
datapath = paste0(here("V2/teaseq/generated_data_v2/"), project_title, "/")

# load processed object and data 
s = readRDS(file = here('V2/teaseq/generated_data_v2/TEA-seq/full_teaseq_r4s4_object_processed.rds'))
adt_background = readRDS(file = here('V2/teaseq/generated_data_v2/TEA-seq/adt_background.rds'))

# set theme for umap 
boxbox = list(
  theme_bw(), 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())
  )

# clr map 
DefaultAssay(s) <- 'CLR'
p = AugmentPlot(DimPlot(s, reduction = 'CLR_wnn_umap', group.by = 'clr_pseudo_knn_res.2',
                         label = TRUE, repel = TRUE, label.size = 7.5)) + 
  boxbox + 
  NoLegend() + 
  xlab('CLR protein RNA multimodal UMAP 1') + 
  ylab('CLR protein RNA multimodal UMAP 2') +
  ggtitle("TEA-seq CLR normalized protein")
ggsave(p, filename = paste0(figpath, "p3_wnncluster_clr.pdf"), width = 3.5, height = 3.5)

# dsb map 
DefaultAssay(s) <- 'CITE'
p = AugmentPlot(DimPlot(s, reduction = 'dsb_wnn_umap', group.by = 'dsb_pseudo_knn_res.2',
                         label = TRUE, repel = TRUE, label.size = 7.5, cols = BuenColors::jdb_palette("corona"))) + 
  boxbox + 
  NoLegend() + 
  xlab('dsb protein RNA multimodal UMAP 1') + 
  ylab('dsb protein RNA multimodal UMAP 2') +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) + 
  ggtitle("TEA-seq dsb normalized protein")
ggsave(p, filename = paste0(figpath, "p3_wnncluster_dsb.pdf"), width = 3.5, height = 3.5)

# set theme for biaxial plots
mg_theme = list( 
  theme_bw(),
  theme(axis.title.x =element_text(size = 18), axis.title.y = element_text(size = 18)), 
  geom_bin2d(bins = 200, show.legend = FALSE),
  viridis::scale_fill_viridis(option = "B"), 
  geom_vline(xintercept = 0, linetype = 'dashed'),
  geom_hline(yintercept = 0, linetype = 'dashed')
)

# dsb 
d1 = cbind(s@meta.data, as.data.frame(t(s@assays$CITE@data)))
p = ggplot(d1, aes(x = CD19, y = CD3)) +
  mg_theme + 
  geom_hline(yintercept = 3.5, color = 'red') + 
  geom_vline(xintercept = 3.5, color = 'red')
ggsave(p, filename = paste0(figpath, 'dsb_cd19_cd3.pdf'), width = 3, height = 3)
p = ggplot(d1, aes(x = CD4, y = CD14)) + 
  mg_theme + 
  geom_hline(yintercept = 3.5, color = 'red') + 
  geom_vline(xintercept = 3.5, color = 'red')
ggsave(p, filename = paste0(figpath, 'dsb_cd4_cd14.pdf'), width = 3, height = 3)

# clr 
d2 = cbind(s@meta.data, as.data.frame(t(s@assays$CLR@data)))
p = ggplot(d2, aes(x = CD19, y = CD3)) + mg_theme
ggsave(p, filename = paste0(figpath, 'clr_cd19_cd3.pdf'), width = 3, height = 3)
p = ggplot(d2, aes(x = CD4, y = CD14)) + mg_theme
ggsave(p, filename = paste0(figpath, 'clr_cd4_cd14.pdf'), width = 3, height = 3)

# library size normalization
ln = CreateAssayObject(counts = s@assays$CITE@counts)
s[["lognorm"]] = ln
s = NormalizeData(s, assay = "lognorm",normalization.method = "LogNormalize")
d3 = cbind(s@meta.data, as.data.frame(t(as.matrix(s@assays$lognorm@data))))
p = ggplot(d3, aes(x = CD19, y = CD3)) + mg_theme
ggsave(p, filename = paste0(figpath, 'Lognorm_cd19_cd3.pdf'), width = 3, height = 3)
p = ggplot(d3, aes(x = CD4, y = CD14)) + mg_theme
ggsave(p, filename = paste0(figpath, 'Lognorm_cd4_cd14.pdf'), width = 3, height = 3)


## dsb clusters dsb vs clr values 
d = cbind(as.data.frame(t(s@assays$CITE@data)), s@meta.data) %>% 
  group_by(dsb_pseudo_knn_res.2) %>% 
  gather(prot, count, CD10:`TCR-g/d`) %>% 
  group_by(prot, dsb_pseudo_knn_res.2) %>% 
  summarize(median_dsb = median(count), 
  )

clr_df = cbind(as.data.frame(t(s@assays$CLR@data)), s@meta.data) %>% 
  group_by(dsb_pseudo_knn_res.2) %>% 
  gather(prot, count, CD10:`TCR-g/d`) %>% 
  group_by(prot, dsb_pseudo_knn_res.2) %>% 
  summarize(median_clr = median(count), 
  )

# add neg median for each protein  to dsb_df
q98 = function(x){quantile(x, probs = 0.98)}
neg_med = data.frame(median_neg = apply(log10(adt_background + 1), 1, median) ) %>% rownames_to_column("prot")
neg_q98 = data.frame(q98_neg = apply(log10(adt_background + 1), 1, q98)) %>% rownames_to_column("prot")

# merge summary stats 
d = full_join(d, neg_med, by = "prot")
d = full_join(d, neg_q98, by = "prot")
d$median_clr = clr_df$median_clr

# plot dsb vs clr 
p = ggplot(d, aes(x = median_clr, y = median_dsb, label = prot)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ dsb_pseudo_knn_res.2, nrow = 3) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  geom_vline(xintercept = 1, color = "red") + 
  xlab(" median clr normalized across cells within cluster") + 
  ylab(" dsb normalized median within cluster ") + 
  ggrepel::geom_text_repel(data = d %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 2.4, force = 2) + 
  ggrepel::geom_text_repel(data = d %>% filter(median_clr > 1.5 & median_dsb < 3.5), 
                           segment.size = 0.5, size = 2.4, color = 'red', force = 2) + 
  ggtitle('dsb WNN clusters')
ggsave(plot = p, filename = paste0(figpath, 'dsb_clr_dsb_clusters.pdf'), width = 12, height = 6, device = "pdf")

# cluster 14 highlight 
c14 = d %>% filter(dsb_pseudo_knn_res.2 %in% c('14'))
p = ggplot(c14, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point() + 
  geom_point(data = c14 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = 'deepskyblue3')+
  theme_bw() + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  theme(strip.text = element_text(size = 12 )) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab("empty drop median log10+1") + 
  ylab("dsb normalized median") + 
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 3, force = 2) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_clr > 1 & median_dsb < 3.5), 
                           nudge_x = 0.4, segment.color = 'grey', segment.size = 0.5, size = 2, color = 'red', force = 6) + 
  ggtitle('dsb WNN clusters')
p
ggsave(p, filename = paste0(figpath, 'c14_dsb_emptymedian.pdf'), width = 3.5, height = 3.5)

# dsb vs empty drop q98
p = ggplot(c14, aes(x = q98_neg, y = median_dsb, label = prot)) + 
  geom_point() + 
  geom_point(data = c14 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = 'deepskyblue3')+
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', label.x.npc = 0.6, label.y.npc = 0.85) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  theme(strip.text = element_text(size = 12 )) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab("background drop 98th percentile") + 
  ylab("dsb normalized median") + 
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 4, force = 2) 
  # ggrepel::geom_text_repel(data = c14 %>% filter(median_clr > 1 & median_dsb < 3.5),
  #                          nudge_x = 0.4, segment.color = 'grey', segment.size = 0.5, size = 2, color = 'red', force = 6) + 
ggsave(p, filename = paste0(figpath, 'c14_dsb_empty_Q98.pdf'), width = 3.5, height = 3.5)

# clr vs empty droplets 
p = ggplot(c14, aes(x = q98_neg, y = median_clr, label = prot)) + 
  geom_point() + 
  geom_point(data = c14 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = "grey60")+
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', label.x.npc = 0.5, label.y.npc = 0.8) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  theme(strip.text = element_text(size = 12 )) + 
  geom_hline(yintercept = 1, color = "red") + 
  xlab("background drop 98th percentile log10+1") + 
  ylab("CLR transformed median ") + 
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 4, force = 2, seed = 1990) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_clr > 1 & median_dsb < 3.5),  seed = 1990,
                           nudge_x = 0.2, segment.color = 'grey', segment.size = 0.5, size = 2, color = 'red', force = 6) 
ggsave(p, filename = paste0(figpath, 'c14_q89_clr.pdf'), width = 3.5, height = 3.5)

# CLR vs dsb 
p = ggplot(c14, aes(x = median_clr, y = median_dsb, label = prot)) + 
  geom_point() + 
  geom_point(data = c14 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = 'deepskyblue3')+
  theme_bw() + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  theme(strip.text = element_text(size = 12 )) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  geom_vline(xintercept = 1, color = "red") + 
  xlab("CLR transformed median ") + 
  ylab("dsb normalized median") + 
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 3, force = 2) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_clr > 1 & median_dsb < 3.5), 
                           nudge_x = 0.1, segment.color = 'grey', segment.size = 0.5, size = 2, color = 'red', force = 6) 
ggsave(p, filename = paste0(figpath, 'c14_dsb_clr.pdf'), width = 3.5, height = 3.5)


# de genes c14 vs others
Idents(s) = 'dsb_pseudo_knn_res.2'
DefaultAssay(s) = 'RNA'
c14_de = FindMarkers(object = s, ident.1 = '14', test.use = 'roc')
c14_de = c14_de %>% rownames_to_column('gene') 

p = ggplot(c14_de %>% filter(avg_log2FC > 0), aes(x = avg_log2FC, y = myAUC, label = gene)) + 
  theme_bw() + 
  ggrepel::geom_text_repel(data = c14_de %>% 
                             filter(myAUC > 0.75 ), segment.size = 0.3, size = 3, force = 2) +
  geom_point() + 
  xlab("average log2 fold change ") +
  ylab('AUC Cluster 14 classifier')
p
ggsave(p, filename = paste0(figpath, 'c14_degenes_full.pdf'), width = 3.5, height = 3.5)

# highlight fold change contingency in park et al  genes 
pg = data.table::fread(file = here('V2/teaseq/park_scirep_maitgene/PMID31068647_table1gene.txt'))
pg = pg %>% filter(gene %in% c14_de$gene)
c14_de_sub = c14_de %>% filter(gene %in% pg$gene)
dg = full_join(c14_de_sub, pg, by = 'gene')

# plot 
p = ggplot(dg, aes(x = avg_log2FC, y = myAUC, label = gene)) + 
  theme_bw() + 
  geom_point(shape = 16, size = 3) + 
  xlab("c14 vs other clusters log2FC") +
  ylab("Sorted MAIT vs. TCR-Va7.2 negative T cells") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  ggrepel::geom_text_repel(data = dg, segment.size = 0.5, size = 5, force = 2) 
ggsave(p, filename = paste0(figpath, 'c14_degenes_park.pdf'), width = 3.5, height = 3.5)
p



ggplot(dsb_rna_m %>% filter(cluster == '3' & avg_log2FC > 0),
       aes(x = avg_log2FC, y = myAUC, label = gene)) +
  theme_bw() + 
  ggrepel::geom_text_repel(data = dsb_rna_m %>% filter( cluster == '3' & myAUC > 0.9),
                           segment.size = 0.5, size = 3, force = 2, max.overlaps = 20) + 
  geom_point() + 
  xlab('average mRNA log2 fold change') + 
  ylab('AUC Cluster 3 classifier')

# CLR vs dsb cluster contingency table 
pheatmap::pheatmap(log(as.matrix(table(s@meta.data$clr_pseudo_knn_res.2, s@meta.data$dsb_pseudo_knn_res.2))+1), 
                   color = viridis::inferno(n = 8),
                   fontsize_row = 10, fontsize_col = 10, main = 'clr (y) dsb (x)',
                   treeheight_row = 0, treeheight_col = 0,width = 5, height = 5,
                   filename = paste0(figpath, 'dsb_clr_contingency.pdf'))

# dsb heatmap 
prots = rownames(s@assays$CITE@data)
dh = cbind(as.data.frame(t(s@assays$CITE@data)), s@meta.data) %>% 
  group_by(dsb_pseudo_knn_res.2) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  remove_rownames() %>% 
  column_to_rownames('dsb_pseudo_knn_res.2') %>% 
  as.matrix()
mat1 = t(dh)
rownames(mat1)[rownames(mat1) == 'IgG1-K-Isotype-Control']  = 'IgG1-K-Isotype'

ph = pheatmap::pheatmap(mat1, treeheight_row = 10, treeheight_col = 10,
                        color = viridis::inferno(10),
                        breaks = seq(from = 0, to = 20, length.out = 11),
                        width = 4.5, height = 5.8, 
                        border_color = NA, filename = paste0(figpath, 'dsb_heatmap.pdf'))

# extract dsb orders 
porder = ph$tree_row$labels[ph$tree_row$order]

# clr heatmap 
ch = cbind(as.data.frame(t(s@assays$CLR@data)), s@meta.data) %>% 
  group_by(clr_pseudo_knn_res.2) %>% 
  summarize_at(.vars = prots,.funs = mean) %>% 
  remove_rownames() %>% 
  column_to_rownames('clr_pseudo_knn_res.2') %>% 
  as.matrix()
mat = t(ch)
rownames(mat)[rownames(mat) == 'IgG1-K-Isotype-Control']  = 'IgG1-K-Isotype'
pheatmap::pheatmap(mat[porder, ],
                   cluster_rows = FALSE,
                   treeheight_row = 10, treeheight_col = 10,fontsize_col = 8,
                   color = viridis::inferno(n = 10),
                   width = 4.5, height = 5.8,
                   border_color = NA, filename = paste0(figpath, 'clr_heatmap.pdf'))

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
#   [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
# [8] tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       cowplot_1.1.1      magrittr_2.0.1    
# [15] ggridges_0.5.3     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        rprojroot_2.0.2       mclust_5.4.7         
# [7] fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0   farver_2.0.3          leiden_0.3.7          listenv_0.8.0        
# [13] bit64_4.0.5           ggrepel_0.9.1         fansi_0.4.2           RSpectra_0.16-0       lubridate_1.7.9.2     xml2_1.3.2           
# [19] codetools_0.2-18      splines_4.0.5         polyclip_1.10-0       jsonlite_1.7.2        broom_0.7.5           ica_1.0-2            
# [25] cluster_2.1.1         dbplyr_2.1.0          png_0.1-7             pheatmap_1.0.12       uwot_0.1.10           shiny_1.6.0          
# [31] sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5        httr_1.4.2            backports_1.2.1       assertthat_0.2.1     
# [37] Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2        cli_2.3.0             limma_3.46.0          later_1.1.0.1        
# [43] htmltools_0.5.1.1     tools_4.0.5           igraph_1.2.6          gtable_0.3.0          glue_1.4.2            RANN_2.6.1           
# [49] reshape2_1.4.4        Rcpp_1.0.6            scattermore_0.7       cellranger_1.1.0      vctrs_0.3.6           nlme_3.1-152         
# [55] lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10             miniUI_0.1.1.1        lifecycle_1.0.0      
# [61] irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1         zoo_1.8-8             scales_1.1.1         
# [67] spatstat.core_2.0-0   hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0  parallel_4.0.5        RColorBrewer_1.1-2   
# [73] reticulate_1.18       pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15          stringi_1.5.3         rlang_0.4.10         
# [79] pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41       ROCR_1.0-11           tensor_1.5            labeling_0.4.2       
# [85] patchwork_1.1.1       htmlwidgets_1.5.3     bit_4.0.4             tidyselect_1.1.0      parallelly_1.23.0     ggsci_2.9            
# [91] RcppAnnoy_0.0.18      plyr_1.8.6            R6_2.5.0              generics_0.1.0        DBI_1.1.1             withr_2.4.1          
# [97] pillar_1.4.7          haven_2.3.1           mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10       abind_1.4-5          
# [103] future.apply_1.7.0    hdf5r_1.3.3           modelr_0.1.8          crayon_1.4.1          utf8_1.1.4            KernSmooth_2.23-18   
# [109] spatstat.geom_2.0-1   plotly_4.9.3          viridis_0.5.1         grid_4.0.5            readxl_1.3.1          data.table_1.14.0    
# [115] reprex_1.0.0          digest_0.6.27         xtable_1.8-4          httpuv_1.5.5          munsell_0.5.0         viridisLite_0.3.0   