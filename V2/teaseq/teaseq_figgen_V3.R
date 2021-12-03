# R4 Seurat 4 
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

# cluster contingency between the joint clustering models with different 
# input adt data (dsb vs CLR normalized)
cs = chisq.test(table(s@meta.data$dsb_pseudo_knn_res.2,s@meta.data$clr_pseudo_knn_res.2))
# Pearson's Chi-squared test
# data:  table(s@meta.data$dsb_pseudo_knn_res.2, s@meta.data$clr_pseudo_knn_res.2)
# X-squared = 507683, df = 575, p-value < 2.2e-16
pheatmap::pheatmap(table(s@meta.data$dsb_pseudo_knn_res.2,s@meta.data$clr_pseudo_knn_res.2),
                   treeheight_col = 15, treeheight_row = 15,
                   fontsize_col = 5, fontsize_row = 5, 
                   filename = paste0(figpath,'heat_clustercontingency.pdf'),
                   width = 2.6, height = 2.6) 
# lower contingency area plot 
library(gplots)
pdf(file =  paste0(figpath,'clustercontingency.pdf'),  width = 4.5, height = 4)
balloonplot(table(s@meta.data$clr_pseudo_knn_res.2,s@meta.data$dsb_pseudo_knn_res.2),
            main = '',
            dotcolor = 'black',
            xlab ="", 
            ylab="", rowmar = 10, colmar = 10, text.size = 0.5,
            label = FALSE, 
            show.margins = FALSE)
dev.off()

# set theme 
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
col_use = BuenColors::jdb_palette("corona")
length(col_use)
col_use = col_use[1:24]
col_use[15] = '#6978fd'
col_use[24] = 'orange'
DefaultAssay(s) <- 'CITE'
p = AugmentPlot(DimPlot(s, reduction = 'dsb_wnn_umap', group.by = 'dsb_pseudo_knn_res.2',
                         label = TRUE, repel = TRUE, label.size = 9.5, pt.size = 1.2, cols = col_use)) + 
  boxbox + 
  NoLegend() + 
  xlab('protein + RNA WNN UMAP 1') + 
  ylab('protein + RNA WNN UMAP 2') +
  scale_color_manual(values = BuenColors::jdb_palette("corona")) + 
  ggtitle("dsb normalized TEA-seq protein data")
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

p = ggplot(d1 %>% filter(CD14 > -5), aes(x = CD4, y = CD14)) + 
  mg_theme + 
  geom_hline(yintercept = 3.5, color = 'red') + 
  geom_vline(xintercept = 3.5, color = 'red') 
ggsave(p, filename = paste0(figpath, 'dsb_cd4_cd14.pdf'), width = 3.5, height = 3.5)

# clr 
d2 = cbind(s@meta.data, as.data.frame(t(s@assays$CLR@data)))
p = ggplot(d2, aes(x = CD4, y = CD14)) + mg_theme
ggsave(p, filename = paste0(figpath, 'clr_cd4_cd14.pdf'), width = 3, height = 3)

# library size normalization
ln = CreateAssayObject(counts = s@assays$CITE@counts)
s[["lognorm"]] = ln
s = NormalizeData(s, assay = "lognorm",normalization.method = "LogNormalize")
d3 = cbind(s@meta.data, as.data.frame(t(as.matrix(s@assays$lognorm@data))))
p = ggplot(d3, aes(x = CD4, y = CD14)) + mg_theme
ggsave(p, filename = paste0(figpath, 'Lognorm_cd4_cd14.pdf'), width = 3, height = 3)


# visualization of cluster 14 TCR proteins 
# manual gate theme 
plottheme = list(theme_bw(),
                 geom_vline(xintercept = 0),
                 geom_hline(yintercept = 0)
                 )

p1 = 
  ggplot(d1 %>% filter(dsb_pseudo_knn_res.2 == '14'), aes(x = `TCR-a/b`,y = `TCR-Va7.2`)) + 
  plottheme + 
  geom_vline(xintercept = 3.5, color = 'red') + 
  geom_hline(yintercept = 3.5, color = 'red') +
  geom_density_2d() +
  xlim(c(-5,33)) + ylim(-5,33) + 
  theme(aspect.ratio = 1) 
  
p2 = 
  ggplot(d2 %>% filter(dsb_pseudo_knn_res.2 == '14'), aes(x = `TCR-a/b`,y = `TCR-Va7.2`)) + 
  plottheme + 
  geom_density_2d() +
  xlim(c(-1,3)) + ylim(-1,3) + 
  theme(aspect.ratio = 1) 
  
p3 = p1 | p2
ggsave(p3,filename = paste0(figpath ,'C14_TCRs.pdf'), width = 5.5, height = 3)


# combined plots 
# comparisons of values normalized by dsb vs clr by cluster 
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
d = full_join(d, clr_df)

# cluster 14 highlight 
# dsb vs empty drop q98
c14 = d %>% filter(dsb_pseudo_knn_res.2 == '14')
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
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 5, force = 2) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(prot %in% c("TCR-a/b","TCR-g/d")),  seed = 1990, box.padding = 0.5,
                         nudge_x = 0.2, nudge_y = -0.2, segment.color = 'grey', segment.size = 0.5, size = 5, color = 'red', force = 6) 
ggsave(p, filename = paste0(figpath, 'c14_dsb_empty_Q98.pdf'), width = 3.5, height = 3.5)

# clr c14 from dsb model vs empty droplets 
p = ggplot(c14, aes(x = q98_neg, y = median_clr, label = prot)) + 
  geom_point() + 
  geom_point(data = c14 %>% filter(median_dsb > 3.5), shape = 21, size = 2.5, fill = "grey60")+
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', label.x.npc = 0.5, label.y.npc = 0.8) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  theme(strip.text = element_text(size = 12 )) + 
  xlab("background drop 98th percentile log10+1") + 
  ylab("CLR transformed median ") + 
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(median_dsb > 3.5), segment.size = 0.5, size = 5, force = 2, seed = 1990, nudge_x = -0.2, nudge_y = -0.1) + 
  ggrepel::geom_text_repel(data = c14 %>% filter(prot %in% c("TCR-a/b","TCR-g/d")),  seed = 1990, box.padding = 0.7,
                           nudge_x = 0.2, nudge_y = 0.2, segment.color = 'grey', segment.size = 0.5, size = 5, color = 'red', force = 6) 
ggsave(p, filename = paste0(figpath, 'c14_q89_clr.pdf'), width = 3.5, height = 3.5)

###### differential expression of mait cluster
# identify genes for cluster 14 vs other clusters 
Idents(s) = 'dsb_pseudo_knn_res.2'
DefaultAssay(s) = 'RNA'
c14_de = FindMarkers(object = s, ident.1 = '14', test.use = 'roc',
                     logfc.threshold = 0,
                     min.pct = 0,
                     min.cells.feature = 0)
c14_de = c14_de %>% rownames_to_column('gene') 
saveRDS(c14_de, file = paste0(datapath,'c14_de.rds'))

# (loads saved differential expression results)
# c14_de = readRDS(file = here('V2/teaseq/generated_data_v2/TEA-seq/c14_de.rds'))
p = ggplot(c14_de %>% filter(myAUC > 0.6 & avg_log2FC > 0), aes(x = avg_log2FC, y = myAUC, label = gene)) + 
  theme_bw() + 
  ggrepel::geom_text_repel(data = c14_de %>% 
                             filter(myAUC > 0.75 ), segment.size = 0.3, size = 3.3, force = 2, seed = 1990) +
  geom_point() + 
  geom_bin2d(data = c14_de %>% filter(myAUC<0.6 & avg_log2FC > 0), bins = 200, fill = 'black') + 
  xlab("average log2 fold change ") +
  ylab('AUC Cluster 14 classifier')
p
ggsave(p, filename = paste0(figpath, 'c14_degenes_full.pdf'), width = 3.5, height = 3.5)

# Gene Set Enrichment of mait sig vs top Cluster 14 differentially expressed genes
# load mait signature formed from PMID: 31068647 Supplementary Figure 1a  
# https://www.nature.com/articles/s41598-019-43578-9
# genes extracted from pdf supplementary fig 1a
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-43578-9/MediaObjects/41598_2019_43578_MOESM1_ESM.pdf
# format mait signature 
msig = suppressWarnings(read.csv(file = here('data/revision_data/tea_seq/maitsig.csv'))) %>% 
  colnames() %>% 
  as.character()
maitsig = list(msig)
names(maitsig) = 'maitsig'
maitsig

# format Cluster 14 genes for fgsea;
# rank genes by by log fold change in C14 vs other clusters 
derank = c14_de %>% arrange(desc(avg_log2FC))
gene_ranks = derank$avg_log2FC
names(gene_ranks) = derank$gene


# run fast gene set enrichment 
library(fgsea)
gs = fgsea(maitsig, gene_ranks)
# print results 
gs
# pathway  pval  padj log2err        ES      NES size                                    leadingEdge
# 1: maitsig 1e-10 1e-10      NA 0.9023002 2.631651  154 SLC7A5,PHACTR2,ADAM12,SLC4A10,ZBTB16,IL4I1,...

# plot enrichment 
enrline = list(geom_line(color = "#6978fd", size = 2 )) 
p = plotEnrichment(pathway = maitsig[[1]], stats = gene_ranks, gseaParam = 1, ticksSize = 0.1)
p = p + enrline + ggtitle(paste0('Flow-sorted TCRVa7.2+ cell signature \nCluster 14 vs other cells', '\ngene set enrichment padj = ' ,gs$padj)) 
ggsave(p, filename = paste0(figpath,'maitsig_c14.pdf'), width = 4, height = 3.8)

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
#   [1] fgsea_1.16.0       gplots_3.1.1       here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0       
# [9] tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    magrittr_2.0.1     SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       fastmatch_1.1-0       BuenColors_0.5.6      plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2       
# [8] splines_4.0.5         BiocParallel_1.24.1   listenv_0.8.0         scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1     viridis_0.5.1        
# [15] fansi_0.4.2           tensor_1.5            cluster_2.1.2         ROCR_1.0-11           openxlsx_4.2.3        globals_0.14.0        modelr_0.1.8         
# [22] matrixStats_0.58.0    spatstat.sparse_2.0-0 colorspace_2.0-0      rvest_0.3.6           ggrepel_0.9.1         haven_2.3.1           crayon_1.4.1         
# [29] jsonlite_1.7.2        spatstat.data_2.1-0   survival_3.2-10       zoo_1.8-8             glue_1.4.2            polyclip_1.10-0       gtable_0.3.0         
# [36] leiden_0.3.7          car_3.0-10            future.apply_1.7.0    abind_1.4-5           scales_1.1.1          pheatmap_1.0.12       DBI_1.1.1            
# [43] rstatix_0.7.0         miniUI_0.1.1.1        Rcpp_1.0.6            isoband_0.2.3         viridisLite_0.3.0     xtable_1.8-4          reticulate_1.18      
# [50] spatstat.core_2.0-0   foreign_0.8-81        htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2            
# [57] pkgconfig_2.0.3       farver_2.0.3          uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10         utf8_1.1.4            tidyselect_1.1.0     
# [64] labeling_0.4.2        rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0         cellranger_1.1.0      tools_4.0.5          
# [71] cli_2.5.0             generics_0.1.0        broom_0.7.5           ggridges_0.5.3        fastmap_1.1.0         goftest_1.2-2         fs_1.5.0             
# [78] fitdistrplus_1.1-3    zip_2.1.1             caTools_1.18.1        RANN_2.6.1            pbapply_1.4-3         future_1.21.0         nlme_3.1-152         
# [85] mime_0.10             xml2_1.3.2            compiler_4.0.5        rstudioapi_0.13       plotly_4.9.3          curl_4.3              png_0.1-7            
# [92] ggsignif_0.6.0        spatstat.utils_2.1-0  reprex_1.0.0          stringi_1.5.3         lattice_0.20-41       Matrix_1.3-2          vctrs_0.3.6          
# [99] pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18      data.table_1.14.0     cowplot_1.1.1        
# [106] bitops_1.0-6          irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       R6_2.5.0              promises_1.2.0.1      KernSmooth_2.23-18   
# [113] gridExtra_2.3         rio_0.5.16            parallelly_1.23.0     codetools_0.2-18      MASS_7.3-53.1         gtools_3.8.2          assertthat_0.2.1     
# [120] rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           parallel_4.0.5        hms_1.0.0             grid_4.0.5           
# [127] rpart_4.1-15          carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0          shiny_1.6.0           lubridate_1.7.9.2    

