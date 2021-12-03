set.seed(1)
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# read combined dataframe for visualization
df = readRDS(here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_sng_metadata_umapdim_prots_dataframe.rds"))
cu = pals::kelly(n = 22); cu[1] = "dodgerblue" 

# plot umap 
centers = df %>% dplyr::group_by(p3_dist_3) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p = ggplot(df, aes(x = -1*UMAP1, y = UMAP2)) + 
  theme_void() + 
  geom_point(mapping = aes(fill = p3_dist_3), size = 1, shape = 21, stroke = 0, alpha = 0.6, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data = centers, size = 8.5, mapping = aes(label = p3_dist_3), show.legend = FALSE) +
  scale_fill_manual(values = cu) + 
  theme(legend.title =  element_blank()) + 
  theme(legend.text = element_text(colour="black", size=8, face="bold")) 
ggsave(p, filename = paste0(figpath,"h1_umap_DSB.png"),width = 8, height = 7)

# plot distribution of lineage proteins on umap
prot_plot = c("CD3_PROT", "CD4_PROT", "CD8_PROT", "CD56_PROT", "CD16_PROT",
              "CD14_PROT", "CD19_PROT", "IgD_PROT", "CD303_PROT")
index1 = prot_plot[1]
index2 = prot_plot[length(prot_plot)]
df2 = df %>% select(UMAP1, UMAP2, prot_plot) %>% gather(prot, dsb, index1:index2) %>% filter(dsb > -10)
df2$prot = factor(df2$prot, levels = prot_plot)
p = ggplot(df2 %>% filter(prot %in% prot_plot)) + aes(x = -1*UMAP1,y =UMAP2, color = dsb) + 
  geom_point(show.legend = TRUE, size = 0.2, shape = 16, alpha = 0.7) + 
  scale_color_viridis_c(option = "B", limits = c(-5,30)) + 
  theme_void() + 
  theme(strip.background = element_blank()) + 
  theme(legend.key.size  = unit(1, units = "cm")) + 
  theme(strip.text = element_text(color = "black", size = 18, face = "bold"), 
        legend.text = element_text(color = "black", size = 18, face = "bold"),
        legend.title = element_text(color = "black", size = 18, face = "bold")) + 
  facet_wrap(~prot, nrow = 3)
ggsave(p, filename = paste0(figpath,"all_lineage.png"), height = 9, width = 12)


## Figure 2 comparison of NK cell noise vs true staining proteins based on DSB in normalized data and neg drops. 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
raw = h1@assay$CITE@raw.data
dsb = h1@assay$CITE@data
md = h1@meta.data

# get NK cell barcode vectors 
nk_cells = md %>% filter(celltype_label_3 == "CD16++ NK") %$% barcode_check

# subset raw NK data 
raw = as.data.frame(as.matrix(t(raw[ , nk_cells ])))
dsb = as.data.frame(as.matrix(t(dsb[ , nk_cells ])))

# load negative drops 
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))
neg_adt = do.call(cbind, neg_adt)
neg = as.data.frame(as.matrix(t(neg_adt)))

# add clr values 
# import CLR across cells 
cell_clr = readRDS(file = here('V2/dsb_normalize_cluster_pipeline/generated_data/h1_CLR_acroass_cells_matrix.rds'))
nkc = as.data.frame(as.matrix(t(cell_clr[ , nk_cells ])))

# calc summary stats for highlight cluster and background drops 
q98 = function(x){quantile(x, probs = 0.98)}
df1 = data.frame(median_c4_raw_log10 = apply(log10(raw + 1), 2, median), 
                 median_neg_log10 = apply(log10(neg + 1), 2, median),
                 q98_neg_log10 = apply(log10(neg + 1), 2, q98),
                 median_c4_dsb = apply(dsb, 2, median), 
                 median_c4_clr = apply(nkc, 2, median)) %>%
  as.data.frame() %>% 
  rownames_to_column("protein") %>% 
  mutate(protein = str_remove(protein, pattern = "_PROT")) 

# add variable for subset of proteins to label
df1 = df1 %>% mutate(noise_ = if_else(median_c4_dsb < 3.5 & median_c4_raw_log10 > 1, true = '1',false = '0'))

# log prot 
p1 = ggplot(df1, aes( x = median_neg_log10, y = median_c4_raw_log10, label = protein, color = noise_)) + 
  geom_point(shape = 16 , alpha  = 0.7, show.legend = FALSE) + 
  geom_abline() + 
  theme_bw() + 
  theme(axis.title.y = element_text(size = 14)) + 
  scale_color_manual(values = c("black", "red")) +
  labs(x = "background drops median log10+1",  y = "cluster 4 median log10+1") + 
  xlab("") +
  ggrepel::geom_text_repel(data = df1 %>% filter(median_c4_raw_log10 > 1), 
                           size = 2.7, segment.color = "grey",   show.legend = FALSE) 
ggsave(p1,filename = paste0(figpath,"nk_log_vs_empty_2.pdf"), width = 2.2, height = 3)

# dsb
p2 = ggplot(df1, aes(x = median_neg_log10, y = median_c4_dsb, label = protein, color = noise_)) + 
  geom_point(shape = 16 , alpha  = 0.7, show.legend = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = 14)) + 
  labs(x = "background drops median log10 +1", y = "cluster 4 median dsb") +
  xlab("") +
  geom_hline(yintercept = 3.5, color = "red3") + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df1 %>% filter(median_c4_raw_log10 > 1),
                           segment.color = "grey", 
                           size = 2.7, nudge_y = -0.5, show.legend = FALSE) 
ggsave(p2,filename = paste0(figpath,"nk_dsb_vs_empty_2.pdf"), width = 2.2, height = 3)

pn = ggplot(df1, aes(x = median_neg_log10, y = median_c4_clr, label = protein, color = noise_)) + 
  geom_point(shape = 16 , alpha  = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + 
  geom_abline() + 
  theme(axis.title.y = element_text(size = 14)) + 
  labs(x = "background drops median log10 +1", y = "cluster 4 median CLR") +
  xlab("") + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df1 %>% filter(median_c4_raw_log10 > 1),
                           segment.color = "grey", 
                           size = 2.7, nudge_y = 0.1, nudge_x = -0.1, show.legend = FALSE) 
ggsave(pn,filename = paste0(figpath,"nk_clr_vs_empty_2.pdf"), width = 2.2, height = 3)



# tidy clr value across cells for NK cells 
nk_clr = cell_clr[ ,nk_cells]
nk_clr = nk_clr %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('bc') %>% 
  gather(protein, CLR, `AnnexinV-PROT`:`CD20-PROT`) %>% 
  mutate(protein = str_remove(protein, pattern = "-PROT")) 

# tidy dsb nk cells -- data 'dsb' above is subset to only include the nkcells 
nk_dsb = dsb %>% 
  rownames_to_column('bc') %>% 
  gather(protein, dsb, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(protein = str_remove(protein, pattern = "_PROT")) 

# full join dsb and clr 
d = full_join(nk_clr, nk_dsb)
d = d %>% gather(norm_method, value, dsb:CLR)

# calc mean clr for plot threshold 
psub = d %>% filter(norm_method == "CLR") %>% 
  group_by(protein) %>% 
  summarize(mean_clr = mean(value)) %>% 
  filter(mean_clr > 0.55) %$% 
  protein
d$norm_method[d$norm_method == "CLR"] = "CLR (across cells)"
p = ggplot(d %>% filter(value > -5 & value < 30 & protein %in% psub), 
           aes(x = value, y = reorder(protein,value) , fill = norm_method, color = norm_method)) + 
  theme_bw() + 
  facet_wrap(~norm_method, scales = "free_x") + 
  ggridges::geom_density_ridges2(show.legend = FALSE, alpha = 0.8, size = 0.3) + 
  scale_fill_manual(values = c("grey60", "deepskyblue3")) + 
  scale_color_manual(values = c("grey60", "deepskyblue3")) + 
  theme(strip.background = element_blank()) + 
  theme(axis.text.y =  element_text(size = 10, color = "black")) + 
  xlab("") + ylab("") + 
  geom_vline(xintercept = 0,  size = 0.7, linetype = 'dashed') 
ggsave(p, filename = paste0(figpath, "2_dsb_clr_nkcell_distribution_topmedian.pdf"), width = 4.4, height = 5)



# protein annotaiton plots vs empty drops
# version 3 with smaller main figure and rest in supplement
dsb_df = cbind(as.data.frame(t(h1@assay$CITE@data)), md) %>% 
  group_by(celltype_label_3, p3_dist_3) %>% 
  gather(prot, val, AnnexinV_PROT:CD20_PROT) %>% 
  group_by(prot, p3_dist_3, celltype_label_3) %>% 
  summarize(median_dsb = median(val)) %>% 
  ungroup() 

# add neg median for each protein  to dsb_df
dsb_df$cluster_celltype = paste(dsb_df$p3_dist_3, dsb_df$celltype_label_3, sep = ": ")
neg_med = data.frame(median_neg = apply(log10(neg + 1), 2, median)) %>% rownames_to_column("prot")
dsb_df = full_join(dsb_df, neg_med, by = "prot") %>% 
  mutate(prot = str_remove(prot, pattern = "_PROT"))

#### supplemental figure 
marker_highlight = c("CD1d", "CD1c", "CD14","CD103", "CD16", "CD3", "CD4", "CD8", "CD28", "CD161",
                     "CD45RO", "CD45RA", "CD33", "CD56", "CD71", "CD27", "CD244", "KLRG1",
                     "CD195", "CD38", "CD127",  "CD16", "CD34")

# remove the example B cell and dcs that are shown in the main figure (top row from v28)
highl = c("13: Unswitched B cells", "5: Transitional B cells", "12: Switched B cells",  "16: pDC", "14: mDC")
subs = dsb_df %>% filter(!cluster_celltype %in% highl)
p = ggplot(subs, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point(shape = 16 , alpha  = 0.5) + 
  geom_point(data = subs %>% filter(prot %in% marker_highlight), shape = 21, size = 2.5, fill = 'deepskyblue3')+
  ggrepel::geom_text_repel(data = subs %>% filter(prot %in% marker_highlight & median_dsb > 3.5),
                           segment.size = 0.5, segment.color = "grey", size = 4.7, nudge_y = 2, box.padding = 0.5, force = 5) +
  theme_bw() +
  facet_wrap(~ cluster_celltype, nrow = 3) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 16)) + 
  theme(strip.text = element_text(size = 10)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab(" median log + 1 protein in empty droplets") + 
  ylab(" median dsb normalized  protein ") 
ggsave(p, filename = paste0(figpath,"V2supplneg_vs_dsb_highlight_2.pdf"), width = 18, height = 8)


# smaller version for main figure with the B and DC subsets highlighted
cts = c("13: Unswitched B cells", "5: Transitional B cells", "12: Switched B cells")
mh = c("CD20", "CD19", "IgD", "IgM")
for (i in 1:length(cts)) {
  
  psub =  dsb_df %>% filter(cluster_celltype == cts[i])
  p = ggplot(psub, 
             aes(x = median_neg, y = median_dsb, label = prot)) + 
    geom_point(shape = 16 , alpha  = 0.5) + 
    geom_point(data = psub %>% filter(prot %in% mh), shape = 21, size = 2.5, fill = 'deepskyblue3') +
    ggrepel::geom_text_repel(data = psub %>% filter(prot %in% mh & median_dsb > 3.5),
                             segment.size = 0.5, segment.color = "grey", size = 5, force = 10, box.padding = 1, nudge_y = 1) +
    ggrepel::geom_text_repel(data = psub %>% filter(prot %in% mh & median_dsb < 3.5),
                             segment.size = 0.5, segment.color = "grey", size = 3, force = 3, color = 'red', nudge_y = -0.4) +
    theme_bw() + 
    ylim(c(-2,22)) + 
    ggtitle(cts[i]) + 
    #geom_segment(data = psub %>% filter(prot %in% mh ), aes(yend=median_dsb, xend=median_neg)) +
    theme(strip.background = element_blank(), axis.title = element_text(size = 0)) + 
    geom_hline(yintercept = 3.5, color = "red") + 
    xlab(" median log10 + 1 protein in empty droplets ") + 
    ylab(" median dsb normalized  protein ") 
  p
  ggsave(p, filename = paste0(figpath,cts[i],".pdf"), width = 3, height = 3)
}

cts = c("10: Non-classical Monocytes", "2: Classical Monocytes","16: pDC", "14: mDC", "21: HSC")
mh = c("CD34", "CD14", "CD16", "CD303",  "CD1c", "CD1d")
for (i in 1:length(cts)) {
  
  psub =  dsb_df %>% filter(cluster_celltype == cts[i])
  p = ggplot(psub, 
             aes(x = median_neg, y = median_dsb, label = prot)) + 
    geom_point(shape = 16 , alpha  = 0.5) + 
    geom_point(data = psub %>% filter(prot %in% mh), shape = 21, size = 2.5, fill = 'deepskyblue3') +
    ggrepel::geom_text_repel(data = psub %>% filter(prot %in% mh & median_dsb > 3.5),
                             segment.size = 0.5, segment.color = "grey", size = 5, force = 10, box.padding = 1, nudge_y = 1) +
    ggrepel::geom_text_repel(data = psub %>% filter(prot %in% mh & median_dsb < 3.5),
                             segment.size = 0.5, segment.color = "grey", size = 3, force = 3, color = 'red', nudge_y = -0.4) +
    theme_bw() + 
    ggtitle(cts[i]) + 
    ylim(c(-2,23)) + 
    #geom_segment(data = psub %>% filter(prot %in% mh ), aes(yend=median_dsb, xend=median_neg)) +
    theme(strip.background = element_blank(), axis.title = element_text(size = 0)) + 
    geom_hline(yintercept = 3.5, color = "red") + 
    xlab(" median log10 + 1 protein in empty droplets ") + 
    ylab(" median dsb normalized  protein ") 
  
  ggsave(p, filename = paste0(figpath,cts[i],".pdf"), width = 3, height = 3)
}


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
#   [1] here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1    
# [8] tidyr_1.0.2     tibble_2.1.1    tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             plyr_1.8.4             
# [6] igraph_1.2.4.1          lazyeval_0.2.2          splines_3.5.3           crosstalk_1.0.0         digest_0.6.25          
# [11] foreach_1.4.4           htmltools_0.3.6         lars_1.2                fansi_0.4.0             gdata_2.18.0           
# [16] checkmate_1.9.3         cluster_2.0.7-1         mixtools_1.1.0          ROCR_1.0-7              modelr_0.1.4           
# [21] R.utils_2.8.0           colorspace_1.4-1        ggrepel_0.8.1           rvest_0.3.4             haven_2.1.0            
# [26] xfun_0.7                crayon_1.3.4            jsonlite_1.6            survival_2.43-3         zoo_1.8-6              
# [31] iterators_1.0.10        ape_5.3                 glue_1.3.1              pals_1.5                gtable_0.3.0           
# [36] webshot_0.5.1           kernlab_0.9-27          prabclus_2.3-1          DEoptimR_1.0-8          maps_3.3.0             
# [41] scales_1.0.0            mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1          Rcpp_1.0.1             
# [46] metap_1.1               dtw_1.20-1              viridisLite_0.3.0       xtable_1.8-4            htmlTable_1.13.1       
# [51] reticulate_1.12         foreign_0.8-71          bit_1.1-14              mapproj_1.2.6           proxy_0.4-23           
# [56] mclust_5.4.5            SDMTools_1.1-221.1      Formula_1.2-3           stats4_3.5.3            tsne_0.1-3             
# [61] htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2      fpc_2.2-1              
# [66] ellipsis_0.3.0          acepack_1.4.1           modeltools_0.2-22       ica_1.0-2               pkgconfig_2.0.2        
# [71] R.methodsS3_1.7.1       flexmix_2.3-15          nnet_7.3-12             utf8_1.1.4              labeling_0.3           
# [76] tidyselect_0.2.5        rlang_0.4.5             manipulateWidget_0.10.0 reshape2_1.4.3          later_0.8.0            
# [81] munsell_0.5.0           cellranger_1.1.0        tools_3.5.3             cli_1.1.0               generics_0.0.2         
# [86] broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0            knitr_1.23              bit64_0.9-7            
# [91] fitdistrplus_1.0-14     robustbase_0.93-5       rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1             
# [96] packrat_0.5.0           pbapply_1.4-0           nlme_3.1-137            mime_0.6                R.oo_1.22.0            
# [101] xml2_1.2.0              hdf5r_1.2.0             compiler_3.5.3          rstudioapi_0.10         png_0.1-7              
# [106] lsei_1.2-0              stringi_1.4.3           lattice_0.20-38         ggsci_2.9               vctrs_0.2.4            
# [111] pillar_1.4.1            lifecycle_0.1.0         Rdpack_0.11-0           lmtest_0.9-37           data.table_1.12.2      
# [116] bitops_1.0-6            irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0               
# [121] latticeExtra_0.6-28     promises_1.0.1          KernSmooth_2.23-15      gridExtra_2.3           codetools_0.2-16       
# [126] dichromat_2.0-0         MASS_7.3-51.1           gtools_3.8.1            assertthat_0.2.1        rprojroot_1.3-2        
# [131] withr_2.1.2             diptest_0.75-7          parallel_3.5.3          doSNOW_1.0.16           hms_0.4.2              
# [136] grid_3.5.3              rpart_4.1-13            class_7.3-15            segmented_0.5-4.0       Rtsne_0.15             
# [141] ggpubr_0.2              shiny_1.3.2             lubridate_1.7.4         base64enc_0.1-3  