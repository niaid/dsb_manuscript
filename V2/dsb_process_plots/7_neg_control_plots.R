# ambient correction with background (empty) droplets robustness assessment 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here)) 
set.seed(1)

# savepaths 
figpath = here("V2/dsb_process_plots/figures/")
datapath = here("V2/dsb_process_plots/generated_data/")

# define a random subset of proteins across range of background mean values to label 
plabel = c("CD56_PROT", "CD294_PROT", "CD40_PROT", "HLA-DR_PROT", "CD278_PROT", 
           "CD197_PROT", "CD25_PROT","CD62L_PROT", "CD28_PROT", "CD5_PROT", 
           "CD45RA_PROT", "CD8_PROT", "HLA-ABC_PROT", "CD18_PROT", "IgM_PROT", 
           "CD244_PROT", "CD16_PROT", "CD38_PROT", "TCRgd_PROT", "IgD_PROT")

#  load unstained control cell data 
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))

# summarize unstained control data protein stats. 
uadt = 
  un@assay$CITE@raw.data %>% 
  as.matrix() %>% t %>%
  as.data.frame() %>% 
  rownames_to_column("barcode_check") %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_unstained_control = median(count),
            mean_unstained_control = mean(count), 
            var_unstained_control = var(count)
            ) 

# full ADT output (NOT hashing) based subsetting of negative drops. 
adt_neg_full_qc = readRDS(file = here("data/V2_Data/background_data/adt_neg_full_qc.rds"))

#  Empty / background droplets (full) summarize protein stats. 
adt_neg_full_summary =
  adt_neg_full_qc %>% 
  rename(barcode_check = bc) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_negative_droplet = median(count),
            mean_negative_droplet = mean(count),
            var_negative_droplet = var(count)
            ) 

# raw log adt across stained cells summary stats 
raw = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")
cells_logadt = raw@assay$CITE@raw.data %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode_check") %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_cells = median(count),
            mean_cells = mean(count), 
            var_cells = var(count)
  ) 

# calculate each protein background distribution (mean A Supplementary Fig 1)
adtraw = raw@assay$CITE@raw.data
adtlog = log10(adtraw + 1)
library(mclust)
# gaussian mixture across rows (proteins)
prot_model = apply(adtlog, 1, function(x) { mclust::Mclust(x, G=2, warn = TRUE , verbose = TRUE) } )

# extract mean and variance for each protein within each gaussian population 
# mean 1 
p_mean1 = lapply(prot_model, function(x){x$parameters$mean[1]}) %>% unlist()
p_var1 = lapply(prot_model, function(x){x$parameters$variance$sigmasq[1]}) %>% unlist() 

# mean 2 
p_mean2 = lapply(prot_model, function(x){x$parameters$mean[2]}) %>% unlist()
p_var2 = lapply(prot_model, function(x){x$parameters$variance$sigmasq[2]}) %>% unlist() 

# merge and save per protein background 
log_stats = as.data.frame(cbind(p_var1, p_var2, p_mean1, p_mean2)) %>% rownames_to_column("prot")
saveRDS(log_stats,file = paste0(datapath, "logstats.rds"))
# log_stats = readRDS(here("V2/dsb_process_plots/generated_data/logstats.rds"))

# merge data by protein.   
# merge full background with unstained controls 
merged = full_join(adt_neg_full_summary, uadt, by = "prot")
# merge stained cells 
merged = full_join(merged, cells_logadt, by = "prot")
# merge per cell modeled background 
merged = full_join(merged, log_stats, by = "prot")

## plot background mean by unstained and empty 
# plot theme 
corp = list(geom_point(shape = 21, size = 2.8, fill = "red3"), 
  geom_abline(slope = 1, linetype = "dashed"),  
  geom_smooth(method = "lm", color = "black"), 
  ggpubr::stat_cor(method = "pearson"), 
  theme_bw()
)

# mean A vs negative droplet mean 
p1 = ggplot(merged, aes(x = p_mean1, y = mean_negative_droplet)) + 
  corp + 
  xlab("mean of background: stained cells") + 
  ylab(" Empty droplet mean log10 + 1 protein") + 
  ggtitle("negative stained cell mean A vs empty droplets") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
  aes(label = prot), size = 2.5,force = TRUE,fontface = "bold",  segment.size = 0.2, box.padding = 0.5) 

# mean A vs unstained control 
p2 = ggplot(merged, aes(x = p_mean1, y = mean_unstained_control))  +
  corp + 
  xlab("mean of background: stained cells") + 
  ylab("unstained controls mean log10 + 1 protein") + 
  ggtitle("negative stained cell mean A vs unstained control") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
  aes(label = prot), size = 2.5,force = TRUE,fontface = "bold",  segment.size = 0.2, box.padding = 0.5) 
pg = plot_grid(p2,p1, ncol = 1)
ggsave(pg,filename = paste0(figpath, "mean1_cellestimate_vs_neg_unstained.pdf"), height = 8.5, width = 4.5)


# mean in unstained control cells vs mean in empty droplets. 
ndropsfull = ndrops = adt_neg_full_qc$bc %>% unique %>% length()

# mean unstained vs mean empty defined by protein library size distribution for main Fig 1B 
corp2 = list(geom_point(shape = 21, size = 2.8, fill = "red3"), 
             geom_smooth(method = "lm", color = "black"), 
             ggpubr::stat_cor(method = "pearson"), theme_bw()
)
pmain = ggplot(merged, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp2 + 
  ylab("Unstained controls mean log10+1") + 
  xlab("Empty droplets mean log10+1") + 
  theme(axis.title.x = element_text(size = 17)) + 
  theme(axis.title.y = element_text(size = 1)) + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 3,force = TRUE,
                           fontface = "bold", segment.size = 0.2, box.padding = 0.5)  
ggsave(pmain, filename = paste0(figpath, "/mean_stainedvs_unstained.pdf"), height = 4.2, width = 4.2)

# calculate slope 
summary(lm(mean_unstained_control ~ mean_negative_droplet, data = merged))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.02609    0.01152   2.264   0.0261 *  
#   mean_negative_droplet  1.23799    0.01966  62.968   <2e-16 ***


# mean unstained vs mean emty drops for supplemental figure comparing definitions of empty drops 
# comparison of the correlation structure between unstained controls and hashing negatives 
# vs the other definitions of background defined above

# load Hashing Negatives (empty drops defined by hashing)
adt_neg_dmx = readRDS(file = here("data/V2_Data/background_data/adt_neg_dmx.rds")) %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("bc")

#  Hashing Negatives (empty drops defined by hashing) protein stats. 
adt_neg_sub_summary =
  adt_neg_dmx %>% 
  rename(barcode_check = bc) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  mutate(count = log10(count + 1)) %>% 
  group_by(prot) %>% 
  summarize(median_negative_droplet = median(count), 
            mean_negative_droplet = mean(count), 
            var_negative_droplet = var(count)) 

# new merged data with hashing definition of background 
merged2 = full_join(adt_neg_sub_summary, uadt, by = "prot")

# mean in unstained control cells vs mean in empty droplets. 
ndrops = adt_neg_dmx$bc %>% unique %>% length()
p2 = ggplot(merged2, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp + 
  ylab("unstained controls mean log10 + 1 protein") + 
  xlab(" Empty droplet mean log10 + 1 protein") + 
  ggrepel::geom_text_repel(data = merged2 %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 2.5,force = TRUE,fontface = "bold", 
                           segment.size = 0.2, box.padding = 0.5) +  
  ggtitle(paste0(ndrops, " background drops defined by hashing"))


# Compare above with Library size distribution defined background (from adt_neg_full_qc)
# note plotting merged not merged2 here 
p1 = ggplot(merged, aes(x = mean_negative_droplet, y = mean_unstained_control)) + 
  corp + 
  ylab("unstained controls mean log10 + 1 protein") + 
  xlab(" Empty droplet mean log10 + 1 protein") + 
  ggrepel::geom_text_repel(data = merged %>% filter(prot %in% plabel), 
                           aes(label = prot), size = 2.5,force = TRUE,fontface = "bold", 
                           segment.size = 0.2, box.padding = 0.5) +  
  ggtitle(paste0(ndropsfull, " total background drops"))

# combine plot 
p5 = plot_grid(p2, p1, ncol = 1)
ggsave(p5, filename = paste0(figpath, "/mean_stainedvs_unstained_BKGRND_COMPARE.pdf"), height = 8.5, width = 4.5)

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
#   [1] mclust_5.4.5    here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1    
# [9] tidyr_1.0.2     tibble_2.1.1    tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    ellipsis_0.3.0      class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2    
# [8] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        ggpubr_0.2          npsurv_0.4-0        ggrepel_0.8.1      
# [15] flexmix_2.3-15      bit64_0.9-7         fansi_0.4.0         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16   
# [22] splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6       
# [29] packrat_0.5.0       broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0        
# [36] compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0           lars_1.2           
# [43] acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1         
# [50] reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137       
# [57] iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0    
# [64] irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2          
# [71] doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13       
# [78] segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2       
# [85] Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6       
# [92] lattice_0.20-38     ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4         
# [99] R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1       
# [106] foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3         
# [113] modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0         utf8_1.1.4          KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3         
# [120] data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0      
# [127] viridisLite_0.3.0  


