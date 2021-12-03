suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# set save path 
figpath = here("V2/dsb_process_plots/figures/")
datapath = here("V2/dsb_process_plots/generated_data/")

#######################################
# figure generation 
# read processed data generated above (norm_adt is step 1 only ambient corrected values)
df = readRDS(file = here("V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds"))
norm_adt = readRDS(file = here('V2/dsb_process_plots/generated_data/dsb_norm_adt_mtx.rds'))
iso  = readRDS(file = here("V2/dsb_process_plots/generated_data/isotype_values_dsb.rds"))


# plot noise vector vs lib size
# facet by main lineage (Fig 1)
celltypes = df$celltype_label_1 %>% unique() %>% sort()
tcell = celltypes[c(2,3,4,5,10)]
myeloid = celltypes[c(6,8,9)]
bcell = celltypes[c(1)]
nk = celltypes[c(7)]
plot_sub = df %>%
  mutate(lineage = 
  if_else(celltype_label_1 %in% tcell, "T Cells",
  if_else(celltype_label_1 %in% myeloid, "Myeloid",
  if_else(celltype_label_1 %in% bcell, "B Cells", false = "NK")))) 
plot_sub$lineage = factor(plot_sub$lineage, levels = c("T Cells","Myeloid","B Cells","NK" ))
p = ggplot(plot_sub, aes(x = log10(nUMI_Prot), y = -1*(noise_vector))) + 
  facet_wrap(~lineage, scales = "free", nrow = 1) + 
  theme_bw() + 
  geom_bin2d(bins = 80, show.legend = TRUE) + 
  theme(legend.position = c(0.85, 0.75), legend.key.size = unit(0.3, units = "cm")) + 
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede") + 
  theme(strip.background = element_blank()) + 
  xlab("log10 protein library size") + 
  ylab("dsb Technical Component") + 
  theme(strip.text.x = element_text(size = 12)) + 
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) 
ggsave(p, filename = paste0(figpath, "row_lineage_noise_vector_vs_libsize.pdf"), width = 5.4, height = 3.05)  


# for a supplemental figure (3a) plot the noise vector by each p3 dist cell type 
p = ggplot(df, aes(x = log10(nUMI_Prot), y = -1*(noise_vector))) + 
  facet_wrap(~p3_dist_3, scales = "free", nrow = 2) + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  theme(legend.position = "bottom") +
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede", method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', label.x.npc = 0.1, label.y.npc = 0.9, aes(label = ..r.label..) )+ 
  theme(strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 10)) + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4)) 
ggsave(p, filename = paste0(figpath, "celltype_noise_vector_vs_libsize.pdf"), width = 15.5, height = 5.5)  

# isotype control vs µ1 (Fig 1)
p2 = ggplot(data = df, aes(x = iso_mean, y = mean.1))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.29, units = "cm")) + 
  ggpubr::stat_cor(aes(label = ..r.label..), method = 'pearson') + 
  geom_smooth(color = "black", method = "lm", se = TRUE) + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" mean of isotype controls ") + 
  ylab(" background mean µ1  ")
ggsave(p2, filename = paste0(figpath,"meanisovs_meanbackground.pdf"), width = 4.2, height = 4.3)

# mean 1 mean 2 
# isotype control vs µ1 
p2 = ggplot(data = df, aes(x = mean.1, y = mean.2))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.25, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson") + 
  geom_smooth(color = "black", method = "lm") + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" background mean µ1 ") + 
  ylab(" positive mean µ2  ")
ggsave(p2, filename = paste0(figpath,"mean1_mean2.pdf"), width = 3.3, height = 3.3)

# isotype control vs µ2
p2 = ggplot(data = df, aes(x = iso_mean, y = mean.2))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.25, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson") + 
  geom_smooth(color = "black", method = "lm") + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" mean of isotype controls ") + 
  ylab(" positive mean µ2  ")
ggsave(p2, filename = paste0(figpath,"meanisovs_mean2.pdf"), width = 3.3, height = 3.3)

### distribution of mixture model µ1 and µ2 across cells 
df3 = df %>% gather(protein_population, value, mean.1:mean.2)
p = ggplot(df3, aes(x = value, fill = protein_population)) + 
  geom_histogram(position =  "identity", bins = 100, show.legend = FALSE) + 
  theme_bw() +
  scale_fill_manual(values = c("#3e8ede", "red")) + 
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  xlab(label = "2 component mixture model means") + 
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 15)) +
  ylab(label = "") + 
  xlim(c(-3,15)) 
ggsave(p, filename = paste0(figpath, "FULL_mixture_model_param_dist.pdf"), width = 5, height = 5)

# correlation across latent component.
isotype_value = iso %>% t() %>% as.data.frame()
stopifnot(isTRUE(all.equal(rownames(df), rownames(isotype_value))))
iso_mu = cbind(isotype_value, df$mean.1)
cplot = cor(iso_mu, method = "pearson")
pcor = Hmisc::rcorr(as.matrix(iso_mu), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "Isotype 4", "µ1")

# plot correlation matrix 
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
pdf(file = paste0(figpath,"background_correlation_plot_spearman.pdf"), width = 4,height = 4)
corrplot::corrplot(cplot,
                   method="color", 
                   col=col(200),  
         type="upper",
         addCoef.col = "black", 
         cl.pos = "b", 
         tl.col="black", tl.srt=45, 
         cl.lim = c(0,1),
         diag=FALSE, addgrid.col = "ghostwhite"
)
dev.off()

############### additional correlations discussed in text 
# correlation across latent component.
isotype_value = iso %>% t %>% as.data.frame()
isomu = Matrix::rowMeans(isotype_value)
iso_mu = cbind(isotype_value,
               'isotypemean' = isomu, 
               'mu1' = df$mean.1,
               'plibsize' = df$nUMI_Prot)
cplot = cor(iso_mu, method = "pearson")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "Isotype 4",
                                      "Isotypes Mean", "µ1", "library size")
pdf(file = paste0(figpath,"background_correlation_plot_spearman_extra.pdf"), width = 5,height = 5)
corrplot::corrplot(cplot,
                   method="color", 
                   col=col(200),  
                   type="upper",
                   addCoef.col = "black", 
                   cl.pos = "b", 
                   tl.col="black", tl.srt=45, 
                   cl.lim = c(0,1),
                   diag=F, addgrid.col = "ghostwhite"
)
dev.off()

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
#   [1] Rtsne_0.15          colorspace_1.4-1    ellipsis_0.3.0      class_7.3-15        modeltools_0.2-22  
# [6] ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5        htmlTable_1.13.1    base64enc_0.1-3    
# [11] rstudioapi_0.10     proxy_0.4-23        ggpubr_0.2          npsurv_0.4-0        flexmix_2.3-15     
# [16] bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16   
# [21] splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23         
# [26] Formula_1.2-3       jsonlite_1.6        packrat_0.5.0       broom_0.5.2         ica_1.0-2          
# [31] cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3     
# [36] httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0          
# [41] lars_1.2            acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1     
# [46] gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1         
# [51] cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137       
# [56] iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7           
# [61] rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8     
# [66] MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16      
# [71] parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0       gridExtra_2.3      
# [76] rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       corrplot_0.84      
# [81] foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0      
# [86] SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1     
# [91] bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3    
# [96] bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          R6_2.4.0            generics_0.0.2     
# [101] snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         mgcv_1.8-27         haven_2.1.0        
# [106] pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0     
# [111] survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4       
# [116] hdf5r_1.2.0         KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2  
# [121] metap_1.1           digest_0.6.25       diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3       
# [126] munsell_0.5.0       viridisLite_0.3.0  


