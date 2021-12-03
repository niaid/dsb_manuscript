# Evaluate per cell Gaussian mixture models with k = 1 - 6 components 
set.seed(1)
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(mclust))
'%ni%' = Negate('%in%')

# set save paths 
figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)

# load  batch 1 meta data
b1met = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")@meta.data %>% filter(batch == "1")
# load ambient-corrected (post dsb step I only) data from batch 1
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))

##### mixture model parameter evaluation 
# Fit a gaussian mixture with G =  1 through 6 mixing components to each cell 
# the purpose of this analysis is to evaluate G (i.e. k) that maximizes (using mclust definition of BIC) the BIC
# define cell model data (cmd) list 
cmd = list()
for (i in 1:6) {
  cmd[[i]]  = apply(norm_adt1, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = TRUE)  
    return(g) ; gc()
  })
}

# tidy mixture model parameter estimates for each cell 
# define model result (mr) list
mr = list()
for (i in 1:length(cmd)) {
  # for each k component mixture (i), extract model information including BIC value and
  # format into a list of data frames indexed by k (mixing components)
  print(i)
  mr[[i]]  = lapply(cmd[[i]], broom::glance)
  stopifnot(isTRUE(all.equal(names(mr[[i]]) , colnames(norm_adt1))))
  mr[[i]] = do.call(rbind, mr[[i]])
  # add barcode name 
  mr[[i]]$barcode_check = colnames(norm_adt1)
  # merge this with cell metadata 
  mr[[i]] = full_join(mr[[i]], b1met, by = "barcode_check")
}
# merge into model result (mr) dataframe (df)
# mrdf has 169374 rows corresponding to 6 models fit to 28229 cells
mrdf = do.call(rbind, mr)
saveRDS(mrdf, file = paste0(datapath, "multi_component_model_comparison.rds"))
mrdf = readRDS(file = here('V2/dsb_process_plots/generated_data/multi_component_model_comparison.rds'))

##### plot an arbitrary single fit (the first cell in the matrix) from the 2 component mixture 
# This is inset in Figure 1c
cell = cmd[[2]][[1]]
cell = MclustDR(cell)
pdf(file = paste0(figpath, "examplecell.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()
#### 

# create a reduced result data frame with the mixing component k that 
# optimized (maximized by the convention in mclust) the BIC 
test = mrdf %>% 
  select(BIC, G,barcode_check, p3_dist_3) %>% 
  as.data.frame() %>% 
  arrange(barcode_check) %>% 
  group_by(barcode_check) %>% 
  mutate(bic_ = max(BIC)) %>% 
  filter(bic_ == BIC)
test_best = test %>% select(barcode_check, G_BEST = G)

for (i in 1:6) {
  print(
  (test_best %>% filter(G_BEST == i) %>% nrow) / nrow(test_best)
  )
}
# pct of cells with g = i best fit. 
# [1] 0.0001062737
# [2] 0.8162882
# [3] 0.1724114
# [4] 0.01055652
# [5] 0.000566793
# [6] 7.084913e-05


# extract out all of the fitted means from each Gaussian mixture model k=1 to 6 from each cell
dmean = list()
for (i in 1:length(cmd)) {
  dmean[[i]] = lapply(cmd[[i]], function(x) {x$parameters$mean }) 
  dmean[[i]] = do.call(rbind, dmean[[i]])
  colnames(dmean[[i]]) = paste0("G_",i,"_mu_",colnames(dmean[[i]]))
}
dmean_df = do.call(cbind, dmean) %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode_check")
saveRDS(dmean_df, file = paste0(datapath, "dmean_df.rds"))
dmean_df = readRDS(file = here("V2/dsb_process_plots/generated_data/dmean_df.rds"))

# add BIC values 
BICdf = mrdf %>% 
  select(BIC, G,barcode_check) %>%
  spread(G, BIC) %>% 
  as.data.frame()
colnames(BICdf) = c("barcode_check", "BIC_G1", "BIC_G2", "BIC_G3", "BIC_G4", "BIC_G5", "BIC_G6")
dmean_df = full_join(dmean_df, BICdf, by = 'barcode_check')

# add the optimal BIC 
dmean_test = full_join(dmean_df, test_best, by = "barcode_check")



###############################
# figure generation 

# plot BIC distribution 
g3best = dmean_test %>% filter(G_BEST == 3) %>% 
  mutate(g3_g2_delta = BIC_G3 - BIC_G2) %>%
  mutate(g2_g1_delta = BIC_G2 - BIC_G1)
p1 = ggplot(g3best, aes(x = g3_g2_delta)) + 
  theme_bw() + 
  geom_histogram(bins = 30, fill = "grey", color = "black") + 
  xlab("k=3 - k=2 model") + 
  theme(axis.title.x = element_text(size = 18))
p2 = ggplot(g3best, aes(x = g2_g1_delta )) + 
  theme_bw() + 
  geom_histogram(bins = 30, fill = "grey", color = 'black') + 
  xlab("k=2 - k=1 model") + 
  theme(axis.title.x = element_text(size = 18))
p3 = cowplot::plot_grid(p1,p2, ncol = 2) 
ggsave(p3, filename = paste0(figpath, 'BIC_DELTA_cells_g3best.pdf'), width = 5, height = 2.5)

# models with 3 component gaussian fits:
# plot mean distributions for the 3 components 
dsub = dmean_test %>% filter(G_BEST == 3) %>% 
  select(G_2_mu_1,G_2_mu_2, G_3_mu_1 ,  G_3_mu_2 , G_3_mu_3) %>% 
  gather(parameter, mean, G_2_mu_1:G_3_mu_3) %>% 
  mutate(model = str_sub(parameter, 1,3 ))
p = ggplot(dsub %>% filter(mean < 15), aes(x = mean, fill = parameter)) + 
  theme_bw() + 
  geom_density() +
  theme(strip.background = element_blank()) + 
  facet_wrap(~model) + 
  ggsci::scale_fill_d3(alpha = 0.8)  
ggsave(p, filename = paste0(figpath, "cells_with_G3best_g3vg2.pdf"), width = 4.5, height = 2.5)

# models with 4 component gaussian fits:
# plot mean distributions for the 4 components  
dsub = dmean_test %>% filter(G_BEST == 4) %>% 
  select(G_2_mu_1,G_2_mu_2, G_4_mu_1 ,  G_4_mu_2 , G_4_mu_3, G_4_mu_4) %>% 
  gather(parameter, mean, G_2_mu_1:G_4_mu_4) %>% 
  mutate(model = str_sub(parameter, 1,3 ))
p = ggplot(dsub %>% filter(mean < 15), aes(x = mean, fill = parameter)) + 
  theme_bw() + 
  geom_density() +
  theme(strip.background = element_blank()) + 
  facet_wrap(~model) + 
  ggsci::scale_fill_d3(alpha = 0.8)  
ggsave(p, filename = paste0(figpath, "cells_with_G4_best_g4vg2.pdf"), width = 4.5, height = 2.5)

## best fit distribution per cluster 
# get cell type meta data ctmd for each barcode 
ctmd = b1met %>% select(barcode_check, celltype_label_3, p3_dist_3)
merge = full_join(dmean_test, ctmd, by = "barcode_check")
merge$G_BEST = as.character(merge$G_BEST)

# distribution of optimal BIC values per cluster for supplement 
p = ggplot(data = merge, aes(x = p3_dist_3, fill = G_BEST )) + 
  geom_bar(position = "fill") + ggsci::scale_fill_d3() + 
  theme_minimal()
ggsave(p ,filename = paste0(figpath, "maxbic_celtype.pdf"), width = 5, height = 2.1)

# proportion plot (Global best BIC) for Fig 1 
merge$dv = "dv" # dummy variable for plot 
p = ggplot(merge, aes(x = dv, fill = G_BEST )) + 
  geom_bar(position = "fill", show.legend = FALSE) + 
  ggsci::scale_fill_d3() + 
  ylab("proportion with optimal BIC") + 
  theme_minimal() + 
  theme(axis.text.y =  element_text(size = 6)) + 
  scale_y_continuous(position = "right") + 
  theme(axis.text.x =  element_blank(), axis.title.x = element_blank())
ggsave(p ,filename = paste0(figpath, "maxbic_GLOBAL.pdf"), width = 0.8, height = 3) 

# plot example cell with optimal 3 component BIC 3 mean distributions
cellsub = test_best %>% filter(G_BEST == '3')
cellsub = cellsub$barcode_check[1]
cell = MclustDR(cmd[[3]][[cellsub]])
pdf(file = paste0(figpath,cell, "G3.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()

# plot the 2 component model distributions from the same cell
cell = MclustDR(cmd[[2]][[cellsub]])
pdf(file = paste0(figpath,cell, "G2.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()


## plot global fits showing G=2 is optimal  across cell types (Supplementary Fig 2)
# Mixture model components vs BIC plot 
# saved data 
#mrdf = readRDS("V2/dsb_process_plots/generated_data/multi_component_model_comparison.rds")
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3", "4", "5", "6"))
p = ggplot(mrdf , aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,  show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  facet_wrap(~ p3_dist_3, nrow = 2) + 
  xlab(" number of mixing components in Gaussian mixture model ") + 
  ylab(" BIC ") + 
  theme(axis.title.x =element_text(size = 18)) + 
  theme(axis.title.y =element_text(size = 18)) + 
  theme(strip.text = element_text(size = 14))
ggsave(p, filename = paste0(figpath, "mixture_component_vs_bic.pdf"), width = 15 ,height = 5)  

# same plot, show only global distribution (Fig 1D)
p = ggplot(mrdf , aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,  show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size = 7)) + 
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 0)) 
p
ggsave(p, filename = paste0(figpath, "GLOBALmixture_component_vs_bic.pdf"), width = 2, height = 3.28)  



## mu 1 for G = 2 vs mu1 for G = 3 
# address reviewer comment: test the concordance between the dsb technical component 
# when using a 2-component mixture model vs a 3-component mixture model to define µ1 

# extract the background mean µ1 from the 2 and 3 component fits 
mu1_g2 = lapply(cmd[[2]], function(x){ x$parameters$mean[1] }) %>% unlist()
mu1_g3 = lapply(cmd[[3]], function(x){ x$parameters$mean[1] }) %>% unlist()
d = data.frame( mu1_g2, mu1_g3)
rownames(d) = str_sub(rownames(d), 1,-3)

# regress out the dsb technical component from the ambient corrected data (from dsb step I)
# with 2 versions of the dsb technical component, each defined using k = 2 or 3 component Gaussian mixtures 
# the ambient corrected data is loaded in the session at the top of the script (path below)
# norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))

# define isotype controls 
isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
                             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT" )

# construct noise matrix from both µ1 parameters 
# noise matrix with 2 component model µ1 parameter 
noise_matrix_g2 = t(rbind(norm_adt1[isotype.control.name.vec, ], mu1g2 = d$mu1_g2))
# noise matrix with 3 component model µ1 parameter 
noise_matrix_g3 = t(rbind(norm_adt1[isotype.control.name.vec, ], mu1g3 = d$mu1_g3))

# pc1 regression 
g2pc1 = prcomp(noise_matrix_g2, scale = TRUE)$x[ ,1]
g3pc1 = prcomp(noise_matrix_g3, scale = TRUE)$x[ ,1]
denoised_adt1 = limma::removeBatchEffect(norm_adt, covariates = g2pc1)
denoised_adt2 = limma::removeBatchEffect(norm_adt, covariates = g3pc1)

cormat = cor(t(denoised_adt1), t(denoised_adt2), method = 'pearson')
range(diag(cormat))
# 0.9536114 0.9999648
mean(diag(cormat))
# 0.9918949
median(diag(cormat))
# 0.9949126

d2 = data.frame(g2pc1, g3pc1)

# visualize the distribution and correlaiton of the technical component calculated both ways
p = 
  ggplot(d2, aes(x = g2pc1, y = g3pc1)) + 
  theme_bw() +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(option = "B") + 
  geom_abline(slope = 1, linetype = 'dashed') + 
  xlab('isotype controls combined with \n 2 component mixture µ1') +  
  ylab('isotype controls combined with \n 3 component mixture µ1') + 
  theme(axis.title = element_text(size =15)) + 
  ggpubr::stat_cor(aes(label = ..r.label..)) + 
  theme(legend.position = c(0.85,0.26), legend.key.size = unit(0.5, units = 'cm'))
ggsave(p, filename = paste0(figpath,'tech_component_g2_vs_g3.pdf'), width = 3.7, height = 4.2)

# for future reference, save cmd 
saveRDS(cmd,file = paste0(datapath,'cmd.rds'))


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
# [8] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        ggpubr_0.2          npsurv_0.4-0        flexmix_2.3-15     
# [15] bit64_0.9-7         fansi_0.4.0         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16    splines_3.5.3      
# [22] R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6        packrat_0.5.0      
# [29] broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3     
# [36] httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      limma_3.38.3        cli_1.1.0           lars_1.2           
# [43] acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1         
# [50] reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137       
# [57] iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0    
# [64] irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2          
# [71] doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13       
# [78] segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2       
# [85] Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6       
# [92] lattice_0.20-38     ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    ggsci_2.9          
# [99] plyr_1.8.4          R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0        
# [106] pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12        
# [113] tsne_0.1-3          modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0         utf8_1.1.4          KernSmooth_2.23-15  readxl_1.3.1       
# [120] grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3       
# [127] munsell_0.5.0       viridisLite_0.3.0  