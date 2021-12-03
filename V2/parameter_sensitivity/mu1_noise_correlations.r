suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here)) 
suppressMessages(library(Seurat)) 
suppressMessages(library(mclust)) 
set.seed(1)
####### 

figpath = here("V2/parameter_sensitivity/figures_mu1noise/"); dir.create(figpath)
datapath = here("V2/parameter_sensitivity/generateddata_mu1noise/"); dir.create(datapath)

# load dsb step 1 ambient corrected data from b1 
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))
h = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")
isotypes = c("MouseIgG1kappaisotype_PROT", "MouseIgG2akappaisotype_PROT",
             "Mouse IgG2bkIsotype_PROT",  "RatIgG2bkIsotype_PROT")

##### apply per cell mixture model  
cmd = apply(norm_adt1, 2, function(x){ g = Mclust(x, G=2, warn = TRUE , verbose = TRUE) }) 

# extract vector of all mu1 and mu2 proteins for each cell
m2p = lapply(cmd, function(x) x$classification[x$classification == 2] %>% names())
m1p = lapply(cmd, function(x) x$classification[x$classification == 1] %>% names())
names(m1p) = names(cmd)

# get 100 samples of k=4 proteins for each cell i 
psamples = list()
for (i in 1:length(m1p)) {
  
  # define vector of µ1 proteins from cell i, not including isotypes 
  protvec = m1p[[i]]
  protvec = setdiff(protvec, isotypes)
  
  # 100 samples of k=4 proteins from the population of proteins comprising µ1 
  psamples[[i]] = replicate(100, expr = sample(x = protvec,size = 4, replace = FALSE), simplify = FALSE)
}

# get the mean of the k protein samples from mu 1 for each i for each n  
mx = list()
for (i in 1:length(psamples)) {
  pv = psamples[[i]]
  cell_n = names(cmd)[i]
  mx[[i]] = lapply(pv, function(x){mean(norm_adt1[x,cell_n])})
}

# vector of the resampled mu1 from the k=4 subsamples
mu1_ls = lapply(mx, unlist)
mu1_df = do.call(cbind, mu1_ls) %>% as.matrix() %>% as.data.frame()
colnames(mu1_df) = names(cmd)
saveRDS(mu1_df, file = paste0(datapath, "mu1_df.rds"))

# load cell metadata containing true values mu 2 
md = readRDS(file = "V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds")
mdsub = md %>% filter(barcode_check %in% names(cmd)) 


# correlate mu1 protein subsample means (randomly sampled 0100 times) 
# with µ2 for each of the 100 sub samples across single cells 
rho2 = list()
yy2 = mdsub$mean.2
for (i in 1:nrow(mu1_df)) {
  xx = mu1_df[i,] %>% as.matrix() %>% as.vector
  rho2[[i]] = cor.test(x = xx,  y = yy2, method = "pearson")$estimate
}
rv2 = unlist(rho2, use.names = FALSE) %>% as.vector()
saveRDS(rv2, file = paste0(datapath, "rv2.rds"))

# correlate isotype control mean with full mu 1
yy = Matrix::colMeans(norm_adt1[isotypes, ]) %>% as.matrix() %>% as.vector()
m2 = mdsub$mean.2
stopifnot(isTRUE(all.equal(colnames(norm_adt1), mdsub$barcode_check)))
ct2 =  cor.test(m2, yy, method = "pearson")$estimate
ct2
# cor 
# 0.2576368 


# visualization and comparison with mu2 and isotype controls
pdf(file = paste0(figpath, "pearson_mu1_cor_m2.pdf"),height = 4, width = 4.5)
rethinking::dens(rv2, 
                 col = '#3e8ede', 
                 adj = 0.5,
                 lwd = 3, 
                 show.HPDI = 0.50, 
                 font.main= 1, cex.main = 0.8,
                 xlim = c(0.24, 0.48), xlab ='Pearson Correlation',
                 main = '100 random samples of k = 4 µ1 proteins from n = 28229 cells: \n Distribution of µ1 k-sample means Pearson correlation w/ µ2 (blue)', 
)
grid(lty = 1, lwd = 1, col = "grey")
abline(v = ct2, col="red", lwd=3, lty=2)
legend('bottom', col = 'red',lty = 2, lwd = 2, cex = 0.7,
       legend = 'Pearson correlation \n isotype controls and µ2')
dev.off()


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
#   [1] mclust_5.4.5    Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5    
# [10] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1         snow_0.4-3           backports_1.1.4      Hmisc_4.2-0          plyr_1.8.4           igraph_1.2.4.1       lazyeval_0.2.2      
# [8] splines_3.5.3        inline_0.3.15        digest_0.6.25        foreach_1.4.4        htmltools_0.3.6      lars_1.2             rethinking_2.12     
# [15] gdata_2.18.0         checkmate_1.9.3      cluster_2.0.7-1      mixtools_1.1.0       ROCR_1.0-7           modelr_0.1.4         matrixStats_0.54.0  
# [22] R.utils_2.8.0        prettyunits_1.0.2    colorspace_1.4-1     rvest_0.3.4          haven_2.1.0          xfun_0.7             callr_3.2.0         
# [29] crayon_1.3.4         jsonlite_1.6         survival_2.43-3      zoo_1.8-6            iterators_1.0.10     ape_5.3              glue_1.3.1          
# [36] gtable_0.3.0         pkgbuild_1.0.3       kernlab_0.9-27       rstan_2.19.3         shape_1.4.4          prabclus_2.3-1       DEoptimR_1.0-8      
# [43] scales_1.0.0         mvtnorm_1.0-10       bibtex_0.4.2         Rcpp_1.0.1           metap_1.1            dtw_1.20-1           htmlTable_1.13.1    
# [50] reticulate_1.12      foreign_0.8-71       bit_1.1-14           proxy_0.4-23         SDMTools_1.1-221.1   Formula_1.2-3        StanHeaders_2.21.0-1
# [57] stats4_3.5.3         tsne_0.1-3           htmlwidgets_1.3      httr_1.4.0           gplots_3.0.1.1       RColorBrewer_1.1-2   fpc_2.2-1           
# [64] acepack_1.4.1        modeltools_0.2-22    ica_1.0-2            loo_2.3.1            pkgconfig_2.0.2      R.methodsS3_1.7.1    flexmix_2.3-15      
# [71] nnet_7.3-12          tidyselect_0.2.5     rlang_0.4.5          reshape2_1.4.3       munsell_0.5.0        cellranger_1.1.0     tools_3.5.3         
# [78] cli_1.1.0            generics_0.0.2       broom_0.5.2          ggridges_0.5.1       npsurv_0.4-0         processx_3.3.1       knitr_1.23          
# [85] bit64_0.9-7          fitdistrplus_1.0-14  robustbase_0.93-5    caTools_1.17.1.2     RANN_2.6.1           packrat_0.5.0        pbapply_1.4-0       
# [92] nlme_3.1-137         R.oo_1.22.0          xml2_1.2.0           hdf5r_1.2.0          compiler_3.5.3       rstudioapi_0.10      png_0.1-7           
# [99] lsei_1.2-0           stringi_1.4.3        ps_1.3.0             lattice_0.20-38      vctrs_0.2.4          pillar_1.4.1         lifecycle_0.1.0     
# [106] Rdpack_0.11-0        lmtest_0.9-37        data.table_1.12.2    bitops_1.0-6         irlba_2.3.3          gbRd_0.4-11          R6_2.4.0            
# [113] latticeExtra_0.6-28  KernSmooth_2.23-15   gridExtra_2.3        codetools_0.2-16     MASS_7.3-51.1        gtools_3.8.1         assertthat_0.2.1    
# [120] rprojroot_1.3-2      withr_2.1.2          diptest_0.75-7       parallel_3.5.3       doSNOW_1.0.16        hms_0.4.2            grid_3.5.3          
# [127] rpart_4.1-13         coda_0.19-2          class_7.3-15         segmented_0.5-4.0    Rtsne_0.15           lubridate_1.7.4      base64enc_0.1-3  


