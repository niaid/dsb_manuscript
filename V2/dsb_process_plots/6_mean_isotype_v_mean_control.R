suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# apply dsb in steps and save intermediate model fits and other results for correlation among latent variables. 
# set save path 
figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)
project_title = "Healthy Donor prevaccination"

# load raw data 
h1 = readRDS(here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds"))
neg = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# separate by batch 
negadt1 = neg[[1]]; negadt2 = neg[[2]]
b1cells = h1@meta.data %>% filter(batch == "1") %$% barcode_check
b2cells = h1@meta.data %>% filter(batch == "2") %$% barcode_check
adt = h1@assay$CITE@raw.data %>% as.matrix()
adt1 = adt[ ,b1cells]
adt2 = adt[ ,b2cells]

# DSB normalize in steps over the 2 batches 
# Step I ambient correction 
pseudocount.use = 10
adtu_log1 = log(negadt1 + pseudocount.use) 
adt_log1 = log(adt1 + pseudocount.use)
adtu_log2 = log(negadt2 + pseudocount.use) 
adt_log2 = log(adt2 + pseudocount.use) 

# batch 1 rescaling 
mu_u1 = apply(adtu_log1, 1 , mean)
sd_u1 = apply(adtu_log1, 1 , sd)
norm_adt1 = apply(adt_log1, 2, function(x) (x  - mu_u1) / sd_u1) 

# batch 2 rescaling 
mu_u2 = apply(adtu_log2, 1 , mean)
sd_u2 = apply(adtu_log2, 1 , sd)
norm_adt2 = apply(adt_log2, 2, function(x) (x  - mu_u2) / sd_u2) 

# merge adt normalized by batch 
norm_adt = cbind(norm_adt1, norm_adt2)

# save ambient corrected (but not yet step II denoised) data 
saveRDS(norm_adt, file = paste0(datapath, "dsb_norm_adt_mtx.rds"))
saveRDS(norm_adt1, file = paste0(datapath, "b1_dsb_norm_adt_mtx.rds"))
saveRDS(norm_adt2, file = paste0(datapath, "b2_dsb_norm_adt_mtx.rds"))

# visualize the relationship between the per cell technical component and protein library zise 
# run per-cell gaussian mixture model (on each batch)
library(mclust)
cellwise_model1 = apply(norm_adt1, 2, function(x) {
			g = Mclust(x, G=2, warn = TRUE , verbose = TRUE)  
			return(g) 
		})
cellwise_model2 = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = TRUE)  
  return(g) 
})

# create tm tidy model fit data (extract out the mean of each population)
cm = c(cellwise_model1, cellwise_model2)
tm  = lapply(cm, function(x){broom::tidy(x)[ 2, 5:6]}) %>% bind_rows()
# create model result mr data frame with BIC for each cell 
mr  = lapply(cm, broom::glance) %>% bind_rows()
tm$barcode_check = names(cm); mr$barcode_check = names(cm)
md1 = full_join(tm, mr, by = "barcode_check")

# merge model results () with metadata 
md = full_join(md1,h1@meta.data,  by = "barcode_check")

# calculate latent component noise vector 
cellwise_background_mean = lapply(cm, function(x) {x$parameters$mean[1] })
cellwise_background_mean = unlist(cellwise_background_mean, use.names = FALSE)
cellwise_positive_mean = lapply(cm, function(x) {x$parameters$mean[2] })
cellwise_positive_mean = unlist(cellwise_positive_mean, use.names = FALSE)

# background mean is in tidy model data above tm 
all.equal(cellwise_background_mean, tm$mean.1)

# dsb step II 
# define pc1 through isotypes and background protein as a latent variable 
isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
                             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT" )
noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
get_noise_vector = function(noise_matrix) { 
  g = prcomp(t(noise_matrix), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
df = full_join(md, PC, by = "barcode_check")

# add library size 
nUMI_Prot = h1@assay$CITE@raw.data %>% colSums()
df$nUMI_Prot = nUMI_Prot

# add isotype control means 
isotypes = c("MouseIgG1kappaisotype_PROT","MouseIgG2akappaisotype_PROT", 
             "Mouse IgG2bkIsotype_PROT", "RatIgG2bkIsotype_PROT")
iso = norm_adt[ isotypes,  ]
iso_mean = apply(iso, 2, mean)

# mean comparison 
df = cbind(iso_mean, df)
saveRDS(df ,file = paste0(datapath, "mixture_model_metadata_merged.rds"))
saveRDS(iso, file = paste0(datapath, "isotype_values_dsb.rds"))

# these are the values in mixture_model_metadata_merged.rds
# [1] "iso_mean"                       "mean.1"                         "mean.2"                         "barcode_check"                 
# [5] "model"                          "n"                              "G"                              "BIC"                           
# [9] "logLik"                         "df"                             "hypvol"                         "nGene"                         
# [13] "nUMI"                           "orig.ident"                     "pctMT"                          "tenx_lane"                     
# [17] "batch"                          "hto_classification"             "hto_classification_global"      "DEMUXLET.RD.PASS"              
# [21] "DEMUXLET.N.SNP"                 "DMX_GLOBAL_BEST"                "DEMUXLET.BARCODE"               "sampleid"                      
# [25] "joint_classification_global"    "dmx_hto_match"                  "total_features_by_counts"       "log10_total_features_by_counts"
# [29] "total_counts"                   "log10_total_counts"             "pct_counts_in_top_50_features"  "p3_dist_1"                     
# [33] "p3_dist_2"                      "p3_dist_3"                      "p3_dist_4"                      "K0"                            
# [37] "K1"                             "K2"                             "K3"                             "celltype_label_3"              
# [41] "celltype_label_1"               "noise_vector"                   "nUMI_Prot"  

# save another version with a better name for the variable `noise_vector` which is the 
# the dsb technical component as calculated by dsb 
df2 = df
df2 = df2 %>% rename(dsb_technical_component = noise_vector) %>% 
  select(dsb_technical_component, nUMI_Prot, everything()) 
saveRDS(df2,file = paste0(datapath, "pbmc_metadata_modelfits_2component.rds"))

# write stats table for table (Fig 2I)  
df = readRDS(file = here('V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds'))
iso = readRDS(file = here('V2/dsb_process_plots/generated_data/isotype_values_dsb.rds'))
d = data.frame(
  dataset = project_title,
  nprot = nrow(norm_adt), 
  ncells = ncol(norm_adt), 
  n_background = ncol(cbind(neg[[1]], neg[[2]])), 
  cor_tech_size = cor(log10(df$nUMI_Prot), -1*(df$noise_vector)), 
  cor_isotype_background = cor(df$iso_mean, df$mean.1)
)
data.table::fwrite(d, file = paste0(datapath, project_title, "TABLE_STATS.txt"), sep = "\t")

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
#   [1] here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    tidyverse_1.2.1
# [11] Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5        htmlTable_1.13.1   
# [9] base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7         lubridate_1.7.4     mvtnorm_1.0-10     
# [17] xml2_1.2.0          codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3      
# [25] jsonlite_1.6        packrat_0.5.0       broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0        
# [33] compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0           lars_1.2            acepack_1.4.1      
# [41] htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1         
# [49] cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11        
# [57] lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1      
# [65] zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0      
# [73] gridExtra_2.3       rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2   
# [81] bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6       
# [89] lattice_0.20-38     ROCR_1.0-7          htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    plyr_1.8.4          R6_2.4.0            generics_0.0.2     
# [97] snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14
# [105] mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0         KernSmooth_2.23-15 
# [113] readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25       diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3       
# [121] munsell_0.5.0      