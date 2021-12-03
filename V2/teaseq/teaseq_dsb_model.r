# R4 Seurat 4
'%ni%' = Negate('%in%')
suppressMessages(library(dsb))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# test applicability of the ambient component and technical component with multiple model fitting comparisons
# project info 
project_title = "TEA-seq"
figpath = paste0(here("V2/teaseq/figures_v2/"), project_title, "/")
datapath = paste0(here("V2/teaseq/generated_data_v2/"), project_title, "/")

# load object 
s = readRDS(file = here("V2/teaseq/generated_data_v2/TEA-seq/full_teaseq_r4s4_object_processed.rds"))

# load raw data from cells and backgrond drops 
cells_mtx_rawprot = as.matrix(s@assays$CITE@counts)
negative_mtx_rawprot = readRDS(file = here('V2/teaseq/generated_data_v2/TEA-seq/adt_background.rds'))

########################
##DSB process plots 
pseudocount.use = 10
adtu_log = log(negative_mtx_rawprot + pseudocount.use) 
adt_log = log(cells_mtx_rawprot + pseudocount.use)

# apply ambient correction
mu_u = apply(adtu_log, 1, mean)
sd_u = apply(adtu_log, 1, sd)
norm_adt2 = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u) 

# gaussian mixture on ambient corrected values 
library(mclust)
cellwise_model = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = FALSE)  
  return(g) 
})

# tidy Gaussian mixture model data 
tm = lapply(cellwise_model, function(x){broom::tidy(x)})
stopifnot(isTRUE(all.equal(colnames(adt_log), names(tm))))
tm = do.call(rbind, tm)
tm = tm %>% filter(component == 1)
tm$barcode_check = colnames(adt_log)
# get model result (mr); with BIC values 
mr = lapply(cellwise_model, broom::glance)
mr = do.call(rbind, mr)
mr$barcode_check = colnames(adt_log)

# merge model results with metadata 
d = s@meta.data
d$barcode_check = rownames(d)
#df_dsb$barcode_check = rownames(d)
md = full_join(d, mr, by = "barcode_check")
fullmd = full_join(md, tm, by = "barcode_check")

# calculate latent component noise vector 
cellwise_background_mean = lapply(cellwise_model, function(x) {x$parameters$mean[1] }) %>% unlist(use.names = FALSE) 
stopifnot(isTRUE(all.equal(cellwise_background_mean, tm$mean)))

# define latent noise component 
noise_matrix = rbind(norm_adt2['IgG1-K-Isotype-Control', ], cellwise_background_mean)
get_noise_vector = function(x) { 
  g = prcomp(t(x), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
fullmd = full_join(fullmd, PC, by = "barcode_check")

# noise vector vs protein library size clusters 
p = ggplot(fullmd, aes(x = log10(nCount_CITE), y = -1*(noise_vector))) + 
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', aes(label = ..r.label..)) + 
  geom_bin2d(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(strip.background = element_blank()) + 
  geom_smooth(color = "#3e8ede", method = 'lm') + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) + 
  facet_wrap(~dsb_pseudo_knn_res.2, scales = 'free')

ggsave(p, filename = paste0(figpath, project_title, "ALTBINnoise_vector_vs_libsize.pdf") , width = 6, height = 6)

# isotype control vs mixture mean 1 
iso = norm_adt2['IgG1-K-Isotype-Control',  ]

# add isotype mean to metadata.  
# here this is not as applicable a comparison because there is not 
# shared variation acoss the isotypes there is only isotype control 
# The mean just represents the single isotype correlation with µ1
fullmd = cbind(fullmd, iso)
fullmd = fullmd %>% rename(mean.1 = mean)

# Save cell metadata created 
data.table::fwrite(fullmd, file = paste0(datapath, project_title, '_cellmetadata.txt'),sep = "\t")

 #####################
# mixture model parameter evaluation 
cmdl = list()
for (i in 1:3) {
  cmdl[[i]] = apply(norm_adt2, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = FALSE)  
    return(g) ; gc()
  })
}
stopifnot(all.equal(names(cmdl[[1]]), colnames(adt_log)))
mr2 = list()
for (i in 1:length(cmdl)) {
  mr2[[i]] = lapply(cmdl[[i]], broom::glance)
  mr2[[i]] = do.call(rbind, mr2[[i]])
  mr2[[i]]$barcode_check = colnames(adt_log)
}
mrdf = do.call(rbind, mr2)

# add cluster metadata from clustering above 
cmd = d %>% select(barcode_check, dsb_pseudo_knn_res.2)
mrdf = full_join(mrdf, cmd, by = "barcode_check")

# save mixture fit dataframe 
write_delim(mrdf, file  = paste0(datapath, project_title, "_mixture_model_comparison.txt"),delim = "\t")

# calculate percent of cells with G = 2 best fit 
test = mrdf %>% 
  select(BIC, G,barcode_check) %>% 
  as.data.frame() %>% 
  arrange(barcode_check) %>% 
  group_by(barcode_check) %>% 
  mutate(bic_ = max(BIC)) %>% 
  filter(bic_ == BIC)
test_best = test %>% select(barcode_check, G_BEST=  G)
table(test_best$G_BEST)
# 1     2     3 
# 81 22226  6654 


# plot g vs BIC 
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3"))

###############################
# write stats table for s table.  
d = data.frame(
  dataset = project_title,
  nprot = nrow(cells_mtx_rawprot), 
  ncells = ncol(cells_mtx_rawprot), 
  n_background = ncol(negative_mtx_rawprot), 
  cor_tech_size = cor(log10(fullmd$nCount_CITE), -1*(fullmd$noise_vector)), 
  cor_isotype_background = cor(fullmd$iso, fullmd$mean.1)
)
data.table::fwrite(d, file = paste0(datapath, project_title, "TABLE_STATS.txt"), sep = "\t")

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
#   [1] mclust_5.4.7       here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
# [9] tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       magrittr_2.0.1     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2        splines_4.0.5         listenv_0.8.0        
# [8] scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1     fansi_0.4.2           tensor_1.5            cluster_2.1.2         ROCR_1.0-11          
# [15] openxlsx_4.2.3        limma_3.46.0          globals_0.14.0        modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0 colorspace_2.0-0     
# [22] rvest_0.3.6           ggrepel_0.9.1         haven_2.3.1           crayon_1.4.1          jsonlite_1.7.2        spatstat.data_2.1-0   survival_3.2-10      
# [29] zoo_1.8-8             glue_1.4.2            polyclip_1.10-0       gtable_0.3.0          leiden_0.3.7          car_3.0-10            future.apply_1.7.0   
# [36] abind_1.4-5           scales_1.1.1          DBI_1.1.1             rstatix_0.7.0         miniUI_0.1.1.1        Rcpp_1.0.6            viridisLite_0.3.0    
# [43] xtable_1.8-4          reticulate_1.18       spatstat.core_2.0-0   foreign_0.8-81        htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2   
# [50] ellipsis_0.3.1        ica_1.0-2             pkgconfig_2.0.3       farver_2.0.3          uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10        
# [57] utf8_1.1.4            tidyselect_1.1.0      labeling_0.4.2        rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0        
# [64] cellranger_1.1.0      tools_4.0.5           cli_2.5.0             generics_0.1.0        broom_0.7.5           ggridges_0.5.3        fastmap_1.1.0        
# [71] goftest_1.2-2         fs_1.5.0              fitdistrplus_1.1-3    zip_2.1.1             RANN_2.6.1            pbapply_1.4-3         future_1.21.0        
# [78] nlme_3.1-152          mime_0.10             xml2_1.3.2            compiler_4.0.5        rstudioapi_0.13       plotly_4.9.3          curl_4.3             
# [85] png_0.1-7             ggsignif_0.6.0        spatstat.utils_2.1-0  reprex_1.0.0          stringi_1.5.3         lattice_0.20-41       Matrix_1.3-2         
# [92] ggsci_2.9             vctrs_0.3.6           pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18     
# [99] data.table_1.14.0     cowplot_1.1.1         irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       R6_2.5.0              promises_1.2.0.1     
# [106] KernSmooth_2.23-18    gridExtra_2.3         rio_0.5.16            parallelly_1.23.0     codetools_0.2-18      MASS_7.3-53.1         assertthat_0.2.1     
# [113] rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           parallel_4.0.5        hms_1.0.0             grid_4.0.5           
# [120] rpart_4.1-15          carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0          shiny_1.6.0           lubridate_1.7.9.2   
