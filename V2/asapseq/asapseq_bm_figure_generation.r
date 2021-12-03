suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(here))
suppressMessages(library(dsb))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

# asap seq data 
figpath = here('V2/asapseq/figures/'); dir.create(figpath)
datapath = here('V2/asapseq/generated_data/'); dir.create(datapath)

# project 
project_title = "ASAPseq"

# read seurat object  
s = readRDS(here("V2/asapseq/generated_data/s_asapseqbm_processed_protein_seurat4.rds"))
cells_mtx_rawprot = readRDS(file = here('V2/asapseq/generated_data/cells_mtx_rawprot.rds'))
negative_mtx_rawprot = readRDS(file = here('V2/asapseq/generated_data/negative_mtx_rawprot.rds'))

# labeled umap 
p = AugmentPlot(DimPlot(s, reduction = 'umapdsb',
                        repel = TRUE, label.size = 10.5, pt.size = 0.2,
                        cols = BuenColors::jdb_palette(name = 'lawhoops'), 
                        group.by = 'CITE_snn_res.1.5', label = TRUE)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank()) +
  NoLegend() + 
  xlab('dsb UMAP 1') + 
  ylab('dsb UMAP 2') + 
  ggtitle('ASAP-seq Bone Marrow')
p
ggsave(p, filename = paste0(figpath, project_title, "clusters.pdf"), width = 3.3, height = 3.3)

# MG plots to show zero centering  
d = cbind(s@meta.data, as.data.frame(t(s@assays$CITE@data)),  s@reductions$umapdsb@cell.embeddings)
d2 = cbind(s@meta.data, as.data.frame(t(s@assays$CLR@data)),  s@reductions$umapdsb@cell.embeddings)
p1 = ggplot(d, aes( x =  `CD3-2` ,y = `CD4-1` )) + 
  theme_bw() +
  geom_bin2d(bins = 300, show.legend = F) +
  geom_vline(xintercept = 0) + geom_vline(xintercept = 3.5, color = "red")  + 
  geom_hline(yintercept =  0) + geom_hline(yintercept = 3.5, color = "red")  + 
  scale_fill_viridis_c(option = "B") + 
  ggtitle('dsb normalized')
ggsave(p1, filename = paste0(figpath, project_title, "CD4_CD3.pdf"), width = 3, height = 3.2)
p2 = ggplot(d2, aes( x =  `CD3-2` ,y = `CD4-1` )) + 
  theme_bw() + 
  geom_bin2d(bins = 300, show.legend = F) +
  geom_vline(xintercept = 0)  +
  geom_hline(yintercept =  0) +
  scale_fill_viridis_c(option = "B") + 
  ggtitle("CLR normalized")
ggsave(p2, filename = paste0(figpath, project_title, "clr_CD4_CD3.pdf"), width = 3, height = 3.2)

########################
##DSB process plots 
pseudocount.use = 10
adtu_log = log(negative_mtx_rawprot + pseudocount.use) 
adt_log = log(cells_mtx_rawprot + pseudocount.use)

# apply ambient correction 
mu_u = apply(adtu_log, 1, mean)
sd_u = apply(adtu_log, 1, sd)
norm_adt2 = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u) 

# gaussian mixture on background ambient corrected matrix 
library(mclust)
cellwise_model = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = FALSE)  
  return(g) 
})

# tidy Gaussian mixture model data 
# slight difference in this code for mclust version 
# in R4.0.5 vs 3.5.1 as done in the dsb_process scripts 
# because of how broom formats the mclust model data in R4.0.5
tm  = lapply(cellwise_model, function(x){broom::tidy(x)})
tm = do.call(rbind, tm)
tm = tm %>% filter(component == 1)
tm$barcode_check = colnames(adt_log)
# format model result (mr) data with BIC values 
mr  = lapply(cellwise_model, broom::glance)
mr = do.call(rbind, mr)
mr$barcode_check = colnames(adt_log)

# merge model results with metadata and normalized protein values in `d`
# repeat for background mean in tm 
d$barcode_check = rownames(d)
md = full_join(d, mr, by = "barcode_check")
fullmd = full_join(md, tm, by = "barcode_check")

# calculate latent component noise vector 
cellwise_background_mean = lapply(cellwise_model, function(x) {x$parameters$mean[1] }) %>% unlist(use.names = FALSE) 
stopifnot(isTRUE(all.equal(cellwise_background_mean, tm$mean)))

# define latent noise component 
isotypes = rownames(norm_adt2)[grepl(x = rownames(norm_adt2),pattern =  "sotype")]
noise_matrix = rbind(norm_adt2[isotypes, ], cellwise_background_mean)
get_noise_vector = function(x) { 
  g = prcomp(t(x), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
fullmd = full_join(fullmd, PC, by = "barcode_check")

# noise vector vs protein library size clusters 
p = ggplot(fullmd, aes(x = log10_prot_size, y = noise_vector)) + 
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', aes(label = ..r.label..)) + 
  geom_bin2d(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(strip.background = element_blank()) + 
  theme(strip.text = element_text(size=16)) + 
  geom_smooth(color = "#3e8ede", method = 'lm') + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) + 
  facet_wrap(~seurat_clusters, scales = 'free', nrow = 2)
ggsave(p, filename = paste0(figpath, project_title, "ALTBINnoise_vector_vs_libsize.pdf") , width = 12, height = 5)

# isotype control vs mixture mean 1 
iso = norm_adt2[isotypes,  ]
iso_mean = data.frame(iso_mean = apply(iso, 2, mean)) %>% rownames_to_column('barcode_check')

# add isotype mean to metadata.  
fullmd = full_join(fullmd, iso_mean, by = 'barcode_check')
fullmd = fullmd %>% rename(mean.1 = mean)

# plot 
p = ggplot(data = fullmd, 
           aes(x = iso_mean, y = mean.1)) + 
  geom_bin2d(bins = 180, show.legend = FALSE) + 
  ggpubr::stat_cor(method = "pearson") + 
  scale_fill_viridis_c(option = "B") + 
  theme_bw() + 
  xlab("isotype means") + 
  ylab(" µ1 background mean ") + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title,"meanisovs_meanbackground.pdf"), width = 3, height = 3)


# correlation across latent component. 
isotype_value = iso %>% t() %>% as.data.frame() %>% rownames_to_column("bc")
mean1df  = fullmd %>% select(mean.1, bc = barcode_check)
dcor = full_join(isotype_value, mean1df, by = "bc") %>% column_to_rownames("bc")
cplot = cor(dcor, method = "spearman")
pcor = Hmisc::rcorr(as.matrix(dcor), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "Isotype 4", "µ1")
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
pdf(file = paste0(figpath,project_title,"background_correlation_plot_spearman.pdf"), width = 4,height = 4)
corrplot::corrplot(cplot,method="color", col=col(200),  
                   cl.lim = c(0,1),
                   type="upper",
                   order = "original",
                   addCoef.col = "black", 
                   cl.pos = "b", 
                   tl.col="black", tl.srt=45, 
                   diag=FALSE, addgrid.col = "ghostwhite"
)
dev.off()

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
# format new model results 
mr2 = list()
for (i in 1:length(cmdl)) {
  mr2[[i]] = lapply(cmdl[[i]], broom::glance)
  mr2[[i]] = do.call(rbind, mr2[[i]])
  mr2[[i]]$barcode_check = colnames(adt_log)
}
mrdf = do.call(rbind, mr2)
# add cluster metadata cmd from clustering above 
cmd = d %>% select(barcode_check, seurat_clusters)
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
test_best$G_BEST %>% table 
9386 / (41 +1500+ 9386)
#0.8589732 prop with k3 optimal fit in theis dataset with many proetins 
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3"))


#  test correlation with the background mean of the 3 component model 
mu1_g3 = lapply(cmdl[[3]], function(x){ x$parameters$mean[1] }) %>% unlist()
names(mu1_g3) = str_replace(string = names(mu1_g3), pattern = ".1", replacement = "")
noise_matrix2 = rbind(norm_adt2[isotypes, ], mu1_g3)
get_noise_vector = function(x) { 
  g = prcomp(t(x), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector2 = get_noise_vector(noise_matrix2)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector2) %>% rownames_to_column("barcode_check")

# noise vector vs protein library size clusters 
prot_lib_size = colSums(cells_mtx_rawprot) %>% as.data.frame() 
colnames(prot_lib_size) = "nUMI_Prot"
prot_lib_size = prot_lib_size %>% rownames_to_column("barcode_check")
md = full_join(PC, prot_lib_size, by = "barcode_check")
cor_tech_size2 = cor(log10(md$nUMI_Prot),  md$noise_vector)
cor_tech_size2
# 0.5073412

# above shows the correlation with the library size is lower with 
# the 3 component gaussian mixture derived µ1 parameter for the dsb technical component
# despite the more optimal fit with the g3 model suggesting the 2 component fit is capturing 
# the latent component across cells when combined with the shared variaiton of istype controls 


###############################
# write stats table for s table.  
d = data.frame(
  dataset = project_title,
  nprot = nrow(cells_mtx_rawprot), 
  ncells = ncol(cells_mtx_rawprot), 
  n_background = ncol(negative_mtx_rawprot), 
  cor_tech_size = cor(fullmd$log10_prot_size, fullmd$noise_vector), 
  cor_isotype_background = cor(fullmd$iso_mean, fullmd$mean.1)
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
#   [1] mclust_5.4.7       magrittr_2.0.1     forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0       
# [8] tidyr_1.1.2        tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    dsb_0.1.0          here_1.0.1         Matrix_1.3-2      
# [15] data.table_1.14.0  SeuratObject_4.0.0 Seurat_4.0.1      
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       Hmisc_4.5-0           corrplot_0.84         plyr_1.8.6            igraph_1.2.6         
# [7] lazyeval_0.2.2        splines_4.0.5         listenv_0.8.0         scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1    
# [13] viridis_0.5.1         checkmate_2.0.0       tensor_1.5            cluster_2.1.2         ROCR_1.0-11           openxlsx_4.2.3       
# [19] limma_3.46.0          globals_0.14.0        modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0 jpeg_0.1-8.1         
# [25] colorspace_2.0-0      rvest_0.3.6           ggrepel_0.9.1         xfun_0.21             haven_2.3.1           crayon_1.4.1         
# [31] jsonlite_1.7.2        spatstat.data_2.1-0   survival_3.2-10       zoo_1.8-8             glue_1.4.2            polyclip_1.10-0      
# [37] gtable_0.3.0          leiden_0.3.7          car_3.0-10            future.apply_1.7.0    abind_1.4-5           scales_1.1.1         
# [43] pheatmap_1.0.12       DBI_1.1.1             rstatix_0.7.0         miniUI_0.1.1.1        Rcpp_1.0.6            htmlTable_2.1.0      
# [49] viridisLite_0.3.0     xtable_1.8-4          reticulate_1.18       spatstat.core_2.0-0   foreign_0.8-81        Formula_1.2-4        
# [55] htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2             farver_2.0.3         
# [61] pkgconfig_2.0.3       nnet_7.3-15           uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10         labeling_0.4.2       
# [67] tidyselect_1.1.0      rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0         cellranger_1.1.0     
# [73] tools_4.0.5           cli_2.5.0             generics_0.1.0        broom_0.7.5           ggridges_0.5.3        fastmap_1.1.0        
# [79] goftest_1.2-2         knitr_1.31            fs_1.5.0              fitdistrplus_1.1-3    zip_2.1.1             RANN_2.6.1           
# [85] pbapply_1.4-3         future_1.21.0         nlme_3.1-152          mime_0.10             xml2_1.3.2            compiler_4.0.5       
# [91] rstudioapi_0.13       plotly_4.9.3          curl_4.3              png_0.1-7             ggsignif_0.6.0        spatstat.utils_2.1-0 
# [97] reprex_1.0.0          stringi_1.5.3         RSpectra_0.16-0       lattice_0.20-41       ggsci_2.9             vctrs_0.3.6          
# [103] pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18      cowplot_1.1.1        
# [109] irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       R6_2.5.0              latticeExtra_0.6-29   promises_1.2.0.1     
# [115] KernSmooth_2.23-18    gridExtra_2.3         rio_0.5.16            parallelly_1.23.0     codetools_0.2-18      MASS_7.3-53.1        
# [121] assertthat_0.2.1      rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           parallel_4.0.5       
# [127] hms_1.0.0             grid_4.0.5            rpart_4.1-15          carData_3.0-4         Rtsne_0.15            ggpubr_0.4.0         
# [133] shiny_1.6.0           lubridate_1.7.9.2     base64enc_0.1-3      