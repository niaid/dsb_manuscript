# 10x PBMC 5k v3 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
"%ni%" = Negate("%in%")

# set params 
project_title = "10X PBMC5k V3"
prot_plot = c("CD19_TotalSeqB", "CD3_TotalSeqB", "CD14_TotalSeqB",  "CD4_TotalSeqB" , "CD56_TotalSeqB", "CD8a_TotalSeqB")
neg_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k V3neg_prot2.rds"))
pos_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k V3pos_prot.rds"))
df_dsb = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k V3dsb_merged_result.RDS"))

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/")
datapath = here("V2/10x_analysis/generated_data/")

## make tidy dataframe for plot
index1 = colnames(df_dsb)[14]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# cluster umap plots ; calc centers in umap space 
centers = df_dsb %>% 
  dplyr::group_by(clusters) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# labeled umap 
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +  
  geom_point(mapping = aes(color = clusters), size = 0.5, shape = 16, alpha = 0.3, show.legend = FALSE) + 
  ggsci::scale_color_d3(palette = "category20") + 
  ggrepel::geom_text_repel(data = centers, box.padding = 0.5,
                           mapping = aes(label = clusters, size = 2.3, fontface = "bold"),
                           show.legend = FALSE) 
ggsave(p, filename = paste0(figpath, project_title, "clusters.png"), width = 3.3, height = 3.3)

# protein distributions 
plotsub = dsb %>% filter(DSB > -5) %>% filter(DSB < 40) 
plotsub_spread = plotsub %>% spread(prot, DSB, drop = TRUE)
plotsub_spread =  na.omit(plotsub_spread)

# biaxial plots
mg_theme = list( 
  theme_bw(),
  theme(axis.title.x =element_text(size = 18), axis.title.y = element_text(size = 18)), 
  geom_bin2d(bins = 200, show.legend = FALSE),
  viridis::scale_fill_viridis(option = "B") 
)

# plot manual gate distributions
p = ggplot(plotsub_spread, aes(x = plotsub_spread[[prot_plot[4]]], y = plotsub_spread[[prot_plot[3]]]))  + xlab("CD4") + ylab("CD14") + mg_theme 
ggsave(p, filename = paste0(figpath, project_title, "4_14_mg.pdf"), width = 3, height = 3)  


# protein heatmap by cluster
prots = dsb$prot %>% unique 
adt_plot = plotsub_spread %>% 
  group_by(clusters) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("clusters") %>% 
  t %>% 
  as.data.frame

# average heatmap 
pheatmap::pheatmap(adt_plot, color = viridis::viridis(18, option = "B"), 
                   fontsize_row = 9,border_color = NA, 
                   treeheight_row = 10, treeheight_col = 10,
                   filename = paste0(figpath, project_title, "average_cluster_dsb_heatmap.pdf"),
                   width = 4.5, height = 5)


##########################
## assessment of dsb model 
pseudocount.use = 10
adtu_log = log(neg_cells + pseudocount.use) 
adt_log = log(pos_cells + pseudocount.use)

# apply ambient correction 
mu_u = apply(adtu_log, 1 , mean)
sd_u = apply(adtu_log, 1 , sd)
norm_adt2 = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u) 

# gaussian mixture on background rescaled matrix 
library(mclust)
cellwise_model = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = FALSE)  
  return(g) 
})
stopifnot(isTRUE(all.equal(colnames(adt_log), names(cellwise_model))))

# tidy Gaussian mixture model data 
tm  = lapply(cellwise_model, function(x){broom::tidy(x)[ 2, 5:6]})
tm = do.call(rbind, tm)
tm$barcode_check = colnames(adt_log)
mr  = lapply(cellwise_model, broom::glance)
mr = do.call(rbind, mr)
mr$barcode_check = colnames(adt_log)

# merge model results with metadata 
df_dsb$barcode_check = rownames(df_dsb)
md = full_join(df_dsb, mr, by = "barcode_check")
md = full_join(md, tm, by = "barcode_check")

################
# calculate latent component noise vector 
cellwise_background_mean = lapply(cellwise_model, function(x) {
  x$parameters$mean[1] }) %>%
  unlist(use.names = FALSE)

# define latent noise component 
isotypes = rownames(adt_log)[29:31]; isotypes
noise_matrix = rbind(norm_adt2[isotypes, ], cellwise_background_mean)
get_noise_vector = function(noise_matrix) { 
  g = prcomp(t(noise_matrix), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
md = full_join(md, PC, by = "barcode_check")

# noise vector vs protein library size clusters 
p = ggplot(md, aes(x = log10umiprot, y = -1*(noise_vector))) + 
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
  facet_wrap(~clusters, scales = 'free')
ggsave(p, filename = paste0(figpath, project_title, "ALTBINnoise_vector_vs_libsize.pdf") , width = 4.8, height = 4)


# isotype control vs mixture mean 1 
iso = norm_adt2[isotypes,  ]
iso_mean = apply(iso, 2, mean)

# add isotype mean to metadata.  
md = cbind(md, iso_mean)

# plot 
p = ggplot(data = md, 
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
isotype_value = iso %>% t %>% as.data.frame() %>% rownames_to_column("bc")
mean1df  = md %>% select(mean.1, bc = barcode_check)
dcor = full_join(isotype_value, mean1df, by = "bc") %>% column_to_rownames("bc")
cplot = cor(dcor, method = "spearman")
pcor = Hmisc::rcorr(as.matrix(dcor), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "µ1")
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
data.table::fwrite(md, file = paste0(datapath, project_title, '_cellmetadata.txt'),sep = "\t")


#####################
# mixture model parameter evaluation 
cmd = list()
for (i in 1:3) {
  cmd[[i]]  = apply(norm_adt2, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = FALSE)  
    return(g) ; gc()
  })
}
mr2 = list()
for (i in 1:length(cmd)) {
  mr2[[i]]  = lapply(cmd[[i]], broom::glance)
  mr2[[i]] = do.call(rbind, mr2[[i]])
  mr2[[i]]$barcode_check = colnames(adt_log)
}
mrdf = do.call(rbind, mr2)
cmd = df_dsb %>% select(barcode_check, clusters)
mrdf = full_join(mrdf, cmd, by = "barcode_check")

# save mixture fit dataframe 
write_delim(mrdf, path = paste0(datapath, project_title, "_mixture_model_comparison.txt"),delim = "\t")

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
# 1    2    3 
# 2 3224  255 

# plot g vs BIC 
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3"))
p = ggplot(mrdf, aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,   show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 16)) + 
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title, "GLOBALmixture_model_BICvsG.pdf"), width = 2, height = 3)


###############################
# write stats table for s table.  
d = data.frame(
  dataset = project_title,
  nprot = nrow(pos_cells), 
  ncells = ncol(pos_cells), 
  n_background = ncol(neg_cells), 
  cor_tech_size = cor(md$log10umiprot, -1*(md$noise_vector)), 
  cor_isotype_background = cor(md$iso_mean, md$mean.1), 
  cells_removed_for_visualization = dim(df_dsb)[1] - dim(plotsub_spread)[1]
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
#   [1] mclust_5.4.5    here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2    
# [10] tibble_2.1.1    tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    ellipsis_0.3.0      class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2    
# [8] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        ggpubr_0.2          npsurv_0.4-0        ggrepel_0.8.1      
# [15] flexmix_2.3-15      bit64_0.9-7         fansi_0.4.0         lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16   
# [22] splines_3.5.3       R.methodsS3_1.7.1   lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6       
# [29] packrat_0.5.0       broom_0.5.2         ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0        
# [36] pheatmap_1.0.12     compiler_3.5.3      httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      cli_1.1.0          
# [43] lars_1.2            acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1         
# [50] RANN_2.6.1          reshape2_1.4.3      Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3            
# [57] nlme_3.1-137        iterators_1.0.10    fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4        
# [64] lifecycle_0.1.0     irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0       
# [71] hms_0.4.2           doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.12     pbapply_1.4-0       gridExtra_2.3      
# [78] rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28 stringi_1.4.3       corrplot_0.84       foreach_1.4.4       checkmate_1.9.3    
# [85] caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1  rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1         
# [92] prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7          labeling_0.3        htmlwidgets_1.3     bit_1.1-14         
# [99] tidyselect_0.2.5    ggsci_2.9           plyr_1.8.4          R6_2.4.0            generics_0.0.2      snow_0.4-3          gplots_3.0.1.1     
# [106] Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0     
# [113] survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4        crayon_1.3.4        hdf5r_1.2.0         utf8_1.1.4         
# [120] KernSmooth_2.23-15  viridis_0.5.1       readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1           digest_0.6.25      
# [127] diptest_0.75-7      R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0       viridisLite_0.3.0  
