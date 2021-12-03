# 10x PBMC 5k Nextgem data 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
"%ni%" = Negate("%in%")

# set params 
project_title = "10X PBMC5k NextGem"
prot_plot = c("CD19_TotalSeqB", "CD3_TotalSeqB", "CD14_TotalSeqB",  "CD4_TotalSeqB" , "CD56_TotalSeqB", "CD8a_TotalSeqB")
neg_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k NextGemneg_prot2.rds"))
pos_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k NextGempos_prot.rds"))
df_dsb = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC5k NextGemdsb_merged_result.RDS"))

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/")
datapath = here("V2/10x_analysis/generated_data/")

# annotate
celltypes = c('Monocytes and mDC', 'Memory CD4 T cell', 'Naive CD4 T cell', 'NK cells', 'Bcell',
              'CD62L+ Naive CD4 T cell', 'Memory CD8 T cell', 'Naive CD8 T cell', 'Double Negative T cell')
clusters_ = 0:8
df_dsb$celltype = plyr::mapvalues(x = df_dsb$clusters, from = clusters_, to = celltypes)
df_dsb = df_dsb %>% select(celltype, everything())

## make tidy dataframe for plot
index1 = colnames(df_dsb)[15]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# cluster umap plots ; calc centers in umap space 
centers = df_dsb %>% 
  dplyr::group_by(celltype) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# labeled umap 
cu = pals::cols25(n = length(unique(dsb$celltype))) %>% unlist(use.names = FALSE)
cu = pals::kelly(n = 15) %>% unlist(use.names = FALSE)
cu = unname(cu)
cu = cu[4:12]
borders = list(theme_bw(),
               theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                     axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(), axis.text.y = element_blank())
)
dsb$clus_celltype = paste(dsb$celltype, dsb$clusters, sep = ': ')
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2, color = clus_celltype)) + 
  borders + 
  geom_point(size = 2, shape = 16, alpha = 0.03, show.legend = FALSE) + 
  scale_color_manual(values = cu) 
ggsave(p, filename = paste0(figpath, project_title, "clusters_3.png"), width = 4.5, height = 4.5)

# color map 
my_hist = ggplot(dsb, aes(UMAP_1, fill = clus_celltype)) + 
  geom_bar() + 
  scale_fill_manual(values = cu) + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = 'right') + 
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5,label.position = "left"))  + 
  guides(fill = guide_legend(label.position = "left", label.hjust = 1))
legend = cowplot::get_legend(my_hist)
grid::grid.newpage()
pdf(file = paste0(figpath, 'plotlegen_celltype.pdf'), width = 4, height = 3)
grid::grid.draw(legend)
dev.off()

## next gem improved plots data save 
saveRDS(dsb, file = paste0(datapath,'dsb_fors4.rds'))

# protein distributions 
plotsub = dsb %>% filter(DSB > -10) %>% filter(DSB < 50) 
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
p = ggplot(plotsub_spread, aes(x = plotsub_spread[[prot_plot[4]]], y = plotsub_spread[[prot_plot[3]]])) + 
  xlab("CD4") + ylab("CD14") + mg_theme +
  theme(axis.title.x = element_text(size = 17)) + 
  theme(axis.title.y = element_text(size = 17))
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


# second version 
# average heatmap 
adt_plot2 = adt_plot
rownames(adt_plot2) = str_replace_all(string = rownames(adt_plot2), pattern = '_TotalSeqB', replacement = "")
pheatmap::pheatmap(adt_plot2, color = viridis::viridis(12, option = "B"), 
                   fontsize_row = 9,
                   border_color = NA, 
                   treeheight_row = 0, treeheight_col = 0,
                   filename = paste0(figpath, project_title, "average_cluster_dsb_heatmap2.pdf"),
                   width = 3.5, height = 4)



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

# shorten names 
md2 = md
celltype2 = c('Monoc and mDC', 'Memory CD4', 'Naive CD4', 'NK cells', 'B cell',
              'CD62L+ CD4', 'Memory CD8', 'Naive CD8', 'Double Neg. T')
md2$celltype2 = plyr::mapvalues(x = md2$celltype, from = celltypes,to = celltype2)

# noise vector vs protein library size clusters 
p = ggplot(md2, aes(x = log10umiprot, y = -1*(noise_vector))) + 
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', aes(label = ..r.label..), size = 5) + 
  geom_bin2d(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(strip.background = element_blank()) + 
  geom_smooth(color = "#3e8ede", method = 'lm') + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.text.y = element_text(size = 5)) + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) + 
  facet_wrap(~celltype2, scales = 'free', nrow = 3) + 
  theme(strip.text = element_text(size = 11))
ggsave(p, filename = paste0(figpath, project_title, "ALTBINnoise_vector_vs_libsize.pdf") , width = 4.4, height = 4)


# isotype control vs mixture mean 1 
iso = norm_adt2[isotypes,  ]
iso_mean = apply(iso, 2, mean)

# add isotype mean to metadata.  
md = cbind(md, iso_mean)

# plot 
p = ggplot(data = md, 
           aes(x = iso_mean, y = mean.1)) + 
  geom_bin2d(bins = 180, show.legend = FALSE) + 
  ggpubr::stat_cor(method = 'pearson', aes(label = ..r.label..), size = 5) + 
  scale_fill_viridis_c(option = "B") + 
  theme_bw() + 
  xlab("isotype means") + 
  ylab(" µ1 background mean ") + 
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) 
ggsave(p, filename = paste0(figpath, project_title,"meanisovs_meanbackground.pdf"), width = 2.5, height = 3)

# correlation across latent component. 
isotype_value = iso %>% t %>% as.data.frame() %>% rownames_to_column("bc")
mean1df  = md %>% select(mean.1, bc = barcode_check)
dcor = full_join(isotype_value, mean1df, by = "bc") %>% column_to_rownames("bc")
#dcor = dcor[complete.cases(dcor), ]
cplot = cor(dcor, method = "spearman")
pcor = Hmisc::rcorr(as.matrix(dcor), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "µ1")
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
# plot 
pdf(file = paste0(figpath,project_title,"background_correlation_plot_spearman.pdf"), width = 3.3,height = 3.3)
corrplot::corrplot(cplot,method="color", col=col(200),  
                   cl.lim = c(0,1),cl.cex = 0.6,
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

########

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
# 7 3485  282 


# plot g vs BIC 
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3"))
p = ggplot(mrdf, aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,   show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 6)) + 
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
#   [1] readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             corrplot_0.84           plyr_1.8.4             
# [7] igraph_1.2.4.1          lazyeval_0.2.2          splines_3.5.3           crosstalk_1.0.0         digest_0.6.25           foreach_1.4.4          
# [13] htmltools_0.3.6         viridis_0.5.1           lars_1.2                gdata_2.18.0            checkmate_1.9.3         cluster_2.0.7-1        
# [19] mixtools_1.1.0          ROCR_1.0-7              modelr_0.1.4            R.utils_2.8.0           colorspace_1.4-1        ggrepel_0.8.1          
# [25] rvest_0.3.4             haven_2.1.0             xfun_0.7                crayon_1.3.4            jsonlite_1.6            survival_2.43-3        
# [31] zoo_1.8-6               iterators_1.0.10        ape_5.3                 glue_1.3.1              pals_1.5                gtable_0.3.0           
# [37] webshot_0.5.1           kernlab_0.9-27          prabclus_2.3-1          DEoptimR_1.0-8          maps_3.3.0              scales_1.0.0           
# [43] pheatmap_1.0.12         mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1          Rcpp_1.0.1              metap_1.1              
# [49] dtw_1.20-1              viridisLite_0.3.0       xtable_1.8-4            htmlTable_1.13.1        reticulate_1.12         foreign_0.8-71         
# [55] bit_1.1-14              mapproj_1.2.6           proxy_0.4-23            SDMTools_1.1-221.1      Formula_1.2-3           stats4_3.5.3           
# [61] tsne_0.1-3              htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2      fpc_2.2-1              
# [67] acepack_1.4.1           modeltools_0.2-22       ellipsis_0.3.0          ica_1.0-2               pkgconfig_2.0.2         R.methodsS3_1.7.1      
# [73] flexmix_2.3-15          nnet_7.3-12             labeling_0.3            tidyselect_0.2.5        rlang_0.4.5             manipulateWidget_0.10.0
# [79] reshape2_1.4.3          later_0.8.0             munsell_0.5.0           cellranger_1.1.0        tools_3.5.3             cli_1.1.0              
# [85] generics_0.0.2          broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0            knitr_1.23              bit64_0.9-7            
# [91] fitdistrplus_1.0-14     robustbase_0.93-5       rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1              packrat_0.5.0          
# [97] pbapply_1.4-0           nlme_3.1-137            mime_0.6                R.oo_1.22.0             xml2_1.2.0              hdf5r_1.2.0            
# [103] compiler_3.5.3          rstudioapi_0.10         png_0.1-7               lsei_1.2-0              stringi_1.4.3           lattice_0.20-38        
# [109] ggsci_2.9               vctrs_0.2.4             pillar_1.4.1            lifecycle_0.1.0         Rdpack_0.11-0           lmtest_0.9-37          
# [115] data.table_1.12.2       bitops_1.0-6            irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0               
# [121] latticeExtra_0.6-28     promises_1.0.1          KernSmooth_2.23-15      gridExtra_2.3           codetools_0.2-16        dichromat_2.0-0        
# [127] MASS_7.3-51.1           gtools_3.8.1            assertthat_0.2.1        rprojroot_1.3-2         withr_2.1.2             diptest_0.75-7         
# [133] parallel_3.5.3          doSNOW_1.0.16           hms_0.4.2               grid_3.5.3              rpart_4.1-13            class_7.3-15           
# [139] segmented_0.5-4.0       Rtsne_0.15              ggpubr_0.2              shiny_1.3.2             lubridate_1.7.4         base64enc_0.1-3       
