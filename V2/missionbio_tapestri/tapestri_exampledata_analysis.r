set.seed(1)
suppressMessages(library(tidyverse))
suppressMessages(library(dsb))
suppressMessages(library(Seurat))
suppressMessages(library(here))

# proj title for fig paths 
project_title = "missionbio tapestri"

## data paths 
figpath = here("V2/missionbio_tapestri/figures/"); dir.create(figpath)
datapath = here("V2/missionbio_tapestri/generated_data/"); dir.create(datapath)
test = read.delim(file = here("data/mission_bio_data/AML-4-cell-line-multiomics-adt-counts.tsv"), sep = "\t",header = T)

# make dataframe
test = as.data.frame(test)
prot = test %>% spread(ab_description, raw) 
prot[is.na(prot)] =  0 

# transpose into cells x prot matrix 
prot = prot %>% 
  column_to_rownames("cell_barcode") %>% 
  t %>% 
  as.data.frame()

# calculate library size of droplets to make rough thresholds for cell containing and ambient droplets 
prot_size = log10(Matrix::colSums(prot, na.rm = TRUE))
md = as.data.frame(prot_size)
md$bc = colnames(prot)
hist(md$prot_size, breaks = 100)

# define a vector of background / empty droplet barcodes based on protein library size
background_drops = md[md$prot_size < 2.5 & md$prot_size > 1.4, ]$bc
negative_mtx_rawprot = prot[ , background_drops] %>%  as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size 
positive_cells = md[md$prot_size > 2.7, ]$bc
cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()

# no isotype data available 
# normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot,
  empty_drop_matrix = negative_mtx_rawprot,
  denoise.counts = FALSE,
  use.isotype.control = FALSE
)

##########################
# visualization of drop distribution 
plot_layer = list(theme_bw() , 
                  ggsci::scale_fill_d3(), ggsci::scale_color_d3() ,
                  geom_histogram(aes(y=..count..), alpha=0.5, bins = 50,position="identity"),
                  geom_density(alpha = 0.5), 
                  ylab("Number of Drops"),  xlab("log10 protein library size"), 
                  theme(axis.title.x = element_text(size = 14)),
                  theme(plot.title = element_text(size = 14)),
                  theme(legend.position = c(0.8, 0.7), legend.margin = margin(0,0,0,0))
)
pv = md  %>% filter(bc %in% colnames(cells_mtx_rawprot)) %>% mutate(class = "cell_containing")
nv = md  %>% filter(bc %in% colnames(negative_mtx_rawprot)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = prot_size, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", nrow(test), "\n", 
    "estimated cell containing drops  = ", ncol(cells_mtx_rawprot), "\n",
    "negative droplets = ", ncol(negative_mtx_rawprot)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = prot_size, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 4.5)


##########################
# clustering analysis 

# calculate distance matrix for clustering 
p_dist = dist(t(dsb_norm_prot))
p_dist = as.matrix(p_dist)

# Graph based clustering 
s = CreateSeuratObject(raw.data =  cells_mtx_rawprot, min.cells = 0, min.genes = 0)
s = FindClusters(s, resolution = 0.5, distance.matrix = p_dist)


# heatmap of protein values 
prots = rownames(s@raw.data)
adt_data = cbind(s@meta.data, as.data.frame(t(dsb_norm_prot)))
adt_plot = adt_data %>% 
  group_by(res.0.5) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("res.0.5") %>% 
  t %>% 
  as.data.frame

x = pheatmap::pheatmap(adt_plot, color = viridis::viridis(13, option = "B"), treeheight_row = 10, treeheight_col = 10,
                   filename = paste0(figpath,project_title, "_tapestri_heatmap.pdf"),
                   fontsize_row = 12, fontsize_col = 12, border_color = NA, width = 4, height = 4)

# umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(dsb_norm_prot), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
stopifnot(isTRUE(all.equal(rownames(s@meta.data), rownames(umap_res))))
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(dsb_norm_prot)))
saveRDS(df_dsb, file = paste0(datapath, project_title,"df_dsb.rds"))


# clusters by umap dims 
centers = df_dsb %>% 
  dplyr::group_by(res.0.5) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
p = ggplot(df_dsb, aes(x = UMAP_1, y= UMAP_2)) +
  theme_void() +
  theme(plot.background = element_rect(colour = "black", size = 1)) + 
  geom_point(mapping = aes(color = res.0.5), size = 0.7, shape = 16, alpha = 0.8, show.legend = FALSE) + 
  ggsci::scale_color_d3(palette = "category20") + 
  ggrepel::geom_text_repel(data = centers, box.padding = 0.5,
                           mapping = aes(label = res.0.5, size = 2.3, fontface = "bold"),
                           show.legend = FALSE) 
ggsave(p, filename = paste0(figpath,"tapestri_umap_dsb.png"), width = 3, height = 3)



# table stats 
ncol(cells_mtx_rawprot)
# [1] 1049
ncol(negative_mtx_rawprot)
# [1] 16455

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
#   [1] here_0.1        Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   dsb_0.1.0       forcats_0.4.0   stringr_1.4.0  
# [8] dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_2.1.1    ggplot2_3.1.1   tidyverse_1.2.1
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1        snow_0.4-3          backports_1.1.4     Hmisc_4.2-0         plyr_1.8.4          igraph_1.2.4.1     
# [7] lazyeval_0.2.2      splines_3.5.3       digest_0.6.25       foreach_1.4.4       htmltools_0.3.6     viridis_0.5.1      
# [13] lars_1.2            gdata_2.18.0        magrittr_2.0.1      checkmate_1.9.3     cluster_2.0.7-1     mixtools_1.1.0     
# [19] ROCR_1.0-7          limma_3.38.3        modelr_0.1.4        R.utils_2.8.0       colorspace_1.4-1    ggrepel_0.8.1      
# [25] rvest_0.3.4         haven_2.1.0         xfun_0.7            crayon_1.3.4        jsonlite_1.6        survival_2.43-3    
# [31] zoo_1.8-6           iterators_1.0.10    ape_5.3             glue_1.3.1          gtable_0.3.0        kernlab_0.9-27     
# [37] prabclus_2.3-1      DEoptimR_1.0-8      scales_1.0.0        pheatmap_1.0.12     mvtnorm_1.0-10      bibtex_0.4.2       
# [43] Rcpp_1.0.1          metap_1.1           dtw_1.20-1          viridisLite_0.3.0   htmlTable_1.13.1    reticulate_1.12    
# [49] foreign_0.8-71      bit_1.1-14          proxy_0.4-23        mclust_5.4.5        SDMTools_1.1-221.1  Formula_1.2-3      
# [55] stats4_3.5.3        tsne_0.1-3          htmlwidgets_1.3     httr_1.4.0          gplots_3.0.1.1      RColorBrewer_1.1-2 
# [61] fpc_2.2-1           acepack_1.4.1       modeltools_0.2-22   ica_1.0-2           pkgconfig_2.0.2     R.methodsS3_1.7.1  
# [67] flexmix_2.3-15      nnet_7.3-12         tidyselect_0.2.5    rlang_0.4.5         reshape2_1.4.3      munsell_0.5.0      
# [73] cellranger_1.1.0    tools_3.5.3         cli_1.1.0           generics_0.0.2      broom_0.5.2         ggridges_0.5.1     
# [79] npsurv_0.4-0        knitr_1.23          bit64_0.9-7         fitdistrplus_1.0-14 robustbase_0.93-5   caTools_1.17.1.2   
# [85] RANN_2.6.1          packrat_0.5.0       pbapply_1.4-0       nlme_3.1-137        R.oo_1.22.0         xml2_1.2.0         
# [91] hdf5r_1.2.0         compiler_3.5.3      rstudioapi_0.10     png_0.1-7           lsei_1.2-0          stringi_1.4.3      
# [97] lattice_0.20-38     ggsci_2.9           vctrs_0.2.4         pillar_1.4.1        lifecycle_0.1.0     Rdpack_0.11-0      
# [103] lmtest_0.9-37       data.table_1.12.2   bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11         R6_2.4.0           
# [109] latticeExtra_0.6-28 KernSmooth_2.23-15  gridExtra_2.3       codetools_0.2-16    MASS_7.3-51.1       gtools_3.8.1       
# [115] assertthat_0.2.1    rprojroot_1.3-2     withr_2.1.2         diptest_0.75-7      parallel_3.5.3      doSNOW_1.0.16      
# [121] hms_0.4.2           grid_3.5.3          rpart_4.1-13        class_7.3-15        segmented_0.5-4.0   Rtsne_0.15         
# [127] lubridate_1.7.4     base64enc_0.1-3 
