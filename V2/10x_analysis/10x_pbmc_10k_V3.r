# 10x analysis on 10K dataset including background sensitivity analysis  
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(dsb))
suppressMessages(library(parallelDist))
suppressMessages(library(here))
set.seed(1990)

# set params 
project_title = "10X PBMC 10K V3"
expected_cells = "~10,000"
res = 0.3
kparam = 40

# thresholds 
cellranger_cells = data.table::fread(input = here("data/10x_rds/v3_10k/filtered_feature_bc_matrix/barcodes.tsv.gz"),header = FALSE)$V1
cellranger_cells = str_remove(string = cellranger_cells,pattern = "-1")
max_neg_logprotumi = 2.8
min_cell_logprotumi = 3
genemax = 3000
genemin = 80
mtthresh = 0.1

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/"); dir.create(figpath, recursive = TRUE)
datapath = here("V2/10x_analysis/generated_data/"); dir.create(datapath)

# functions 
source("V2/functions/preprocessing_functions.R")
source("V2/functions/analysis_functions.R")

# read data 
raw = readRDS(file = here("data/10x_rds/10x_pbmc10k_V3.rds"))
prot = raw[grep(rownames(raw), pattern = "Total"), ]
rna = raw[rownames(raw)[rownames(raw) %ni% rownames(prot)], ]

# calculate QC stats (pct mt is proportion of reads that are mitochondrial)
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE)
pctmt = colSums(rna[mtgene, ])/colSums(rna)
umi  = colSums(rna)
log10umi = log10(umi)
umiprot = colSums(prot)
log10umiprot = log10(umiprot)
nGene = colSums(rna > 0)

# check to see if there are protein detected in drops with no RNA 
pdrop = names(umiprot)[names(umiprot) %ni% names(umi)]
stopifnot(length(pdrop) == 0)

# combine into metadata 
md = as.data.frame(cbind(pctmt, umi, log10umi, nGene, umiprot, log10umiprot))

# define negative background and cells 
##### Add additional filters for sensitivity analysis shown in Supplementary Fig 7 
neg_drops1 = md %>%
  rownames_to_column("bc") %>% 
  filter(bc %ni% cellranger_cells) %>% 
  filter(log10umiprot < max_neg_logprotumi & log10umiprot > 0) %>%
  filter(nGene < genemin) %$% bc
neg_drops2 = md %>% 
  rownames_to_column("bc") %>%
  filter(bc %ni% cellranger_cells) %>% 
  filter(log10umiprot < max_neg_logprotumi & log10umiprot > 2)  %>% 
  filter(nGene < genemin) %$% bc
neg_drops3 = md %>% 
  rownames_to_column("bc") %>%
  filter(bc %ni% cellranger_cells) %>% 
  filter(log10umiprot < 2 & log10umiprot > 0) %>% 
  filter(nGene < genemin) %$% bc
neg_drops4 = md %>% 
  rownames_to_column("bc") %>% 
  filter(bc %ni% cellranger_cells) %>% 
  filter(log10umiprot < 2.5 & log10umiprot > 0.9) %>%
  filter(nGene < genemin) %$% bc

# subset protein matrix to 3 defined empty drop subsets 
neg_prot1 = prot[ , neg_drops1] %>% as.matrix()
neg_prot2 = prot[ , neg_drops2] %>%  as.matrix()
neg_prot3 = prot[ , neg_drops3] %>% as.matrix()
neg_prot4 = prot[ , neg_drops4] %>% as.matrix()

# subset out outlier drops from positive protein matrix RNA based; add absolute ceiling for conservative cell estimate. 
positive_cells = md %>% 
  rownames_to_column("bc") %>% 
  filter(bc %in% cellranger_cells) %>% 
  filter(log10umiprot > min_cell_logprotumi) %>% 
  filter(nGene < genemax & nGene > 200) %>% 
  filter(pctmt < mtthresh) %$% 
  bc

# add cell ranger cell call to metadata 
md$droplet_class = ifelse(rownames(md) %in% cellranger_cells, "cell", "background")

# subset protein matrix for cells 
pos_prot = prot[ , positive_cells] %>% as.matrix()

# define non-staining prot 
ns = data.frame(pmax = apply(pos_prot, 1, max)) %>% 
  rownames_to_column('prot') %>%
  arrange(pmax) %>% 
  filter(pmax < 5) %$% 
  prot 
length(ns)
# [1] 0

# define isotype controls for dsb 
isotypes = rownames(pos_prot)[grepl(pattern = "control", x = rownames(pos_prot))]; isotypes 

###########
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
# T1 
pv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(pos_prot)) %>% mutate(class = "cell_containing")
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot1)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 1 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot1)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "T1protein_joint_lib_distribution.pdf"), width = 4, height = 3.5)

# T2 
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot2)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot2)
  )) + 
  plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "T2protein_joint_lib_distribution.pdf"), width = 4, height = 3.5)

# T3 
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot3)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 3 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot3)
  )) + 
  plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "T3protein_joint_lib_distribution.pdf"), width = 4, height = 3.5)

# T4 
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot4)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 4 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot4)
  )) + 
  plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "T4protein_joint_lib_distribution.pdf"), width = 4, height = 3.5)

# save raw outputs 
saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
saveRDS(neg_prot1,file = paste0(datapath, project_title, "neg_prot1.rds"))
saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))
saveRDS(neg_prot3,file = paste0(datapath, project_title, "neg_prot3.rds"))
saveRDS(neg_prot4,file = paste0(datapath, project_title, "neg_prot4.rds"))

##########################
# dsb normalize with different bacground thresholds for sensitivity analysis 
mtx = list()
nplist = list(neg_prot1, neg_prot2, neg_prot3, neg_prot4)
for (i in 1:length(nplist)) {
  mtx[[i]] = DSBNormalizeProtein(cell_protein_matrix = pos_prot,
                                 empty_drop_matrix = nplist[[i]],
                                 denoise.counts = TRUE,
                                 use.isotype.control = TRUE,
                                 isotype.control.name.vec = isotypes
  )
}

# visualize distributions
mgl = list(theme_bw(), 
           geom_bin2d(bins = 200, show.legend = FALSE),
           scale_fill_viridis_c(option = "B"), 
           geom_vline(xintercept = 0, linetype = 'dashed'), 
           geom_hline(yintercept = 0, linetype = 'dashed'), 
           geom_vline(xintercept = 3.5, color = 'red'), 
           geom_hline(yintercept = 3.5, color = 'red'))
for (i in 1:length(mtx)) {
  dfplot = as.data.frame(t(mtx[[i]]))
  p2 = ggplot(dfplot, aes(x = CD4_TotalSeqB, y = CD14_TotalSeqB)) + mgl
  ggsave(p2, filename = paste0(figpath, project_title, "thres_", i, " CD4_CD14_MG.pdf"), width = 3.5, height = 3.5)
}
  
############# 
# cluster 
rna_cells = rna[ ,positive_cells]
md_cells = md[positive_cells, ]
s = CreateSeuratObject(raw.data = rna_cells, min.cells = 40, min.genes = genemin, meta.data = md_cells)
s = SetAssayData(s, assay.type = "CITE", slot = "data",new.data = mtx[[2]])

##### cluster 
prot = rownames(s@assay$CITE@data)
prot_subset = setdiff(prot, isotypes)

# Subset sd background normalized denoised protein 
s2_adt = GetAssayData(s, assay.type = "CITE", slot = "data")
s2_adt3 = s2_adt[prot_subset, ]
p3_dist = parDist(t(s2_adt3))
p3_dist = as.matrix(p3_dist)

# cluster based on surface protein
s = FindClusters(s, 
                 distance.matrix = p3_dist,
                 k.param = kparam,
                 print.output = F, 
                 resolution = res,
                 random.seed = 1,
                 algorithm = 3,
                 modularity.fxn = 1)
s = StashIdent(s, save.name = "clusters")


# run umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(s2_adt3), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(s@assay$CITE@data)))
saveRDS(df_dsb, file = paste0(datapath,project_title, "dsb_merged_result.RDS"))

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
#   [1] umap_0.2.3.1       reticulate_1.12    viridis_0.5.1      viridisLite_0.3.0  here_0.1           parallelDist_0.2.4 dsb_0.1.0         
# [8] magrittr_2.0.1     forcats_0.4.0      stringr_1.4.0      dplyr_0.8.5        purrr_0.3.3        readr_1.3.1        tidyr_1.0.2       
# [15] tibble_2.1.1       tidyverse_1.2.1    Seurat_2.3.4       Matrix_1.2-15      cowplot_0.9.4      ggplot2_3.1.1     
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1        snow_0.4-3          backports_1.1.4     Hmisc_4.2-0         plyr_1.8.4          igraph_1.2.4.1     
# [7] lazyeval_0.2.2      splines_3.5.3       digest_0.6.25       foreach_1.4.4       htmltools_0.3.6     lars_1.2           
# [13] gdata_2.18.0        fansi_0.4.0         checkmate_1.9.3     cluster_2.0.7-1     mixtools_1.1.0      ROCR_1.0-7         
# [19] limma_3.38.3        modelr_0.1.4        RcppParallel_5.0.0  R.utils_2.8.0       askpass_1.1         colorspace_1.4-1   
# [25] rvest_0.3.4         haven_2.1.0         xfun_0.7            crayon_1.3.4        jsonlite_1.6        survival_2.43-3    
# [31] zoo_1.8-6           iterators_1.0.10    ape_5.3             glue_1.3.1          gtable_0.3.0        kernlab_0.9-27     
# [37] prabclus_2.3-1      DEoptimR_1.0-8      scales_1.0.0        pheatmap_1.0.12     mvtnorm_1.0-10      bibtex_0.4.2       
# [43] Rcpp_1.0.1          metap_1.1           dtw_1.20-1          htmlTable_1.13.1    foreign_0.8-71      bit_1.1-14         
# [49] proxy_0.4-23        mclust_5.4.5        SDMTools_1.1-221.1  Formula_1.2-3       stats4_3.5.3        tsne_0.1-3         
# [55] htmlwidgets_1.3     httr_1.4.0          gplots_3.0.1.1      RColorBrewer_1.1-2  fpc_2.2-1           acepack_1.4.1      
# [61] modeltools_0.2-22   ica_1.0-2           pkgconfig_2.0.2     R.methodsS3_1.7.1   flexmix_2.3-15      nnet_7.3-12        
# [67] utf8_1.1.4          tidyselect_0.2.5    labeling_0.3        rlang_0.4.5         reshape2_1.4.3      munsell_0.5.0      
# [73] cellranger_1.1.0    tools_3.5.3         cli_1.1.0           generics_0.0.2      broom_0.5.2         ggridges_0.5.1     
# [79] npsurv_0.4-0        knitr_1.23          bit64_0.9-7         fitdistrplus_1.0-14 robustbase_0.93-5   caTools_1.17.1.2   
# [85] RANN_2.6.1          packrat_0.5.0       pbapply_1.4-0       nlme_3.1-137        R.oo_1.22.0         xml2_1.2.0         
# [91] hdf5r_1.2.0         compiler_3.5.3      rstudioapi_0.10     png_0.1-7           lsei_1.2-0          stringi_1.4.3      
# [97] RSpectra_0.14-0     lattice_0.20-38     ggsci_2.9           vctrs_0.2.4         pillar_1.4.1        lifecycle_0.1.0    
# [103] Rdpack_0.11-0       lmtest_0.9-37       data.table_1.12.2   bitops_1.0-6        irlba_2.3.3         gbRd_0.4-11        
# [109] R6_2.4.0            latticeExtra_0.6-28 KernSmooth_2.23-15  gridExtra_2.3       codetools_0.2-16    MASS_7.3-51.1      
# [115] gtools_3.8.1        assertthat_0.2.1    openssl_1.4         rprojroot_1.3-2     withr_2.1.2         diptest_0.75-7     
# [121] parallel_3.5.3      doSNOW_1.0.16       hms_0.4.2           grid_3.5.3          rpart_4.1-13        class_7.3-15       
# [127] segmented_0.5-4.0   Rtsne_0.15          lubridate_1.7.4     base64enc_0.1-3    
