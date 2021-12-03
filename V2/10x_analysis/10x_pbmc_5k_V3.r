# 10x analysis 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(dsb))
suppressMessages(library(parallelDist))
suppressMessages(library(here))
set.seed(1990)

# set params 
project_title = "10X PBMC5k V3"
expected_cells = "~5,000"
res = 0.8
kparam = 40

# thresholds 
cellranger_cells = data.table::fread(input = here("data/10x_rds/v3_5k/filtered_feature_bc_matrix/barcodes.tsv.gz"),header = FALSE)$V1
cellranger_cells = str_remove(string = cellranger_cells,pattern = "-1")
max_neg_logprotumi = 2.5
min_cell_logprotumi = 2.8
genemax = 3000
genemin = 80
mtthresh = 0.2

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/"); dir.create(figpath, recursive = TRUE)
datapath = here("V2/10x_analysis/generated_data/"); dir.create(datapath)
source("V2/functions/preprocessing_functions.R")
source("V2/functions/analysis_functions.R")

# read raw data
raw = readRDS(file = here("data/10x_rds/10x_pbmc5k_V3.rds"))

# Separate RNA and protein matrix 
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


# subset matrices 
all.equal(names(umiprot), names(umi))
int = intersect(names(umiprot), names(umi))
bcmax = length(int)
# pctmt = pctmt[int]; umi = umi[int]; log10umi = log10umi[int]; nGene = nGene[int]

# combine into metadata 
md = as.data.frame(cbind(pctmt, umi, log10umi, nGene, umiprot, log10umiprot))

# define negative background subset protein matrix to empty drop subset
neg_drops2 = md %>% 
  rownames_to_column("bc") %>% 
  filter(bc %ni% cellranger_cells) %>% 
  filter(log10umiprot < max_neg_logprotumi & log10umiprot > 1.4)  %>% 
  filter(nGene < genemin) %$% bc
neg_prot2 = prot[ , neg_drops2] %>%  as.matrix()

# subset out outlier drops from positive protein matrix RNA based; add absolute ceiling for conservative cell estimate. 
positive_cells = md %>%
  rownames_to_column("bc") %>% 
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
  filter(pmax < 6) %$% 
  prot 

if (length(ns > 0)) {
  # remove non staining prot
  prot_names = rownames(pos_prot)
  pos_prot = pos_prot[!prot_names == ns, ]
  neg_prot2 = neg_prot2[!prot_names == ns, ]
}
# define isotype controls for dsb 
isotypes = rownames(pos_prot)[29:31]; isotypes 


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
pv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(pos_prot)) %>% mutate(class = "cell_containing")
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot2)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title,"\n",  
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets after QC = ", ncol(neg_prot2)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4, height = 3.5)

# save raw 
saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))

#############
# DSB normalize 
mtx2 = DSBNormalizeProtein(cell_protein_matrix = pos_prot,
                           empty_drop_matrix = neg_prot2,
                           denoise.counts = TRUE,
                           use.isotype.control = TRUE,
                           isotype.control.name.vec = isotypes)

# cluster
rna_cells = rna[ ,positive_cells]
md_cells = md[positive_cells, ]
s = CreateSeuratObject(raw.data = rna_cells, min.cells = 40, min.genes = genemin, meta.data = md_cells)
s = SetAssayData(s, assay.type = "CITE", slot = "data",new.data = mtx2)

##### cluster 
prot = rownames(s@assay$CITE@data)
prot_subset = setdiff(prot, isotypes)

# Subset sd background normalized denoised protein 
s2_adt = GetAssayData(s, assay.type = "CITE", slot = "data")
s2_adt3 = s2_adt[prot_subset, ]
p3_dist = parDist(t(s2_adt3))
p3_dist = as.matrix(p3_dist)

# Protein clustering 
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
#   [1] umap_0.2.3.1       reticulate_1.12    viridis_0.5.1      viridisLite_0.3.0  here_0.1           parallelDist_0.2.4 dsb_0.2.0          magrittr_2.0.1    
# [9] forcats_0.4.0      stringr_1.4.0      dplyr_0.8.5        purrr_0.3.3        readr_1.3.1        tidyr_1.0.2        tibble_2.1.1       tidyverse_1.2.1   
# [17] Seurat_2.3.4       Matrix_1.2-15      cowplot_0.9.4      ggplot2_3.1.1     
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    class_7.3-15        modeltools_0.2-22   ggridges_0.5.1      rprojroot_1.3-2     mclust_5.4.5       
# [8] htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10     proxy_0.4-23        npsurv_0.4-0        flexmix_2.3-15      bit64_0.9-7        
# [15] RSpectra_0.14-0     lubridate_1.7.4     mvtnorm_1.0-10      xml2_1.2.0          codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1  
# [22] lsei_1.2-0          robustbase_0.93-5   knitr_1.23          Formula_1.2-3       jsonlite_1.6        packrat_0.5.0       broom_0.5.2        
# [29] ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      png_0.1-7           R.oo_1.22.0         compiler_3.5.3      httr_1.4.0         
# [36] backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      limma_3.38.3        cli_1.1.0           lars_1.2            acepack_1.4.1      
# [43] htmltools_0.3.6     tools_3.5.3         igraph_1.2.4.1      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          reshape2_1.4.3     
# [50] Rcpp_1.0.1          cellranger_1.1.0    vctrs_0.2.4         gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10   
# [57] fpc_2.2-1           gbRd_0.4-11         lmtest_0.9-37       xfun_0.7            rvest_0.3.4         lifecycle_0.1.0     irlba_2.3.3        
# [64] gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-6           scales_1.0.0        hms_0.4.2           doSNOW_1.0.16      
# [71] parallel_3.5.3      RColorBrewer_1.1-2  pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13        segmented_0.5-4.0   latticeExtra_0.6-28
# [78] stringi_1.4.3       foreach_1.4.4       checkmate_1.9.3     caTools_1.17.1.2    bibtex_0.4.2        Rdpack_0.11-0       SDMTools_1.1-221.1 
# [85] rlang_0.4.5         pkgconfig_2.0.2     dtw_1.20-1          prabclus_2.3-1      bitops_1.0-6        lattice_0.20-38     ROCR_1.0-7         
# [92] labeling_0.3        htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5    ggsci_2.9           plyr_1.8.4          R6_2.4.0           
# [99] generics_0.0.2      snow_0.4-3          gplots_3.0.1.1      Hmisc_4.2-0         haven_2.1.0         pillar_1.4.1        foreign_0.8-71     
# [106] withr_2.1.2         fitdistrplus_1.0-14 mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          modelr_0.1.4       
# [113] crayon_1.3.4        hdf5r_1.2.0         KernSmooth_2.23-15  readxl_1.3.1        grid_3.5.3          data.table_1.12.2   metap_1.1          
# [120] digest_0.6.25       diptest_0.75-7      R.utils_2.8.0       openssl_1.4         RcppParallel_5.0.0  stats4_3.5.3        munsell_0.5.0      
# [127] askpass_1.1     
