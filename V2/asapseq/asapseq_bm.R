# asap seq data 
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(here))
suppressMessages(library(dsb))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

# set save path 
figpath = here('V2/asapseq/figures/'); dir.create(figpath)
datapath = here('V2/asapseq/generated_data/'); dir.create(datapath)

# project 
project_title = "ASAPseq"

# Import ADT [function from https://github.com/caleblareau/asap_reproducibility]
import_kite_counts_bm <- function(){
  # set read dir to /data 
  mtx <- fread(paste0("data/revision_data/asapseq/adt/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(here("data/revision_data/asapseq/adt/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("data/revision_data/asapseq/adt/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}
prot <- import_kite_counts_bm()

# define singlets 
cells <- fread(here("data/revision_data/asapseq/barcodes/step3_ADThq.tsv"))[[1]]
cells = str_replace_all(string = cells, pattern = "-1", replacement = "")

# calculate metadata 
log10_prot_size = log10(Matrix::colSums(prot))
prot_size = Matrix::colSums(prot)
nprot = Matrix::colSums(prot > 0)
md = data.frame(log10_prot_size, prot_size, nprot)
md$bc = rownames(md)
md = md %>% filter(nprot > 10)
md$droplet_class = ifelse(test = md$bc %in% cells, yes = 'cell', no = 'background')

# visualize drop metadata 
p = ggplot(md, aes(x = nprot, y = log10(prot_size) )) +
  theme_bw() + 
  geom_bin2d(bins = 100) + 
  theme(strip.background = element_blank()) + 
  scale_fill_viridis_c(option = "C") + 
  facet_wrap(~droplet_class) 
ggsave(p, filename = paste0(figpath,'dropmd.pdf'), width = 5, height = 3)

# define cell and background matrices 
cells_mtx_rawprot = as.matrix(prot[ , cells])
dim(cells_mtx_rawprot); Matrix::nnzero(cells_mtx_rawprot); 
#1847556 / (242*10927) # 0.6986848; data are still 70% non-zero 

# define background 
background_drops = md[md$log10_prot_size > 2 & 
                        md$log10_prot_size < 3 & !md$droplet_class =='cell', ]$bc
negative_mtx_rawprot = as.matrix(prot[ , background_drops])


### Plot droplet distribution 
# plot background distributions 
plot_layer = list(theme_bw() , 
                  ggsci::scale_fill_d3(), ggsci::scale_color_d3() ,
                  geom_histogram(aes(y=..count..), alpha=0.5, bins = 50,position="identity"),
                  geom_density(alpha = 0.5), 
                  ylab("Number of Drops"),  xlab("log10 protein library size"), 
                  theme(axis.title.x = element_text(size = 14)),
                  theme(plot.title = element_text(size = 10)),
                  theme(legend.position = c(0.8, 0.7), legend.margin = margin(0,0,0,0))
)
pv = md %>%filter(bc %in% cells) %>% mutate(class = "cell_containing")
nv = md %>% filter(bc %in% background_drops) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10_prot_size, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, 
    "\n(Protein + chromatin accessibility)", 
    "\ncell containing drops after QC = ", nrow(pv), "\n",
    "negative droplets = ", nrow(nv)
  )) + plot_layer
xtop = cowplot::axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10_prot_size, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = cowplot::insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = cowplot::ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "_joint_lib_distribution.pdf"), width = 4, height = 3.5)


# define isotypes 
isotypes = rownames(cells_mtx_rawprot)[grepl(x = rownames(cells_mtx_rawprot),pattern =  "sotype")]
prots = setdiff(rownames(cells_mtx_rawprot), isotypes)

# normalize adt with dsb using package defaults 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, 
  empty_drop_matrix = negative_mtx_rawprot, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotypes
)

# create cell data 
cells_cite = prot[ ,cells]

# define md for just cells 
cellmd = md[md$droplet_class == "cell", ]

# Initialize Seurat object w **pseudo** rna asssay
s = Seurat::CreateSeuratObject(counts = cells_cite, meta.data = cellmd, assay = "RNA", min.cells = 0)

# add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
s[["CITE"]] = Seurat::CreateAssayObject(data = dsb_norm_prot)
# add CLR 
s[["CLR"]] = Seurat::CreateAssayObject(counts = cells_mtx_rawprot)
s = NormalizeData(object = s,assay = "CLR",normalization.method = "CLR", margin = 2)

# Select proteins for clustering 
d1 = data.frame(pmax = apply(cells_mtx_rawprot, 1, max)) %>% 
  rownames_to_column('prot') %>%
  arrange(pmax) 
prot_sub = d1 %>% filter(pmax > 12) %$% prot 
prot_sub = str_replace_all(string = prot_sub,pattern = "_",replacement = "-")

#### cluster on dsb normalized values 
s = FindNeighbors(object = s, dims = NULL, assay = 'CITE',features = prot_sub, 
                  k.param = 40, verbose = FALSE)
s = FindClusters(object = s, resolution = 1.5, algorithm = 3, 
                 graph.name = 'CITE_snn',verbose = FALSE)
s = RunUMAP(object = s, assay = "CITE",features = prot_sub,
            seed.use = 1990,reduction.name = "umapdsb",
            min.dist = 0.3, n.neighbors = 40, verbose = FALSE)

# save result data 
df_dsb = cbind(s@meta.data, 
               s@reductions$umapdsb@cell.embeddings, 
               as.data.frame(t(s@assays$CITE@data)))
saveRDS(df_dsb, file = paste0(datapath,project_title, "dsb_merged_result.RDS"))
saveRDS(s,file = paste0(datapath,'s_asapseqbm_processed_protein_seurat4.rds'))
saveRDS(cells_mtx_rawprot, file = paste0(datapath, 'cells_mtx_rawprot.rds'))
saveRDS(negative_mtx_rawprot, file = paste0(datapath, 'negative_mtx_rawprot.rds'))

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
#   [1] magrittr_2.0.1     forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
# [8] tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    dsb_0.1.0          here_1.0.1         Matrix_1.3-2       data.table_1.14.0 
# [15] SeuratObject_4.0.0 Seurat_4.0.1      
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
# [49] viridisLite_0.3.0     xtable_1.8-4          reticulate_1.18       spatstat.core_2.0-0   foreign_0.8-81        mclust_5.4.7         
# [55] Formula_1.2-4         htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2            
# [61] farver_2.0.3          pkgconfig_2.0.3       nnet_7.3-15           uwot_0.1.10           dbplyr_2.1.0          deldir_0.2-10        
# [67] labeling_0.4.2        tidyselect_1.1.0      rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1         munsell_0.5.0        
# [73] cellranger_1.1.0      tools_4.0.5           cli_2.5.0             generics_0.1.0        broom_0.7.5           ggridges_0.5.3       
# [79] fastmap_1.1.0         goftest_1.2-2         knitr_1.31            fs_1.5.0              fitdistrplus_1.1-3    zip_2.1.1            
# [85] RANN_2.6.1            pbapply_1.4-3         future_1.21.0         nlme_3.1-152          mime_0.10             xml2_1.3.2           
# [91] compiler_4.0.5        rstudioapi_0.13       plotly_4.9.3          curl_4.3              png_0.1-7             ggsignif_0.6.0       
# [97] spatstat.utils_2.1-0  reprex_1.0.0          stringi_1.5.3         RSpectra_0.16-0       lattice_0.20-41       ggsci_2.9            
# [103] vctrs_0.3.6           pillar_1.4.7          lifecycle_1.0.0       spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18     
# [109] cowplot_1.1.1         irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       R6_2.5.0              latticeExtra_0.6-29  
# [115] promises_1.2.0.1      KernSmooth_2.23-18    gridExtra_2.3         rio_0.5.16            parallelly_1.23.0     codetools_0.2-18     
# [121] MASS_7.3-53.1         assertthat_0.2.1      rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34          
# [127] parallel_4.0.5        hms_1.0.0             grid_4.0.5            rpart_4.1-15          carData_3.0-4         Rtsne_0.15           
# [133] ggpubr_0.4.0          shiny_1.6.0           lubridate_1.7.9.2     base64enc_0.1-3 
