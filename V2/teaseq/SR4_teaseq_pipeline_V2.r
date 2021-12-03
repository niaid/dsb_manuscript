# R4 Seurat 4 (see session info) 
'%ni%' = Negate('%in%')
suppressMessages(library(dsb))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# project info 
project_title = "TEA-seq"
figpath = paste0(here("V2/teaseq/figures_v2/"), project_title, "/"); dir.create(figpath, recursive = TRUE)
datapath = paste0(here("V2/teaseq/generated_data_v2/"), project_title, "/"); dir.create(datapath, recursive = TRUE)

# load metadata data created from geoQuery from preprocess script 1
cell_meta = readRDS(here("data/revision_data/tea_seq_preprocess/generated_data/filteredmetadata_uniqueIDS_barcode.rds"))
unfiltered_meta = readRDS(here("data/revision_data/tea_seq_preprocess/generated_data/UNfilteredmetadata_uniqueIDS_barcode.rds"))
cell_meta$droplet_class = "cell"
unfiltered_meta$droplet_class = ifelse(unfiltered_meta$barcode %in% cell_meta$barcode, 'cell', 'background')

# read RNA data
in_hdf = list.files(path = here('data/revision_data/tea_seq'),recursive = TRUE, pattern = "matrix.h5", full.names = TRUE)
rna = lapply(in_hdf, function(x){ Seurat::Read10X_h5(filename = x)$`Gene Expression`})

# all barcodes are appended with 1 
lapply(rna, function(x) unique(str_sub(colnames(x), -1,-1)))

# now append correct well number onto barcodes 
well_num <- sub(".+W([0-9])_.+","\\1",in_hdf)
for (i in seq_along(well_num)) {
  colnames(rna[[i]]) = str_replace_all(string = colnames(rna[[i]]), 
                                       pattern = '1',
                                       replacement = as.character(well_num[i]))
}

# check barcodes are unique now that the well number is appended 
lapply(rna, function(x) all(Biobase::isUnique(colnames(x))))

# merge barcodes across lanes and sub to include those in filtered metadata 
rna_raw = do.call(cbind, rna)
rna_raw = rna_raw[ , colnames(rna_raw) %in% cell_meta$barcode ]

# read filtered ADTs from script 1 
cell_adt_data = data.table::fread(file = here('data/revision_data/tea_seq_preprocess/generated_data/Swanson_filtered_adt_counts.csv.gz'))
unfiltered_adt_data = data.table::fread(file = here('data/revision_data/tea_seq_preprocess/generated_data/Swanson_UNfiltered_adt_counts.csv.gz'))
filtered_adt_data = unfiltered_adt_data %>% filter(barcode %in% unfiltered_meta$barcode)

# convert to matrix with the correct barcode column with appended lane as rownames 
adt_ = filtered_adt_data %>%
  select(-c('total', 'cell_barcode')) %>%
  tibble::remove_rownames() %>% 
  column_to_rownames('barcode') %>% 
  t()

# calculate qc stats for teaseq data 
adt_size = colSums(adt_) 
adt_nprot = colSums(adt_ > 0)
adt_size_log10 = log10(adt_size)
d = data.frame(adt_size, adt_size_log10, adt_nprot)
d$barcode = rownames(d)

# add adt meta 
unfiltered_meta = full_join(unfiltered_meta, d, by = "barcode")

# plot drop qc cell containing drops defined by authors based on multiple 
# QC stats including RNA, peaks and ADT. 
p = ggplot(unfiltered_meta %>% 
             filter(adt_size_log10 > 0.2 & gex_genes_count > 10 ),
           aes(x = adt_size_log10, y = log10(gex_raw_reads) )) + 
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  geom_bin2d(bins = 100)  + 
  scale_fill_viridis_c(option = "C") + 
  facet_wrap(~droplet_class) + 
  xlab("log10 protein library size") + 
  ylab("log10 RNA library size")
ggsave(paste0(figpath, 'droplet_qc_teaseq.pdf'), width = 4.5, height = 2.3)

# subset major background signal for dsb normalization
background_drops_dsb = 
  unfiltered_meta %>%
  filter(gex_genes_count > 10 & gex_genes_count < 200) %>% 
  filter(adt_size_log10 > 1.6 & adt_size_log10 < 2.5 & log10(gex_raw_reads) < 3.3  & log10(gex_raw_reads) > 2) %>% 
  filter(barcode %ni% cell_meta$barcode) %$% 
  barcode

# subset background drops and cells from adt_ 
adt_background = adt_[ ,background_drops_dsb] %>% as.matrix()
adt_raw = adt_[ , cell_adt_data$barcode] %>% as.matrix()
adt_raw = adt_raw[ ,colnames(adt_raw) %in% colnames(rna_raw)]

# normalize raw ADT data with dsb 
adt_dsb = dsb::DSBNormalizeProtein(cell_protein_matrix = adt_raw,
                                   empty_drop_matrix = adt_background,
                                   denoise.counts = TRUE,
                                   use.isotype.control = TRUE,
                                   isotype.control.name.vec = "IgG1-K-Isotype-Control") 

stopifnot(isTRUE(
  all.equal(colnames(adt_dsb), colnames(rna_raw))
))

# create object with RNA assay (min cells is gene filter)
s = CreateSeuratObject(counts = rna_raw, assay = "RNA", min.cells = 20)

# ADD CITE data in separate slots for CLR and dsb 
cite_assay = CreateAssayObject(counts = adt_raw)
s[["CITE"]] = cite_assay
s[["CLR"]] = cite_assay

# add dsb normalized data to object 
s = SetAssayData(s,slot = "data", assay = "CITE", new.data = adt_dsb)

# process RNA data with vst residual var genes + pca. 
s = NormalizeData(s) %>%
  FindVariableFeatures(selection.method = 'vst', verbose = TRUE) %>% 
  ScaleData() %>% 
  RunPCA()

#############################################################
# WNN analysis comparison with dsb and CLR transformed counts 

# run dsb workflow using normalized proetin values directly 
DefaultAssay(s) <- 'CITE'
# hack to use normalized protein values as a dimensionality reduction object.
# NOT using these scaled values; must run scale() to prevent function stop()
# NOT using principal components; must run RunPCA() to prevent function stop()
s = ScaleData(s, assay = 'CITE')
VariableFeatures(s) <- setdiff(rownames(adt_dsb), "IgG1-K-Isotype-Control")
s = RunPCA(s, reduction.name = 'pdsb', features = VariableFeatures(s))

# make matrix of normalized ADT to add as 'pseudo' dr embeddings
# to use them directly in the WNN algorithm 
pseudo = t(s@assays$CITE@data)
# subset out isotype control 
pseudo = pseudo[ , -40]
pseudo_colnames = paste('pseudo', 1:45, sep = "_")
colnames(pseudo) = pseudo_colnames

# add to object (as pseudo reduction)
s@reductions$pdsb@cell.embeddings = pseudo

# run WNN 
s = FindMultiModalNeighbors(
  object = s, 
  reduction.list = list("pca", "pdsb"),
  weighted.nn.name = "dsb_pseudo_wnn", 
  knn.graph.name = "dsb_pseudo_knn",
  modality.weight.name = "dsb_pseudo_weight",
  snn.graph.name = "dsb_pseudo_snn",
  dims.list = list(1:30, 1:45)
)

# set specific names for weights for first WNN analysis. 
colnames(s@meta.data)[colnames(s@meta.data) == "RNA.weight"] = 'dsb.RNA.weight'
colnames(s@meta.data)[colnames(s@meta.data) == "CITE.weight"] = 'dsb.CITE.weight'

# cluster on pseudo WNN 
s = FindClusters(s, graph.name = "dsb_pseudo_knn", n.start = 5, n.iter = 5,
                 algorithm = 3, resolution = 2, random.seed = 1990, verbose = FALSE)
s = RunUMAP(s, nn.name = "dsb_pseudo_wnn", reduction.name = "dsb_wnn_umap",
            reduction.key = "dsb_wnnUMAP_", seed.use = 1990)


###################################
## clr comparison workflow 
# the code below runs WNN as above using hte same RNA based PCs but different normalozed protein values as input
# the steps above are exactly the same ad the workflow above 
DefaultAssay(s) <- 'CLR'
s = NormalizeData(s, normalization.method = 'CLR', margin = 2) %>% ScaleData() 
VariableFeatures(s) <- setdiff(rownames(adt_dsb), "IgG1-K-Isotype-Control")
s = RunPCA(s, reduction.name = 'pCLR', features = VariableFeatures(s))
pseudo = t(s@assays$CLR@data)
pseudo = pseudo[ , -40]
pseudo_colnames = paste('pseudo', 1:45, sep = "_")
colnames(pseudo) = pseudo_colnames
s@reductions$pCLR@cell.embeddings = pseudo
# run WNN 
s = FindMultiModalNeighbors(
  s, 
  reduction.list = list("pca", "pCLR"),
  weighted.nn.name = "clr_pseudo_wnn", 
  knn.graph.name = "clr_pseudo_knn",
  modality.weight.name = "clr_pseudo_weight",
  snn.graph.name = "clr_pseudo_snn",
  dims.list = list(1:30, 1:45)
)

# set specific names for weights for second WNN analysis. 
colnames(s@meta.data)[colnames(s@meta.data) == "RNA.weight"] = 'CLR.RNA.weight'
colnames(s@meta.data)[colnames(s@meta.data) == "CLR.weight"] = 'CLR.CLR.weight'

# cluster on pseudo WNN 
s = FindClusters(s, graph.name = "clr_pseudo_knn", n.start = 5, n.iter = 5, 
                 algorithm = 3, resolution = 2, random.seed = 1990,  verbose = FALSE)
s = RunUMAP(s, nn.name = "clr_pseudo_wnn", reduction.name = "CLR_wnn_umap",
            reduction.key = "CLR_wnnUMAP_", seed.use = 1990)

# add cell metadata to object 
cell_meta = cell_meta %>% 
  filter(barcode %in% rownames(s@meta.data)) %>% 
  column_to_rownames('barcode')
s = AddMetaData(s,metadata = cell_meta)

# save data 
saveRDS(adt_background,file = paste0(datapath,'adt_background.rds'))
saveRDS(object = s, file = paste0(datapath, 'full_teaseq_r4s4_object_processed.rds'))

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
#   [1] here_1.0.1         forcats_0.5.1      stringr_1.4.0      dplyr_1.0.4        purrr_0.3.4        readr_1.4.0        tidyr_1.1.2       
# [8] tibble_3.0.6       ggplot2_3.3.3      tidyverse_1.3.0    SeuratObject_4.0.0 Seurat_4.0.1       magrittr_2.0.1     dsb_0.1.0         
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1        ggridges_0.5.3        rprojroot_2.0.2      
# [7] mclust_5.4.7          fs_1.5.0              rstudioapi_0.13       spatstat.data_2.1-0   farver_2.0.3          leiden_0.3.7         
# [13] listenv_0.8.0         bit64_4.0.5           ggrepel_0.9.1         RSpectra_0.16-0       lubridate_1.7.9.2     xml2_1.3.2           
# [19] R.methodsS3_1.8.1     codetools_0.2-18      splines_4.0.5         polyclip_1.10-0       jsonlite_1.7.2        broom_0.7.5          
# [25] ica_1.0-2             cluster_2.1.2         dbplyr_2.1.0          R.oo_1.24.0           png_0.1-7             uwot_0.1.10          
# [31] shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.5        httr_1.4.2            backports_1.2.1      
# [37] assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2        cli_2.5.0             limma_3.46.0         
# [43] later_1.1.0.1         htmltools_0.5.1.1     tools_4.0.5           igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
# [49] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6            Biobase_2.50.0        scattermore_0.7       cellranger_1.1.0     
# [55] vctrs_0.3.6           nlme_3.1-152          lmtest_0.9-38         globals_0.14.0        rvest_0.3.6           mime_0.10            
# [61] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2         future_1.21.0         MASS_7.3-53.1        
# [67] zoo_1.8-8             scales_1.1.1          spatstat.core_2.0-0   hms_1.0.0             promises_1.2.0.1      spatstat.utils_2.1-0 
# [73] parallel_4.0.5        RColorBrewer_1.1-2    reticulate_1.18       pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15         
# [79] stringi_1.5.3         BiocGenerics_0.36.1   rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41      
# [85] ROCR_1.0-11           tensor_1.5            labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.3     bit_4.0.4            
# [91] cowplot_1.1.1         tidyselect_1.1.0      parallelly_1.23.0     RcppAnnoy_0.0.18      plyr_1.8.6            R6_2.5.0             
# [97] generics_0.1.0        DBI_1.1.1             withr_2.4.1           pillar_1.4.7          haven_2.3.1           mgcv_1.8-34          
# [103] fitdistrplus_1.1-3    survival_3.2-10       abind_1.4-5           future.apply_1.7.0    hdf5r_1.3.3           modelr_0.1.8         
# [109] crayon_1.4.1          KernSmooth_2.23-18    spatstat.geom_2.0-1   plotly_4.9.3          grid_4.0.5            readxl_1.3.1         
# [115] data.table_1.14.0     reprex_1.0.0          digest_0.6.27         xtable_1.8-4          httpuv_1.5.5          R.utils_2.10.1       
# [121] munsell_0.5.0         viridisLite_0.3.0 
