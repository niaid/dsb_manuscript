Normalizing and denoising protein expression data from droplet-based
single cell profiling
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

### Analysis code to reproduce manuscript results

Mulè MP\* , Martins AJ\* , Tsang JS. [**Normalizing and denoising
protein expression data from droplet based single cell
profiling**](https://www.biorxiv.org/content/10.1101/2020.02.24.963603v3).

This analysis pipeline reproduces all results and figures reported in
the paper above describing control experiments and statistical modeling
used to characterize sources of noise protein expression data from
droplet based single cell experiments (CITE-seq, ASAP-seq, TEA-seq,
REAP-seq, Mission Bio Tapestri etc.). Based on these analysis we
introduced the dsb R package for normalizing and denoising droplet-based
surface protein data. The dsb method was developed in [**John Tsang’s
Lab**](https://www.niaid.nih.gov/research/john-tsang-phd) by Matt Mulè,
Andrew Martins and John Tsang.

The R package dsb is hosted on CRAN  
[**link to latest dsb release on
CRAN**](https://cran.r-project.org/web/packages/dsb/index.html)  
[**link to dsb Github repository**](https://github.com/niaid/dsb)  
**Data used below is available in analysis ready format at this figshare
repository:** **<https://doi.org/10.35092/yhjc.13370915>**

## Table of Contents

1.  [Instructions for analysis workflow](#instructions)
2.  [install packages used in analysis](#software)
3.  [Download starting data and add to data directory](#data)
4.  [dsb normalization on PBMC data from 20 donors](#dsb_1)
5.  [CLR normalization (across cells) on PBMC data from 20 donors
    **R4.0.5**](#clr_cells)
6.  [UMAP based on dsb normalized values](#umap)
7.  [dsb vs CLR nk cluster comparison, manual gating, stained vs
    unstained distribution normalization comparison](#pbmc_analysis)
8.  [dsb technical component robusness assessments, background noise
    comparison and robustness check](#robustness)
9.  [Multi vs single batch normalization, µ1 background resampling
    robustness check](#multibatch)
10. [External 10X genomics data analysis: “NextGem”, “V3”, and “5 Prime”
    assays](#tenx)
11. [dsb normalize protein data from Mission Bio tapestri
    platform](#missionbio)
12. [dsb normalization of TEA-seq data and dsb-based WNN multimodal
    clustering **R4.0.5**](#teaseq)
13. [dsb normalization of ASAP-seq data and dsb-based WNN multimodal
    clustering **R4.0.5**](#asapseq)
14. [dsb vs CLR (acorss cells) Normalization comparison: Differential
    expression, Variance Paratition, Gap Statistic **R4.0.5**](#compare)
15. [dsb vs CLR normalized values as input to WNN multimodal clustering:
    PBMC data from 20 donors **R4.0.5**](#wnn)
16. [dataset summary statistic table](#summarytable)

### Instructions for analysis workflow. <a name="instructions"></a>

All analysis included in this manuscript was run on a laptop with 16GB
RAM. To run the analysis, 1) download the dsb\_manuscript github
repository  
2\) download data from the figshare link above and add the `/data`
folder directly to the root directory containing the .Rproj file; this
directory should now contain dsb\_manuscript.Rproj, the files readme.md
and readme.rmd, the directory `V2` and the directory you just added,
`data`.

One can view the commented code and run each script in each subdirectory
as listed below, or source each R script in the order they appear below.
No file paths need to be specified or changed. Each R script is
self-contained, reading data from the /data folder and writing to
figures or results files within each analysis subdirectory relative to
the root directory using the R package `here`.

Note:  
There are 2 R versions used throughout analysis, R 3.5.3 and R 4.0.5;
analysis using R 4.0.5 are indicated in the table of contents above.
Switching R versions can be done on most HPC systems through separately
installed modules, the [R Switch tool](https://rud.is/rswitch/guide/)
was used throughout analysis to facilitate version switching.

### install packages used in analysis <a name="software"></a>

``` r
# 
pkgs3.5 = list('tidyverse', 'magrittr', 'mclust', 'Seurat', 
               'ggrepel', 'ggridges', 'pals', 'reticulate', 'umap')

lapply(pkgs3.5, function(x){ 
  tryCatch(library(x), 
         error = function(e){
           install.packages(pkgs =  x, repos = 'http://cran.us.r-project.org')
           library(x)
         })
  })

# Note, for the R 3.5 analysis for umap must have python virtual env for this particular pipeline
# python must be installed for example: 
virtualenv_create("r-reticulate")
virtualenv_install("r-reticulate", "umap-learn")
use_virtualenv("r-reticulate")
library(umap)
library(reticulate)
library(umap)
# (Umap can now be called directly from most single cell analysis software). 
```

R packages used in this this analysis: R 4.0.5  
The same packages are required as above (separately installed for R
4.0.5) in addition, the following packages in the R4 envionmnent are
required:

``` r
# switch to R 4.0.5, i.e. using R switch https://rud.is/rswitch/guide/
pkgs4.0 = c(pkgs3.5, 'data.table', 'factoextra', 'GEOquery')
BiocManager::install("variancePartition")
devtools::install_github("caleblareau/BuenColors")
lapply(pkgs4.0, function(x){ 
  tryCatch(library(x), 
         error = function(e){
           install.packages(pkgs =  x, repos = 'http://cran.us.r-project.org')
           library(x)
         })
  })
```

### 3\) Download starting data and add to data directory <a name="data"></a>

After downloading the repository located at
<https://github.com/niaid/dsb_manuscript>, add the data folder to the
repository at the top level (where the file .rproj is located). Prior to
running any analysis confirm the data are in the correct repository with
the scripts below.

``` r
# confirm  *data/V2_Data/* and *data/background_data/*
suppressMessages(library(here))
# all STARTING data 
data_ = c(
  "CITEseq_raw_PMID32094927_seurat2.4.rds", 
  "unstained_control_singlets.rds"
)
stopifnot(data_ %in% list.files(here("data/V2_Data/")))

# background drop data 
data_ = c(
"adt_neg_dmx_list.rds",
"adt_neg_dmx.rds",
"adt_neg_full_list.rds",
"adt_neg_full_qc.rds"
)
stopifnot(data_ %in% list.files(here("data/V2_Data/background_data/")))

### public data directories

# mission bio 
data_ = c("AML-4-cell-line-multiomics-adt-counts.tsv")
stopifnot(data_ %in% list.files(here("data/mission_bio_data/")))
# 10X data 
data_ = c("10x_pbmc5k_V3.rds", "10x_pbmc5k_NextGem.rds",
          "10x_pbmc10k_V3.rds", "10x_pbmc_5prime_5k.rds")
stopifnot(data_ %in% list.files(here("data/10x_rds/")))
tenx_filtered_matrices = c('5prime_5k','v3_10k','v3_5k','nextgem_5k')
# filtereed matrices 
tenx_filtered_dir = list.dirs(path = here('data/10x_rds/'), full.names = FALSE)
stopifnot(tenx_filtered_matrice %in% tenx_filtered_dir)
```

## Analysis

These can be run line by line in an active R session or by sourcing the
script. As described above, file paths do not need to be
changed.

### Run dsb normalization on PBMC data from 20 individuals. <a name="dsb_1"></a>

V2/dsb\_normalize\_cluster\_pipeline/

``` r
# R 3.5.1 
source(here("V2/dsb_normalize_cluster_pipeline/1_dsb_normalize.r"))
```

### Run CLR normalization (across cells) on PBMC data from 20 individuals. <a name="clr_cells"></a>

Here switch to R 4.0.5, CLR transform across cells

``` r
# R 4.0.5 
source(here("V2/dsb_normalize_cluster_pipeline/1.1_calc_clr_across_cells_s4.r"))
```

**The rest of this section uses R 3.5.1**

### UMAP based on dsb normalized values <a name="umap"></a>

reads data
from:  
V2/dsb\_normalize\_cluster\_pipeline/generated\_data/h1\_d0\_singlets\_ann\_Seurat2.4\_dsbnorm.rds

``` r
# switch back to R 3.5
source(here("V2/dsb_normalize_cluster_pipeline/2_run_umap.r"))
```

### dsb vs CLR nk cluster comparison, manual gating, stained vs unstained distribution normalization comparison. <a name="pbmc_analysis"></a>

These scripts compare different normalizations with dsb and produce
various comparison analysis and visualization between methods.

``` r
source(here("V2/dsb_normalize_cluster_pipeline/3_figure_generation.r"))
source(here("V2/dsb_normalize_cluster_pipeline/4_manual_gate_plots.r"))
source(here("V2/dsb_normalize_cluster_pipeline/5_stained_unstained_dsb_distributions.r"))
# an alternate scheme for normalizing stained and unstained cells with highly concordant results (for reference)
source(here("V2/dsb_normalize_cluster_pipeline/5a_alternate_noise_unstained_plots.r"))
```

### dsb technical component robusness assessments, background noise comparison and robustness check <a name="robustness"></a>

This section illustrates the dsb process in steps and models each step
underlying the method. 6 and 6a analyze the correlation structure of
variables comprising the technical component in step II of the dsb
method. 7 models the correlation between multiple measurements of
experimentally and model derived background noise with ‘ground truth’
background from unstained control cells spiked into the cell pool prior
to droplet generation. In 8, the per-cell two component Gaussian mixture
model is compared to k = 1-6 component models and the resulting single
cell fits are analyzed.

``` r
source(here("V2/dsb_process_plots/6_mean_isotype_v_mean_control.R"))
source(here("V2/dsb_process_plots/6a_isotype_figure_generation.r"))
source(here("V2/dsb_process_plots/7_neg_control_plots.R"))
source(here("V2/dsb_process_plots/8_mixture_fits.r"))
```

### Multi vs single batch normalization, µ1 background resampling robustness check <a name="multibatch"></a>

empty\_drop\_threshold\_batch is an assessment of single vs multi batch
normalization and sensitivity of each normalization scheme to defining
background with hashing or library size distribution.
mu1\_noise\_correlations is analysis of the robustness of µ1 background
assessment and correlation with µ2 and isotype control means using 100
random sampling of µ1 proteins from each cell.

``` r
source(here("V2/parameter_sensitivity/empty_drop_threshold_batch.r"))
source(here("V2/parameter_sensitivity/mu1_noise_correlations.r"))
```

### External 10X genomics data analysis: “NextGem”, “V3”, and “5 Prime” assays <a name="tenx"></a>

These scripts are identical for each data set with tuned parameters at
the beginning of each script. Read Cell Ranger raw output, select
negative drops, run DSB normalization, run each normalization modeling
step separately, cluster cells and plot distributions across clusters
and on biaxial gates. Test underlying modeling assumptions for each
external dataset.

``` r
# 10K v3 data 
source(here("V2/10x_analysis/10x_pbmc_10k_V3.r"))
source(here("V2/10x_analysis/10x_pbmc_10k_V3_figure_generation.r"))
# 5k V3 data 
source(here("V2/10x_analysis/10x_pbmc_5k_V3.r"))
source(here("V2/10x_analysis/10x_pbmc_5k_V3_figure_generation.r"))
# Next Gem data 
source(here("V2/10x_analysis/10x_pbmc_NextGem.r"))
source(here("V2/10x_analysis/10x_pbmc_NextGem_figure_generation.r"))
# 5 prime data 
source(here("V2/10x_analysis/10x_pbmc_5prime_5k.r"))
source(here("V2/10x_analysis/10x_pbmc_5prime_5k_figure_generation.r"))
```

### dsb normalize protein data from Mission Bio tapestri platform <a name="missionbio"></a>

The data downloaded from MissionBio are reformatted for dsb and
normalized using dsb step I - ambient correction using empty droplets.

``` r
# tapestri example data dsb normalization 
source(here("V2/missionbio_tapestri/tapestri_exampledata_analysis.r"))
```

### dsb normalization of TEA-seq data and dsb-based WNN multimodal clustering <a name="teaseq"></a>

##### this analysis uses R version 4.0.5 with Seurat version 4.0.1

TEA-seq data were downloaded from
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5123951>. The
files:  
X066-MP0C1W3\_leukopak\_perm-cells\_tea\_200M\_cellranger-arc\_per\_barcode\_metrics.csv  
X066-MP0C1W3\_leukopak\_perm-cells\_tea\_200M\_adt\_counts.csv  
GSM5123951\_X066-MP0C1W3\_leukopak\_perm-cells\_tea\_200M\_cellranger-arc\_filtered\_feature\_bc\_matrix.h5  
downloaded from GEO into directory data/revision\_data/tea\_seq/  
Barcodes are filtered to cells and empty drops based on the provded
metadata is\_cell parameters.

``` r
# run custom mapping script provided by Lucas Graybuck
source(here('V2/teaseq/1_teaseq_preprocess_barcodes.r'))

# run customized joint protein RNA clustering based on Seurat 4 WNN with CLR and dsb for protein norm.   
source(here("V2/teaseq/SR4_teaseq_pipeline_V2.r"))

# figure generation
source(here("V2/teaseq/teaseq_figgen_V2.R"))

# dsb modeling assumptions 
source(here('V2/teaseq/teaseq_dsb_model.r'))
```

### dsb normalization of ASAP-seq data and dsb-based WNN multimodal clustering <a name="asapseq"></a>

##### this analysis uses R version 4.0.5 and Seurat version 4.0.1

ASAP-seq data was downloaded from
<https://github.com/caleblareau/asap_reproducibility/tree/master/bonemarow_asapseq/data>  
See: <https://www.nature.com/articles/s41587-021-00927-2>

``` r
# run dsb on ASAP-seq data 
source(here('V2/asapseq/asapseq_bm.R'))
# Figure generation 
source(here('V2/asapseq/asapseq_bm_figure_generation.r'))
```

### dsb vs CLR (acorss cells) Normalization comparison: Differential expression, Variance Paratition, Gap Statistic <a name="compare"></a>

Comparison of dsb with the updated implementation of CLR normalization
(across cells).

``` r
# Gap statistic for cluster quality 
source(here('V2/si/gap.r'))

# Variance partition analysis 
source(here('V2/si/vpart.r'))

# differential expression of cluster protein markers 
source(here('V2/si/de.r'))
```

### dsb vs CLR normalized values as input to WNN multimodal clustering: PBMC data from 20 donors <a name="wnn"></a>

``` r
# wnn analysis of CITE-seq data normalized wit h CLR vs dsb 
source(here('V2/joint_clustering/SR4_seurat_wnn_dsb.R'))
source(here('V2/joint_clustering/SR4_seurat_wnn_dsb_figgen.r'))
```

### dataset summary statistic table <a name="summarytable"></a>

``` r
# create summary statistics table for analyzed datasets. 
source(here('V2/table/make_table.r'))
```

### public dataset sources

Data are available in analysis-ready format at the figshare link above
for convenience. For 10X data, downloaded only the the *raw, not the
filtered data* in the links below downloaded **feature / cell matrix
(raw)** output. dsb uses the non cell containing empty droplets that are
in the raw output to estimate background - see package documentation.  
Next GEM:
<https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3_nextgem>  
V3 (5K):
<https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3>  
V3 (10K):
<https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3>  
5’
<https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_5gex_protein>  
Mission Bio Tapestri:
<https://portal.missionbio.com/datasets/4-cell-lines-AML-multiomics>  
download the following text file
‘AML-4-cell-line-multiomics-adt-counts.tsv’

### initialization script for public datasets

This code is here for reference. Conversion of public 10X genomics
datasets to .rds files.

``` r
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
source(here("V2/functions/preprocessing_functions.R"))

# Read10X_MPM is equivalent to Read10X in Seurat version 3 (reads compressed cell ranger data)
raw1 = Read10X_MPM("data/10x_data/10x_pbmc10k_V3/raw_feature_bc_matrix/")
raw2 = Read10X_MPM("data/10x_data/10x_pbmc5k_V3/raw_feature_bc_matrix/")
raw3 = Read10X_MPM("data/10x_data/10x_pbmc_5prime_5k/raw_feature_bc_matrix/")
raw4 = Read10X_MPM("data/10x_data/10x_pbmc5k_NextGem/raw_feature_bc_matrix/")

# save 
dir.create(here("data/10x_rds/"))
saveRDS(raw1,file = here("data/10x_rds/10x_pbmc10k_V3.rds"))
saveRDS(raw2,file = here("data/10x_rds/10x_pbmc5k_V3.rds"))
saveRDS(raw3,file = here("data/10x_rds/10x_pbmc_5prime_5k.rds"))
saveRDS(raw4,file = here("data/10x_rds/10x_pbmc5k_NextGem.rds"))
```

### initialization script for healthy donor 53k cell CITE-seq data.

The code below was run outside of this workflow for convenience, it is
shown here for reference. This step reformats the dataset hosted at the
[figshare
repository](https://nih.figshare.com/collections/Data_and_software_code_repository_for_Broad_immune_activation_underlies_shared_set_point_signatures_for_vaccine_responsiveness_in_healthy_individuals_and_disease_activity_in_patients_with_lupus_Kotliarov_Y_Sparks_R_et_al_Nat_Med_DOI_https_d/4753772)
that is associated with the
[manuscript](https://doi.org/10.1038/s41591-020-0769-8). The script
below removes the protein normalized data slot which used an earlier
version of the DSB package for normalization. In addition, normalized
RNA data, metadata, clustering snn graphs and tSNE results are removed
to reduce object size and cell type annotations as reported in the paper
linked above are added to the object as metadata.

``` r
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

datapath = here("data/V2_Data/")

# clear unneeded data slots tsne cell embeddings, unneeded metadata 
h1 = readRDS(file = "data/V2_Data/H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds")
h1@dr$protpca = NULL
h1@dr$tsne_p10 = NULL
h1@dr$tsne_p50 = NULL
h1@dr$tsne_p90 = NULL
h1@dr$tsne_p130 = NULL
h1@dr$tsne_p170 = NULL

h1@assay$HTO = NULL
h1@assay$CITE@data = NULL
h1@assay$CITE@scale.data = NULL
h1@data = NULL
# remove unused cells RNA data from raw slot
h1 = h1 %>% SubsetData(subset.raw = TRUE)

# remove unused additional metadata 
mds = read_delim(file = here("init_data_ignore/md_strip.txt"),delim = "\t")
vrm = mds$vars_
for (i in 1:length(vrm)) {
  vr = vrm[i]
  h1@meta.data[[vr]] = NULL
}

# add cluster labels from Fig 6 of https://www.nature.com/articles/s41591-020-0769-8 for reference in dsb paper.
cmd = read_delim(file = paste0(datapath, "clustree_node_labels_withCellTypeLabels.txt"), delim = '\t')
cmd1 = cmd %>% filter(clustering == "p3_dist_1")
cmd3 = cmd %>% filter(clustering == "p3_dist_3")

mdn = h1@meta.data %>% 
  mutate(celltype_label_1 = plyr::mapvalues(x = h1@meta.data$p3_dist_1, from = cmd1$cluster, to =cmd1$`Cell Type label`)) %>% 
  mutate(celltype_label_3 = plyr::mapvalues(x = h1@meta.data$p3_dist_3, from = cmd3$cluster, to =cmd3$`Cell Type label`)) %>% 
  select(barcode_check, celltype_label_3, celltype_label_1) %>% 
  column_to_rownames("barcode_check")

h1 = h1 %>% AddMetaData(metadata = mdn)

# save this starting object for DSB paper analysis  
saveRDS(h1, file = paste0(datapath, "CITEseq_raw_PMID32094927_seurat2.4.rds"))
```

### Project options and Sessioninfo

R version 3.5.3 attached base packages: stats graphics grDevices utils
datasets methods base

other attached packages: mclust\_5.4.5 reticulate\_1.12 umap\_0.2.3.1
magrittr\_1.5 forcats\_0.4.0 stringr\_1.4.0 dplyr\_0.8.5 purrr\_0.3.3
readr\_1.3.1 tidyr\_1.0.2 tibble\_2.1.1 tidyverse\_1.2.1 Seurat\_2.3.4
Matrix\_1.2-15 cowplot\_0.9.4 ggplot2\_3.1.1 here\_0.1

loaded via a namespace (and not attached): readxl\_1.3.1 snow\_0.4-3
backports\_1.1.4 Hmisc\_4.2-0 plyr\_1.8.4 igraph\_1.2.4.1
lazyeval\_0.2.2 splines\_3.5.3 inline\_0.3.15 digest\_0.6.25
foreach\_1.4.4 htmltools\_0.3.6 lars\_1.2 rsconnect\_0.8.16
gdata\_2.18.0 checkmate\_1.9.3 cluster\_2.0.7-1 mixtools\_1.1.0
ROCR\_1.0-7 modelr\_0.1.4 matrixStats\_0.54.0 R.utils\_2.8.0
askpass\_1.1 prettyunits\_1.0.2 colorspace\_1.4-1 rvest\_0.3.4
haven\_2.1.0 xfun\_0.7 callr\_3.2.0 crayon\_1.3.4 jsonlite\_1.6
survival\_2.43-3 zoo\_1.8-6 iterators\_1.0.10 ape\_5.3 glue\_1.3.1
gtable\_0.3.0 pkgbuild\_1.0.3 kernlab\_0.9-27 rstan\_2.19.3
prabclus\_2.3-1 DEoptimR\_1.0-8 scales\_1.0.0 mvtnorm\_1.0-10
bibtex\_0.4.2 Rcpp\_1.0.1 metap\_1.1 dtw\_1.20-1 htmlTable\_1.13.1
foreign\_0.8-71 bit\_1.1-14 proxy\_0.4-23 SDMTools\_1.1-221.1
Formula\_1.2-3 stats4\_3.5.3 tsne\_0.1-3 StanHeaders\_2.21.0-1
htmlwidgets\_1.3 httr\_1.4.0 gplots\_3.0.1.1 RColorBrewer\_1.1-2
fpc\_2.2-1 acepack\_1.4.1 modeltools\_0.2-22 ica\_1.0-2 pkgconfig\_2.0.2
loo\_2.3.1 R.methodsS3\_1.7.1 flexmix\_2.3-15 nnet\_7.3-12
tidyselect\_0.2.5 rlang\_0.4.5 reshape2\_1.4.3 cellranger\_1.1.0
munsell\_0.5.0 tools\_3.5.3 cli\_1.1.0 generics\_0.0.2 broom\_0.5.2
ggridges\_0.5.1 evaluate\_0.14 yaml\_2.2.0 npsurv\_0.4-0 processx\_3.3.1
knitr\_1.23 bit64\_0.9-7 fitdistrplus\_1.0-14 robustbase\_0.93-5
caTools\_1.17.1.2 RANN\_2.6.1 packrat\_0.5.0 pbapply\_1.4-0
nlme\_3.1-137 R.oo\_1.22.0 xml2\_1.2.0 hdf5r\_1.2.0 compiler\_3.5.3
rstudioapi\_0.10 png\_0.1-7 lsei\_1.2-0 stringi\_1.4.3 ps\_1.3.0
lattice\_0.20-38 vctrs\_0.2.4 pillar\_1.4.1 lifecycle\_0.1.0
Rdpack\_0.11-0 lmtest\_0.9-37 data.table\_1.12.2 bitops\_1.0-6
irlba\_2.3.3 gbRd\_0.4-11 R6\_2.4.0 latticeExtra\_0.6-28
KernSmooth\_2.23-15 gridExtra\_2.3 codetools\_0.2-16 MASS\_7.3-51.1
gtools\_3.8.1 assertthat\_0.2.1 openssl\_1.4 rprojroot\_1.3-2
withr\_2.1.2 diptest\_0.75-7 parallel\_3.5.3 doSNOW\_1.0.16 hms\_0.4.2
grid\_3.5.3 rpart\_4.1-13 class\_7.3-15 rmarkdown\_1.13
segmented\_0.5-4.0 Rtsne\_0.15 lubridate\_1.7.4 base64enc\_0.1-3

### NIAID repository release notes

A review of this code has been conducted, no critical errors exist, and
to the best of the authors knowledge, there are no problematic file
paths, no local system configuration details, and no passwords or keys
included in this code.  
Primary author(s): Matt Mulè  
Organizational contact information: General: john.tsang AT nih.gov,
code: mulemp AT nih.gov \[permanent address mattmule AT gmail\]  
Date of release: Oct 7 2020  
Description: code to reproduce analysis of manuscript  
Usage instructions: Provided in this markdown

### code check

Checked repository for PII and strings containing file paths. Data used
in this analysis does not contain PII.

``` r
library(lintr)
fcn = suppressMessages(list.files(here("functions"), ".r", full.names = TRUE))
pipel = suppressMessages(list.files(here("V2/dsb_normalize_cluster_pipeline/"),pattern = ".r", full.names = TRUE))
process = suppressMessages(list.files(here("V2/dsb_process_plots/"),pattern = ".r", full.names = TRUE))
tenx = suppressMessages(list.files(here("V2/10x_analysis/"),".r", full.names = TRUE))
mb = suppressMessages(list.files(here("V2/missionbio_tapestri/"), ".r", full.names = TRUE))
param = suppressMessages(list.files(here("V2/parameter_sensitivity/"),pattern = ".r", full.names = TRUE))

# code check
scp = c(fcn, tenx, mb, pipel, process, param) %>% as.list()
lt = suppressMessages(lapply(scp, lintr::lint))
```
