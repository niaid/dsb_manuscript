# run umap on DSB normalized values 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# load python virtual environment for python umap learn via reticulate 
suppressMessages(library(reticulate))
use_virtualenv("r-reticulate")
library(umap)

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# DSB protein data matrix without isotype controls 
h1 = readRDS("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds")
proteins = rownames(h1@assay$CITE@data)
proteins = proteins[-c(67:70)]
adt = h1@assay$CITE@data[proteins, ]
adt = t(adt)

# Set parameters and run umap 
config = umap.defaults
config$n_neighbors = 35
config$min_dist = 0.6
ump = umap(adt,config = config)

# extract cell embeddings 
embeddings = ump$layout %>% as.data.frame()
colnames(embeddings) = c("UMAP1", "UMAP2")

# save combined data 
df = cbind(embeddings, h1@meta.data, as.data.frame(t(h1@assay$CITE@data)))
saveRDS(df, file = paste0(datapath, "h1_sng_metadata_umapdim_prots_dataframe.rds"))
