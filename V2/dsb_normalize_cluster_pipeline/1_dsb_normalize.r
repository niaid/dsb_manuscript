suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# load DSB package 
library(dsb)

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/"); dir.create(figpath)
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/"); dir.create(datapath)


# read in final singlets make list of seurat objects indexed by batch. Subset out raw adt data from neg drop subset. 
h1 = readRDS(file = here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")) %>% SetAllIdent(id = "batch")
h1b1 = SubsetData(h1, ident.use = "1")
h1b2 = SubsetData(h1, ident.use = "2")
cells = list(h1b1, h1b2)
rm(h1b1, h1b2); gc()
pos_adt = lapply(cells, function(x){x@assay$CITE@raw.data} %>% as.matrix())

# load ADT data from negative droplet subset   
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# apply denoised scaled by background protein normalization per batch 
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")
dsb_norm = list()
for (i in 1:length(neg_adt)) {
  dsb_norm[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                       empty_drop_matrix = neg_adt[[i]], 
                                       denoise.counts = TRUE, 
                                       use.isotype.control = TRUE, 
                                       isotype.control.name.vec = isotypes)
}

# merge multi batch norm data and ad to Seurat object 
dsb = do.call(cbind, dsb_norm)
h1 = SetAssayData(h1, assay.type = "CITE", new.data = dsb, slot = "data")

# save normalized counts and Seurat object
saveRDS(dsb_norm, file = paste0(datapath,"dsb2_normlog10_list_bybatch.rds"))
saveRDS(h1, file = paste0(datapath,"h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))

# Run DSB normalization without denoising (in comparison figure 1h)
dsb_norm2 = list()
for (i in 1:length(neg_adt)) {
  dsb_norm2[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                        empty_drop_matrix = neg_adt[[i]], 
                                        denoise.counts = FALSE)
}
dsb2 = do.call(cbind, dsb_norm2)
saveRDS(dsb2, file = paste0(datapath,"NonDenoised_dsb_Mtx.rds"))


