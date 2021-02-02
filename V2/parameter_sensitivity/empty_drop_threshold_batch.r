suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(dsb))

# file paths 
figpath = here("V2/parameter_sensitivity/figures/")
datapath = here("V2/parameter_sensitivity/generated_data/")
dir.create(figpath); dir.create(datapath)

######## Change to new versions
neg_adt = readRDS(file = here("data/V2_Data/background_data/adt_neg_full_list.rds"))
negadt_sub = readRDS(file = here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# reformat 
neg_adt[[1]] = neg_adt[[1]][ ,-1] %>% column_to_rownames("bc") %>% t 
neg_adt[[2]] = neg_adt[[2]][ ,-1] %>% column_to_rownames("bc") %>% t 


# list of empty drops at the two thresholds 
em = list(neg_adt[[1]], neg_adt[[2]], negadt_sub[[1]], negadt_sub[[2]])
names(em) = c("batch 1 threshold 1", "batch 2 threshold 1", "batch 1 threshold 2", "batch 2 threshold 2")

# plot empty drp library sizes for each batch 
colvec = pals::tableau20(4)

for (i in 1:length(em)) {
  pdf(file = paste0(figpath, names(em[i]), "hist.pdf"), width = 5, height = 5)
  hist( log10(colSums(em[[i]])), 
        main = paste0(ncol(em[[i]]), " empty drops ", names(em[i])),
        xlim = c(2.5,4), cex.main = 1.3,
        xlab = " droplet  log10 protein library size ",
        breaks = 100 , col = colvec[i])
  dev.off()
}


# multibatch norm with full drop background 
h1 = readRDS(file = here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")) %>% SetAllIdent(id = "batch")
meta = h1@meta.data

# define main lineage same as in fig 1 
celltypes = meta$celltype_label_1 %>% unique() %>% sort()
tcell = celltypes[c(2,3,4,5,10)]
myeloid = celltypes[c(6,8,9)]
bcell = celltypes[c(1)]
nk = celltypes[c(7)]

# add lineage metadata 
meta = meta %>%  
  mutate(lineage = 
          if_else(celltype_label_1 %in% tcell, "T_Cell",
          if_else(celltype_label_1 %in% myeloid, "Myeloid_lineage",
          if_else(celltype_label_1 %in% bcell, "B_Cell",false = "other")
          ))) 

# get raw data for each batch 
h1 = SetAllIdent(h1, id = "batch")  
h1b1 = SubsetData(h1, ident.use = "1")
h1b2 = SubsetData(h1, ident.use = "2")
cells = list(h1b1, h1b2)
rm(h1b1, h1b2,h1 ); gc()
pos_adt = lapply(cells, function(x){x@assay$CITE@raw.data} %>% as.matrix())


# apply DSB using full background for each batch separately.  
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",
             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")
dsb_norm = list()
for (i in 1:length(neg_adt)) {
  dsb_norm[[i]] = 
    DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                        empty_drop_matrix = neg_adt[[i]], 
                        denoise.counts = TRUE, 
                        use.isotype.control = TRUE, 
                        isotype.control.name.vec = isotypes)
}
# merge multi batch norm data into "dsb"
dsb = do.call(cbind, dsb_norm)


## Run single batch norm using full definition of background 
full_neg_merged = do.call(cbind, neg_adt)
pos_adt_merged = do.call(cbind, pos_adt)
dsb_merged_full = DSBNormalizeProtein(cell_protein_matrix = pos_adt_merged, 
                                      empty_drop_matrix = full_neg_merged, 
                                      denoise.counts = TRUE, 
                                      use.isotype.control = TRUE, 
                                      isotype.control.name.vec = isotypes)

#multibatch norm using subset of background is stored in the normalized data object 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
cite_dsb = h1@assay$CITE@data


## single batch norm subset background 
sub_neg_merged = do.call(cbind, negadt_sub)
dsb_merged_sub = DSBNormalizeProtein(cell_protein_matrix = pos_adt_merged, 
                                     empty_drop_matrix = sub_neg_merged, 
                                     denoise.counts = TRUE, 
                                     use.isotype.control = TRUE, 
                                     isotype.control.name.vec = isotypes)



####### 
##### Part II visualization 
## histograms of protein distributions by batch 
proteins = rownames(dsb)
index1 = proteins[1]; index2 = proteins[length(proteins)]

# tidy
full_multibatch = cbind(meta, as.data.frame(t(dsb))) %>% 
  select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

full_mergedbatch = cbind(meta, as.data.frame(t(dsb_merged_full))) %>% 
  select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

### as defined above 
# h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
# cite_dsb = h1@assay$CITE@data
sub_multibatch = cbind(meta, as.data.frame(t(cite_dsb))) %>% 
  select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

sub_mergedbatch = cbind(meta, as.data.frame(t(dsb_merged_sub))) %>% 
  select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

# list of data tidy by normalization procedure 
ml = list(full_multibatch, full_mergedbatch,  sub_multibatch, sub_mergedbatch)
names(ml) = c("full-multi-batch", "full-merged-batch", "subset-multi-batch", "subset-merged-batch" )

# Plot celltype annotation heatmap 
mdf = cbind(meta, as.data.frame(t(cite_dsb))) %>% 
  select(proteins, celltype_label_3) %>% 
  gather(protein, dsb, index1:index2) %>% 
  group_by(celltype_label_3, protein)%>% 
  summarize(dsbmean = mean(dsb)) %>% 
  spread(protein, dsbmean) %>% 
  column_to_rownames("celltype_label_3")
xx=pheatmap::pheatmap(mdf, color = viridis::inferno(15),
                      fontsize_col = 5,fontsize_row = 7, 
                      width = 10, height = 4.5, border_color = NA,
                      filename = paste0(figpath, "dsb_heatmap.pdf"))

pt= xx$tree_col$labels[xx$tree_col$order]

bcp = pt[c(1, 3,4:10,82, 27)]
tcp = pt[c(18,21,22,24,36,64,73:81)]
mcp = pt[c(12:16, 62, 65, 67:71)]

ml = lapply(ml, function(x){ 
  x = x %>% 
    mutate(prot_class = 
          if_else(protein %in% bcp, "bcell_protein", 
          if_else(protein %in% tcp, "tcell_protein", 
          if_else(protein %in% mcp, "myeloid_protein", false = "other"))))
  })

ridge_layers = list(
  geom_vline(xintercept = 0, color = "red"),  
  theme(axis.text.y = element_text(size = 15)), 
  xlab(""),
  ylab(""), 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  
)

for (i in 1:length(ml)) {
  plot_data = ml[[i]] %>% filter(DSB_norm > -5) %>% filter(DSB_norm < 35) 
  
  # T cell subset 
  p1 = ggplot(plot_data %>% filter(lineage == "T_Cell" & prot_class == "tcell_protein"), 
              aes(x = DSB_norm, y = protein, fill = batch)) +
    ggtitle("T cells") + 
    ggridges::geom_density_ridges2(alpha = 0.7, show.legend = FALSE) + 
    ridge_layers
  
  # B cell subset 
  p2 = ggplot(plot_data %>% filter(lineage == "B_Cell" & prot_class == "bcell_protein"), 
              aes(x = DSB_norm, y = protein, fill = batch)) +
    ggtitle("B cells") + 
    ggridges::geom_density_ridges2(alpha = 0.7, show.legend = FALSE) + 
    ridge_layers

  # myeloid subset 
  p3 = ggplot(plot_data %>% filter(lineage == "Myeloid_lineage" & prot_class == "myeloid_protein"), 
              aes(x = DSB_norm, y = protein, fill = batch)) +
    ggtitle("Myeloid Cells") + 
    ggridges::geom_density_ridges2(alpha = 0.7, show.legend = FALSE) + 
    ridge_layers
  
  #merge plots and save 
  p4 = cowplot::plot_grid(plotlist = list(p1,p2,p3),  ncol = 1, rel_heights = c(0.1,0.1, 0.1)) 
  ggsave(p4,filename = paste0(figpath, names(ml[i]), "hist_batch.pdf"), width = 3.5, height = 12)
    
}


