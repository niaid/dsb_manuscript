suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(dsb))

# file paths 
figpath = here("V2/parameter_sensitivity/figures/"); dir.create(figpath)
datapath = here("V2/parameter_sensitivity/generated_data/"); dir.create(datapath)

# load background data 
neg_adt = readRDS(file = here("data/V2_Data/background_data/adt_neg_full_list.rds"))
negadt_sub = readRDS(file = here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# adt_neg_full_list.rds is the full background
# negadt_sub = the hashing negatives 
# each is a list indexed by batch 
print(names(neg_adt))
# "batch 1 Neg drops FULL QCd" "batch 2 Neg drops FULL QCd"
print(names(negadt_sub))
# "batch 1 hash DMX Neg QCd" "batch 2 hash DMX Neg QCd"

# reformat 
neg_adt[[1]] = neg_adt[[1]][ ,-1] %>% column_to_rownames("bc") %>% t() 
neg_adt[[2]] = neg_adt[[2]][ ,-1] %>% column_to_rownames("bc") %>% t() 


# list of empty drops at the two thresholds 
em = list(neg_adt[[1]], neg_adt[[2]], negadt_sub[[1]], negadt_sub[[2]])
names(em) = c("batch 1 threshold 1", "batch 2 threshold 1", "batch 1 threshold 2", "batch 2 threshold 2")


# plot empty drop library sizes for each definition of background split by batch 
colvec = pals::tableau20(4)
colvec = colvec[c(1,3,2,4)]
for (i in 1:length(em)) {
  pdf(file = paste0(figpath, names(em[i]), "hist.pdf"), width = 5, height = 5)
  hist( log10(colSums(em[[i]])), 
        main = paste0(ncol(em[[i]]), " empty drops ", names(em[i])),
        xlim = c(2.5,4), cex.main = 1.3,
        xlab = " droplet  log10 protein library size ",
        breaks = 40 , col = colvec[i])
  dev.off()
}


## Normalization with dsb using split vs merged batch normalization
# both with hashing and library size negatives 


#load raw data 
h1 = readRDS(file = here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")) %>% SetAllIdent(id = "batch")
meta = h1@meta.data

# define main lineage same as in fig 1 (for assessment of protein distribution overlap within lineages)
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

# get raw data ADT counts for each batch 
h1 = SetAllIdent(h1, id = "batch")  
h1b1 = SubsetData(h1, ident.use = "1")
h1b2 = SubsetData(h1, ident.use = "2")
cells = list(h1b1, h1b2)
rm(h1b1, h1b2,h1 ); gc()
pos_adt = lapply(cells, function(x){x@assay$CITE@raw.data} %>% as.matrix())

## 1. Split batch normalize each batch separately with library size (full) negatives from that batch 
# apply dsb using full background for each batch separately.  
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
# merge multi batch norm data into single matrix `dsb`
dsb = do.call(cbind, dsb_norm)



## 2. Merge first then run single normalization with library size (full) negatives from both batches 
## Run single batch norm 
full_neg_merged = do.call(cbind, neg_adt)
pos_adt_merged = do.call(cbind, pos_adt)
dsb_merged_full = DSBNormalizeProtein(cell_protein_matrix = pos_adt_merged, 
                                      empty_drop_matrix = full_neg_merged, 
                                      denoise.counts = TRUE, 
                                      use.isotype.control = TRUE, 
                                      isotype.control.name.vec = isotypes)



## 3. Split batch normalize each batch separately with Hashing negatives from that batch 
#multibatch norm using hashing negatives is stored in the normalized data slot 
# from dsb_normalize_cluster_pipeline/1_dsb_normalize.r
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
cite_dsb = h1@assay$CITE@data



## 4. Merge first then run single normalization with hashing negatives from both batches 
## single batch norm subset background 
sub_neg_merged = do.call(cbind, negadt_sub)
system.time({
dsb_merged_sub = DSBNormalizeProtein(cell_protein_matrix = pos_adt_merged, 
                                     empty_drop_matrix = sub_neg_merged, 
                                     denoise.counts = TRUE, 
                                     use.isotype.control = TRUE, 
                                     isotype.control.name.vec = isotypes)

})

# in order aobve  dsb, dsb_merged_full cite_dsb,  dsb_merged_sub

####### 
##### Part II visualization 
## histograms of protein distributions by batch 
proteins = rownames(dsb)
index1 = proteins[1]; index2 = proteins[length(proteins)]

# analysis of distributions from different versions of normalization 
# full_ prefix refers to 'full' background (library size negatives )
full_multibatch = cbind(meta, as.data.frame(t(dsb))) %>% 
  dplyr::select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

full_mergedbatch = cbind(meta, as.data.frame(t(dsb_merged_full))) %>% 
  dplyr::select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

sub_multibatch = cbind(meta, as.data.frame(t(cite_dsb))) %>% 
  dplyr::select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

sub_mergedbatch = cbind(meta, as.data.frame(t(dsb_merged_sub))) %>% 
  dplyr::select(lineage, proteins, batch, barcode_check) %>% 
  gather(protein, DSB_norm, index1:index2)

# list of data tidy by normalization procedure 
ml = list(full_multibatch, full_mergedbatch,  sub_multibatch, sub_mergedbatch)
names(ml) = c("full-multi-batch", "full-merged-batch", "subset-multi-batch", "subset-merged-batch" )
saveRDS(object = ml,file = paste0(datapath, 'ml.rds'))

# using cite_dsb from main normalization, plot heatmap for Supplemental Figure 6c
# Plot celltype annotation heatmap 
mdf = cbind(meta, as.data.frame(t(cite_dsb))) %>% 
  dplyr::select(proteins, celltype_label_3) %>% 
  gather(protein, dsb, index1:index2) %>% 
  group_by(celltype_label_3, protein)%>% 
  summarize(dsbmean = mean(dsb)) %>% 
  spread(protein, dsbmean) %>% 
  column_to_rownames("celltype_label_3")
colnames(mdf) = str_remove(string = colnames(mdf), pattern = "_PROT")
xx=pheatmap::pheatmap(mdf, color = viridis::inferno(15), 
                      fontsize_col = 6,fontsize_row = 7, 
                      treeheight_row = 12, treeheight_col = 18,
                      width = 9, height = 3.8, border_color = NA,
                      filename = paste0(figpath, "dsb_heatmap.pdf"))

# define lineage defining proteins  
pt = xx$tree_col$labels[xx$tree_col$order]
bcp = pt[c(1, 3,4:10,82, 27)]; bcp = paste(bcp,"_PROT",sep = "")
tcp = pt[c(18,21,22,24,36,64,73:81)]; tcp = paste(tcp,"_PROT",sep = "")
mcp = pt[c(12:16, 62, 65, 67:71)]; mcp = paste(mcp,"_PROT",sep = "")

# add protein category to each normalization 
ml = lapply(ml, function(x){ 
  x = x %>% 
    mutate(prot_class = 
          if_else(protein %in% bcp, "bcell_protein", 
          if_else(protein %in% tcp, "tcell_protein", 
          if_else(protein %in% mcp, "myeloid_protein", false = "other"))))
  })

# set plot theme 
ridge_layers = list(
  geom_vline(xintercept = 0, color = "red"),  
  theme(axis.text.y = element_text(size = 13, color = 'black')), 
  xlab(""),
  ylab(""), 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
)

for (i in 1:length(ml)) {
  
  plot_data = ml[[i]] %>% filter(DSB_norm > -5) %>% filter(DSB_norm < 35) 
  
  # T cell subset 
  p1 = ggplot(plot_data %>% filter(lineage == "T_Cell" & prot_class == "tcell_protein"), 
              aes(x = DSB_norm, y = protein, color = batch , fill = batch)) +
    ggtitle("T cells") + 
    ggridges::geom_density_ridges2(alpha = 0.4, size = 1, show.legend = FALSE) + 
    theme_bw() + 
    ridge_layers +
    scale_fill_manual(values = c(colvec[c(1,2)])) + 
    scale_color_manual(values = c(colvec[c(1,2)]))
  
  
  # B cell subset 
  p2 = ggplot(plot_data %>% filter(lineage == "B_Cell" & prot_class == "bcell_protein"), 
              aes(x = DSB_norm, y = protein, color = batch, fill = batch)) +
    ggtitle("B cells") + 
    ggridges::geom_density_ridges2(alpha = 0.4, size= 1, show.legend = FALSE) + 
    theme_bw() + 
    ridge_layers +
    scale_fill_manual(values = c(colvec[c(1,2)])) + 
    scale_color_manual(values = c(colvec[c(1,2)]))

  # myeloid subset 
  p3 = ggplot(plot_data %>% filter(lineage == "Myeloid_lineage" & prot_class == "myeloid_protein"), 
              aes(x = DSB_norm, y = protein, color = batch, fill = batch)) +
    ggtitle("Myeloid Cells") + 
    ggridges::geom_density_ridges2(alpha = 0.4, size= 1,show.legend = FALSE) + 
    theme_bw() + 
    ridge_layers +
    scale_fill_manual(values = c(colvec[c(1,2)])) + 
    scale_color_manual(values = c(colvec[c(1,2)]))
  
  #merge plots and save 
  p4 = cowplot::plot_grid(plotlist = list(p1,p2,p3),  ncol = 1, rel_heights = c(0.1,0.1, 0.1)) 
  ggsave(p4,filename = paste0(figpath, names(ml[i]), "hist_batch.pdf"), width = 3.5, height = 12)
}

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
#   [1] dsb_0.2.0       here_0.1        magrittr_2.0.1  forcats_0.4.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3     readr_1.3.1    
# [9] tidyr_1.0.2     tibble_2.1.1    tidyverse_1.2.1 Seurat_2.3.4    Matrix_1.2-15   cowplot_0.9.4   ggplot2_3.1.1  
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1            snow_0.4-3              backports_1.1.4         Hmisc_4.2-0             plyr_1.8.4             
# [6] igraph_1.2.4.1          lazyeval_0.2.2          splines_3.5.3           crosstalk_1.0.0         digest_0.6.25          
# [11] foreach_1.4.4           htmltools_0.3.6         viridis_0.5.1           lars_1.2                gdata_2.18.0           
# [16] checkmate_1.9.3         cluster_2.0.7-1         mixtools_1.1.0          ROCR_1.0-7              limma_3.38.3           
# [21] modelr_0.1.4            R.utils_2.8.0           colorspace_1.4-1        rvest_0.3.4             haven_2.1.0            
# [26] xfun_0.7                crayon_1.3.4            jsonlite_1.6            survival_2.43-3         zoo_1.8-6              
# [31] iterators_1.0.10        ape_5.3                 glue_1.3.1              pals_1.5                gtable_0.3.0           
# [36] webshot_0.5.1           kernlab_0.9-27          prabclus_2.3-1          DEoptimR_1.0-8          maps_3.3.0             
# [41] scales_1.0.0            pheatmap_1.0.12         mvtnorm_1.0-10          bibtex_0.4.2            miniUI_0.1.1.1         
# [46] Rcpp_1.0.1              metap_1.1               dtw_1.20-1              viridisLite_0.3.0       xtable_1.8-4           
# [51] htmlTable_1.13.1        reticulate_1.12         foreign_0.8-71          bit_1.1-14              mapproj_1.2.6          
# [56] proxy_0.4-23            mclust_5.4.5            SDMTools_1.1-221.1      Formula_1.2-3           stats4_3.5.3           
# [61] tsne_0.1-3              htmlwidgets_1.3         httr_1.4.0              gplots_3.0.1.1          RColorBrewer_1.1-2     
# [66] fpc_2.2-1               ellipsis_0.3.0          acepack_1.4.1           modeltools_0.2-22       ica_1.0-2              
# [71] pkgconfig_2.0.2         R.methodsS3_1.7.1       flexmix_2.3-15          nnet_7.3-12             labeling_0.3           
# [76] manipulateWidget_0.10.0 tidyselect_0.2.5        rlang_0.4.5             reshape2_1.4.3          later_0.8.0            
# [81] munsell_0.5.0           cellranger_1.1.0        tools_3.5.3             cli_1.1.0               generics_0.0.2         
# [86] broom_0.5.2             ggridges_0.5.1          npsurv_0.4-0            knitr_1.23              bit64_0.9-7            
# [91] fitdistrplus_1.0-14     robustbase_0.93-5       rgl_0.100.30            caTools_1.17.1.2        RANN_2.6.1             
# [96] packrat_0.5.0           pbapply_1.4-0           nlme_3.1-137            mime_0.6                R.oo_1.22.0            
# [101] xml2_1.2.0              hdf5r_1.2.0             compiler_3.5.3          rstudioapi_0.10         png_0.1-7              
# [106] lsei_1.2-0              stringi_1.4.3           lattice_0.20-38         vctrs_0.2.4             pillar_1.4.1           
# [111] lifecycle_0.1.0         Rdpack_0.11-0           lmtest_0.9-37           data.table_1.12.2       bitops_1.0-6           
# [116] irlba_2.3.3             gbRd_0.4-11             httpuv_1.5.1            R6_2.4.0                latticeExtra_0.6-28    
# [121] promises_1.0.1          KernSmooth_2.23-15      gridExtra_2.3           codetools_0.2-16        dichromat_2.0-0        
# [126] MASS_7.3-51.1           gtools_3.8.1            assertthat_0.2.1        rprojroot_1.3-2         withr_2.1.2            
# [131] diptest_0.75-7          parallel_3.5.3          doSNOW_1.0.16           hms_0.4.2               grid_3.5.3             
# [136] rpart_4.1-13            class_7.3-15            segmented_0.5-4.0       Rtsne_0.15              shiny_1.3.2            
# [141] lubridate_1.7.4         base64enc_0.1-3   


