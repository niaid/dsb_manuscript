# S4 R4 
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(mclust))
suppressMessages(library(Seurat))
library(magrittr)
library(variancePartition)

# proj opts 
figpath = here("V2/si/figures/"); dir.create(figpath)
datapath = here("V2/si/generated_data/"); dir.create(datapath)

s = readRDS(file = here('data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds'))
raw = s@assay$CITE@raw.data
rawprot = log10(Matrix::colSums(raw))
dl = readRDS(file = here('V2/si/figures/dat_list.rds'))
d1 = dl$dsb
d2 = dl$clr
met = dl$meta
met$prot_library_size = rawprot

# define main lineage same as in fig 1 
celltypes = met$celltype_label_1 %>% unique() %>% sort()
tcell = celltypes[c(2,3,4,5,10)]
myeloid = celltypes[c(6,8,9)]
bcell = celltypes[c(1)]
nk = celltypes[c(7)]

# add lineage metadata 
met = met %>%  
  mutate(lineage = 
  if_else(celltype_label_1 %in% tcell, "T_Cell",
  if_else(celltype_label_1 %in% myeloid, "Myeloid_lineage",
  if_else(celltype_label_1 %in% bcell, "B_Cell",false = "other")))) 

# split into main lineages 
ml = split(met, f = met$lineage)
cells = lapply(ml, function(x) x$barcode_check )

# list of mxt by cell type 
dl = dl2 = list()
for (i in 1:length(cells)) {
  dl[[i]] = d1[ ,cells[[i]]]
  dl2[[i]] = d2[ ,cells[[i]]]
}

# model formula 
f1 = ~ prot_library_size + (1|sampleid) + (1|batch) + (1|tenx_lane) + (1|celltype_label_3)

# variance partition 
dsb_ct = clr_ct = list()
for (i in 1:length(ml)) {
  dsb_ct[[i]] = fitExtractVarPartModel(exprObj = dl[[i]], data = ml[[i]], formula = f1, REML = FALSE, useWeights = FALSE) 
  clr_ct[[i]] = fitExtractVarPartModel(exprObj = dl2[[i]], data = ml[[i]], formula = f1, REML = FALSE, useWeights = FALSE)
}

for (i in 1:length(cells)) {
  dsb_ct[[i]] = dsb_ct[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column('prot') %>%
    mutate(lineage = names(cells)[i]) %>% 
    mutate(norm = 'dsb') %>% 
    gather(variable, value, batch:Residuals)
  clr_ct[[i]] = clr_ct[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column('prot') %>%
    mutate(lineage = names(cells)[i]) %>% 
    mutate(norm = 'CLR') %>% 
    gather(variable, value, batch:Residuals)
}

# combine data 
dsb = do.call(rbind, dsb_ct)
clr = do.call(rbind, clr_ct)
d = rbind(dsb, clr)
saveRDS(d, paste0(datapath,'d_celltype_varpartmodel_tidy_clr_dsb.rds'))
d = readRDS(file = here('V2/si/generated_data/d_celltype_varpartmodel_tidy_clr_dsb.rds'))

tcp  = c('TCRgd_PROT', 'KLRG1_PROT', 'CD8_PROT', 'CD7_PROT', 'CD5_PROT',
         'CD45RO_PROT', 'CD45RA_PROT', 'CD4_PROT', 'CD3_PROT', 'CD28_PROT', 
         'CD279_PROT', 'CD278 _PROT', 'CD27_PROT', 'CD244_PROT', 'CD127_PROT', 
         'CD38_PROT', "CD39_PROT", "CD69_PROT")

bcp = c('IgM_PROT', 'IgD_PROT', 'CD40_PROT', 'CD32_PROT', 
        'CD21_PROT', 'CD20_PROT', 'CD1c_PROT', 'CD19_PROT', 'CD185_PROT', 
        'CD38_PROT', "CD39_PROT", "CD69_PROT")

mcp = c('CD86_PROT', 'CD64_PROT', 'CD33_PROT', 'CD303_PROT', 'CD1d_PROT',
        'CD163_PROT', 'CD141_PROT', 'CD14_PROT', 'CD13_PROT', 'CD123_PROT',
        'CD11c_PROT', 'CD11b_PROT', 
        'CD38_PROT', "CD39_PROT", "CD69_PROT")
prot_combined = c(tcp, bcp, mcp) %>% unique 

d2 = d 
d2$prot = str_sub(d2$prot, 1 , -6)
mcp = str_sub(mcp, 1 , -6)
tcp = str_sub(tcp, 1 , -6)
bcp = str_sub(bcp, 1 , -6)

###### 
# visualizations 
#####

# myeloid 
porder = d2 %>% filter(norm == "dsb" & lineage == "Myeloid_lineage" & variable == 'prot_library_size') %>%
  arrange(value) %$% prot
d2$prot = factor(d2$prot, levels = porder)
p = 
  ggplot(d2 %>% filter(prot %in% mcp & variable == 'prot_library_size' & lineage == "Myeloid_lineage"),
         aes(x = prot, y = value*100, fill = norm )) +
  theme_bw() + 
  ggtitle(label = 'within Myeloid model') + 
  theme(title = element_text(size = 8)) + 
  ylab("% of variance explained\n by sequencing depth") + 
  geom_col(position = 'dodge') +
  theme(axis.text.y = element_text(color = 'black')) + 
  theme(axis.text.x = element_text(color = 'black')) + 
  theme(legend.position = 'none') +
  theme(axis.title.y = element_blank())  +
  scale_fill_manual(values = c('grey60', 'deepskyblue3')) + 
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, 'myeloid_var_explained.pdf'), width = 2, height = 3.3)


# B cell 
porder = d2 %>% filter(norm == "dsb" & lineage == "B_Cell" & variable == 'prot_library_size') %>%
  arrange(value) %$% prot
d2$prot = factor(d2$prot, levels = porder)
p = 
  ggplot(d2 %>% filter(prot %in% bcp & variable == 'prot_library_size' & lineage == "B_Cell"),
         aes(x = prot, y = value*100, fill = norm )) +
  theme_bw() + 
  ggtitle(label = 'within B cell model') + 
  theme(title = element_text(size = 8)) + 
  ylab("% of variance explained\n by sequencing depth") + 
  geom_col(position = 'dodge') +
  theme(axis.text.y = element_text(color = 'black')) + 
  theme(axis.text.x = element_text(color = 'black')) + 
  theme(legend.position = 'none') +
  theme(axis.title.y = element_blank())  +
  scale_fill_manual(values = c('grey60', 'deepskyblue3')) + 
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, 'Bcell_var_explained.pdf'), width = 2, height = 3.3)


# T cell 
porder = d2 %>% filter(norm == "dsb" & lineage == "T_Cell" & variable == 'prot_library_size') %>%
  arrange(value) %$% prot
d2$prot = factor(d2$prot, levels = porder)
p = 
  ggplot(d2 %>% filter(prot %in% tcp & variable == 'prot_library_size' & lineage == "T_Cell"),
         aes(x = prot, y = value*100, fill = norm )) +
  theme_bw() + 
  ggtitle(label = 'within T cell model') + 
  theme(title = element_text(size = 8)) + 
  ylab("% of variance explained\n by sequencing depth") + 
  geom_col(position = 'dodge') +
  theme(axis.text.y = element_text(color = 'black')) + 
  theme(axis.text.x = element_text(color = 'black')) + 
  theme(legend.position = 'none') +
  theme(axis.title.y = element_blank())  +
  scale_fill_manual(values = c('grey60', 'deepskyblue3')) + 
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, 'Tcell_var_explained.pdf'), width = 2, height = 3.3)


# global model (not as applicable as variation in lib size expected between lineage / celltypes)
dsb_ct = lapply(dl, function(x){ 
  fitExtractVarPartModel(exprObj = x, formula = f1, data = met, 
                         REML = FALSE, useWeights = FALSE)  })

dsb_var = variancePartition::fitExtractVarPartModel(exprObj = d1, formula = f1, data = met, REML = FALSE,useWeights = FALSE)
clr_var = variancePartition::fitExtractVarPartModel(exprObj = d2, formula = f1, data = met, REML = FALSE,useWeights = FALSE)

ddsb = dsb_var
dclr = clr_var

prot_combined = c(tcp, bcp, mcp) %>% unique 

rownames(ddsb) = str_sub(rownames(ddsb), 1, -6)
rownames(dclr) = str_sub(rownames(dclr), 1, -6)

saveRDS(ddsb, paste0(datapath,'ddsb_global_varpartmodel.rds'))
saveRDS(dclr, paste0(datapath,'dclr_global_varpartmodel.rds'))

p = plotPercentBars( ddsb[prot_combined, ] ) + ggtitle('dsb')
ggsave(p,filename = paste0(figpath, 'dsb_global_varpart.pdf'), width = 5, height = 5)
p = plotPercentBars( dclr[prot_combined, ] ) + ggtitle('clr') 
ggsave(p,filename = paste0(figpath, 'clr_global_varpart.pdf'), width = 5, height = 5)

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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] doParallel_1.0.16        iterators_1.0.13         foreach_1.5.1            magrittr_2.0.1          
# [5] SeuratObject_4.0.0       Seurat_4.0.1             mclust_5.4.7             forcats_0.5.1           
# [9] stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4              readr_1.4.0             
# [13] tidyr_1.1.2              tibble_3.0.6             tidyverse_1.3.0          here_1.0.1              
# [17] variancePartition_1.20.0 Biobase_2.50.0           BiocGenerics_0.36.1      scales_1.1.1            
# [21] BiocParallel_1.24.1      limma_3.46.0             ggplot2_3.3.3           
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1          backports_1.2.1       plyr_1.8.6            igraph_1.2.6          lazyeval_0.2.2       
# [6] splines_4.0.5         listenv_0.8.0         scattermore_0.7       digest_0.6.27         htmltools_0.5.1.1    
# [11] viridis_0.5.1         tensor_1.5            cluster_2.1.1         ROCR_1.0-11           globals_0.14.0       
# [16] modelr_0.1.8          matrixStats_0.58.0    spatstat.sparse_2.0-0 prettyunits_1.1.1     colorspace_2.0-0     
# [21] rvest_0.3.6           ggrepel_0.9.1         haven_2.3.1           crayon_1.4.1          jsonlite_1.7.2       
# [26] lme4_1.1-26           spatstat.data_2.1-0   survival_3.2-10       zoo_1.8-8             glue_1.4.2           
# [31] polyclip_1.10-0       gtable_0.3.0          leiden_0.3.7          future.apply_1.7.0    abind_1.4-5          
# [36] pheatmap_1.0.12       DBI_1.1.1             miniUI_0.1.1.1        Rcpp_1.0.6            viridisLite_0.3.0    
# [41] xtable_1.8-4          progress_1.2.2        reticulate_1.18       spatstat.core_2.0-0   htmlwidgets_1.5.3    
# [46] httr_1.4.2            gplots_3.1.1          RColorBrewer_1.1-2    ellipsis_0.3.1        ica_1.0-2            
# [51] farver_2.0.3          pkgconfig_2.0.3       uwot_0.1.10           deldir_0.2-10         dbplyr_2.1.0         
# [56] labeling_0.4.2        tidyselect_1.1.0      rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1        
# [61] munsell_0.5.0         cellranger_1.1.0      tools_4.0.5           cli_2.3.0             generics_0.1.0       
# [66] broom_0.7.5           ggridges_0.5.3        fastmap_1.1.0         goftest_1.2-2         fs_1.5.0             
# [71] fitdistrplus_1.1-3    caTools_1.18.1        RANN_2.6.1            pbapply_1.4-3         future_1.21.0        
# [76] nlme_3.1-152          mime_0.10             xml2_1.3.2            compiler_4.0.5        pbkrtest_0.5-0.1     
# [81] rstudioapi_0.13       plotly_4.9.3          png_0.1-7             spatstat.utils_2.1-0  reprex_1.0.0         
# [86] statmod_1.4.35        stringi_1.5.3         lattice_0.20-41       Matrix_1.3-2          nloptr_1.2.2.2       
# [91] ggsci_2.9             vctrs_0.3.6           pillar_1.4.7          lifecycle_1.0.0       BiocManager_1.30.10  
# [96] spatstat.geom_2.0-1   lmtest_0.9-38         RcppAnnoy_0.0.18      data.table_1.14.0     cowplot_1.1.1        
# [101] bitops_1.0-6          irlba_2.3.3           httpuv_1.5.5          patchwork_1.1.1       colorRamps_2.3       
# [106] R6_2.5.0              promises_1.2.0.1      KernSmooth_2.23-18    gridExtra_2.3         parallelly_1.23.0    
# [111] codetools_0.2-18      boot_1.3-27           MASS_7.3-53.1         gtools_3.8.2          assertthat_0.2.1     
# [116] rprojroot_2.0.2       withr_2.4.1           sctransform_0.3.2     mgcv_1.8-34           hms_1.0.0            
# [121] rpart_4.1-15          grid_4.0.5            minqa_1.2.4           Rtsne_0.15            shiny_1.6.0          
# [126] lubridate_1.7.9.2  
