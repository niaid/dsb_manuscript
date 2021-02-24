suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))

# set save path 
figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)

# load raw data 
h1 = readRDS(here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds"))
neg = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# separate by batch 
negadt1 = neg[[1]]; negadt2 = neg[[2]]
b1cells = h1@meta.data %>% filter(batch == "1") %$% barcode_check
b2cells = h1@meta.data %>% filter(batch == "2") %$% barcode_check
adt = h1@assay$CITE@raw.data %>% as.matrix()
adt1 = adt[ ,b1cells]
adt2 = adt[ ,b2cells]

# DSB normalize in steps over the 2 batches 
pseudocount.use = 10
adtu_log1 = log(negadt1 + pseudocount.use) 
adt_log1 = log(adt1 + pseudocount.use)
adtu_log2 = log(negadt2 + pseudocount.use) 
adt_log2 = log(adt2 + pseudocount.use) 

# batch 1 rescaling 
mu_u1 = apply(adtu_log1, 1 , mean)
sd_u1 = apply(adtu_log1, 1 , sd)
norm_adt1 = apply(adt_log1, 2, function(x) (x  - mu_u1) / sd_u1) 

# batch 2 rescaling 
mu_u2 = apply(adtu_log2, 1 , mean)
sd_u2 = apply(adtu_log2, 1 , sd)
norm_adt2 = apply(adt_log2, 2, function(x) (x  - mu_u2) / sd_u2) 

# merge adt normalized by batch 
norm_adt = cbind(norm_adt1, norm_adt2)

# save raw norm dsb data 
saveRDS(norm_adt, file = paste0(datapath, "dsb_norm_adt_mtx.rds"))
saveRDS(norm_adt1, file = paste0(datapath, "b1_dsb_norm_adt_mtx.rds"))
saveRDS(norm_adt2, file = paste0(datapath, "b2_dsb_norm_adt_mtx.rds"))

# visualize the relationship between the per cell technical component and protein library zise 
# run per-cell gaussian mixture model (on each batch)
library(mclust)
cellwise_model1 = apply(norm_adt1, 2, function(x) {
			g = Mclust(x, G=2, warn = TRUE , verbose = TRUE)  
			return(g) 
		})
cellwise_model2 = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = TRUE)  
  return(g) 
})

# tidy model fit data 
cm = c(cellwise_model1, cellwise_model2)
tm  = lapply(cm, function(x){broom::tidy(x)[ 2, 5:6]}) %>% bind_rows()
mr  = lapply(cm, broom::glance) %>% bind_rows()
tm$barcode_check = names(cm); mr$barcode_check = names(cm)
md1 = full_join(tm, mr, by = "barcode_check")

# merge model results with metadata 
md = full_join(md1,h1@meta.data,  by = "barcode_check")

# calculate latent component noise vector 
cellwise_background_mean = lapply(cm, function(x) {x$parameters$mean[1] })
cellwise_background_mean = unlist(cellwise_background_mean, use.names = FALSE)
cellwise_positive_mean = lapply(cm, function(x) {x$parameters$mean[2] })
cellwise_positive_mean = unlist(cellwise_positive_mean, use.names = FALSE)

# define pc1 through isotypes and background protein as a latent variable 
isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
                             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT" )
noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
get_noise_vector = function(noise_matrix) { 
  g = prcomp(t(noise_matrix), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
df = full_join(md, PC, by = "barcode_check")

# add library size 
nUMI_Prot = h1@assay$CITE@raw.data %>% colSums()
df$nUMI_Prot = nUMI_Prot

# add isotype control means 
isotypes = c("MouseIgG1kappaisotype_PROT","MouseIgG2akappaisotype_PROT", 
             "Mouse IgG2bkIsotype_PROT", "RatIgG2bkIsotype_PROT")
iso = norm_adt[ isotypes,  ]
iso_mean = apply(iso, 2, mean)

# mean comparison 
df = cbind(iso_mean, df)
saveRDS(df ,file = paste0(datapath, "mixture_model_metadata_merged.rds"))
saveRDS(iso, file = paste0(datapath, "isotype_values_dsb.rds"))
