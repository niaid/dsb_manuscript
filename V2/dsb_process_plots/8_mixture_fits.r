# Evaluate per cell Gaussian mixture model 
set.seed(1)
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(mclust))
'%ni%' = Negate('%in%')

# set save paths 
figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)

# load normalized B1 data and meta data
b1met = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")@meta.data %>% filter(batch == "1")
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))

##### mixture model parameter evaluation 
# Fit a gaussian mixture with g =  1 : 6 mixing components to evaluate G that globally maximizes BIC
cmd = list()
for (i in 1:6) {
  cmd[[i]]  = apply(norm_adt1, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = TRUE)  
    return(g) ; gc()
  })
}

# tidy mixture model parameter estimates for each cell 
mr = list()
for (i in 1:length(cmd)) {
  mr[[i]]  = lapply(cmd[[i]], broom::glance)
  mr[[i]] = do.call(rbind, mr[[i]])
  mr[[i]]$barcode_check = colnames(norm_adt1)
  mr[[i]] = full_join(mr[[i]], b1met, by = "barcode_check")
}
mrdf = do.call(rbind, mr)
saveRDS(mrdf, file = paste0(datapath, "multi_component_model_comparison.rds"))
mrdf = readRDS(file = here('V2/dsb_process_plots/generated_data/multi_component_model_comparison.rds'))

##### plot an arbitrary single fit (the first cell in the matrix) from the 2 component mixture. 
cell = cmd[[2]][[1]]
cell = MclustDR(cell)
pdf(file = paste0(figpath, "examplecell.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()
#### 

# compare mixture g = 1,3,4,5,6 to G = 2 difference in BIC
test = mrdf %>% select(BIC, G,barcode_check, p3_dist_3) %>% spread(G, BIC) 
t2 = test %>% as.data.frame()
colnames(t2) = c('barcode_check', 'p3_dist_3', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6')
ls = list()
mixing_components = c("G1",  "G3", "G4", "G5", "G6")
for (i in 1:length(mixing_components)) {
  ls[[i]] =  t2 %>% filter(G2 < t2[[mixing_components[i]]] ) 
  ls[[i]]$bic_min = mixing_components[i]
}
saveRDS(ls, file = paste0(datapath, "cell_BIC_subset_with_G2_notoptimal.rds"))
ls = readRDS(file = here("V2/dsb_process_plots/generated_data/cell_BIC_subset_with_G2_notoptimal.rds"))
# filter to the selected BIC 
test = mrdf %>% 
  select(BIC, G,barcode_check, p3_dist_3) %>% 
  as.data.frame() %>% 
  arrange(barcode_check) %>% 
  group_by(barcode_check) %>% 
  mutate(bic_ = max(BIC)) %>% 
  filter(bic_ == BIC)
test_best = test %>% select(barcode_check, G_BEST=  G)

# pct of cells with g = 2 best fit. 
(test_best %>% filter(G_BEST == 2) %>% nrow) / nrow(test_best)
# 0.8162882

(test_best %>% filter(G_BEST == 4) %>% nrow) / nrow(test_best)

#cellwise_background_mean = cellwise_positive_mean = list()
dmean = list()
for (i in 1:length(cmd)) {
  dmean[[i]] = lapply(cmd[[i]], function(x) {x$parameters$mean }) 
  dmean[[i]] = do.call(rbind, dmean[[i]])
  colnames(dmean[[i]]) = paste0("G_",i,"_mu_",colnames(dmean[[i]]))
}
dmean_df = do.call(cbind, dmean) %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode_check")
saveRDS(dmean_df, file = paste0(figpath, "dmean_df.rds"))
dmean_df = readRDS(file = here("V2/dsb_process_plots/figures/dmean_df.rds"))

# add BIC values 
BICdf = mrdf %>% 
  select(BIC, G,barcode_check) %>%
  spread(G, BIC) %>% 
  as.data.frame()
colnames(BICdf) = c("barcode_check", "BIC_G1", "BIC_G2", "BIC_G3", "BIC_G4", "BIC_G5", "BIC_G6")
dmean_df = full_join(dmean_df, BICdf)

# add the optimal BIC 
dmean_test = full_join(dmean_df, test_best, by = "barcode_check")

###############################
# figure generation 

# plot BIC distribution 
g3best = dmean_test %>% filter(G_BEST == 3) %>% 
  mutate(g3_g2_delta = BIC_G3 - BIC_G2) %>%
  mutate(g2_g1_delta = BIC_G2 - BIC_G1)
p1 = ggplot(g3best, aes(x = g3_g2_delta)) + 
  theme_bw() + 
  xlim(c(0,200)) + 
  geom_histogram(bins = 30, fill = "grey", color = "black") + 
  xlab("k=3 - k=2 model") + 
  theme(axis.title.x = element_text(size = 18))
p2 = ggplot(g3best, aes(x = g2_g1_delta )) + 
  theme_bw() + 
  geom_histogram(bins = 30, fill = "grey", color = 'black') + 
  xlab("k=2 - k=1 model") + 
  theme(axis.title.x = element_text(size = 18))
p3 = cowplot::plot_grid(p1,p2, ncol = 1) 
ggsave(p3, filename = paste0(figpath, 'BIC_DELTA_cells_g3best.pdf'), width = 2.5, height = 5)

# g3 sub models 
dsub = dmean_test %>% filter(G_BEST == 3) %>% 
  select(G_2_mu_1,G_2_mu_2, G_3_mu_1 ,  G_3_mu_2 , G_3_mu_3) %>% 
  gather(parameter, mean, G_2_mu_1:G_3_mu_3) %>% 
  mutate(model = str_sub(parameter, 1,3 ))
#plot distributions
p = ggplot(dsub %>% filter(mean < 15), aes(x = mean, fill = parameter)) + 
  theme_bw() + 
  geom_density() +
  theme(strip.background = element_blank()) + 
  facet_wrap(~model) + 
  ggsci::scale_fill_d3(alpha = 0.8)  
ggsave(p, filename = paste0(figpath, "cells_with_G3best_g3vg2.pdf"), width = 4.5, height = 2.5)

# g4 sub models 
dsub = dmean_test %>% filter(G_BEST == 4) %>% 
  select(G_2_mu_1,G_2_mu_2, G_4_mu_1 ,  G_4_mu_2 , G_4_mu_3, G_4_mu_4) %>% 
  gather(parameter, mean, G_2_mu_1:G_4_mu_4) %>% 
  mutate(model = str_sub(parameter, 1,3 ))
#plot distributions
p = ggplot(dsub %>% filter(mean < 15), aes(x = mean, fill = parameter)) + 
  theme_bw() + 
  geom_density() +
  theme(strip.background = element_blank()) + 
  facet_wrap(~model) + 
  ggsci::scale_fill_d3(alpha = 0.8)  
ggsave(p, filename = paste0(figpath, "cells_with_G4_best_g4vg2.pdf"), width = 4.5, height = 2.5)

# g5 sub models 
dsub = dmean_test %>% filter(G_BEST == 4) %>% 
  select(G_2_mu_1,G_2_mu_2, G_5_mu_1 ,  G_5_mu_2 , G_5_mu_3, G_5_mu_4, G_5_mu_5) %>% 
  gather(parameter, mean, G_2_mu_1:G_5_mu_5) %>% 
  mutate(model = str_sub(parameter, 1,3 ))
#plot distributions
p = ggplot(dsub %>% filter(mean < 15), aes(x = mean, fill = parameter)) + 
  theme_bw() + 
  geom_density() +
  theme(strip.background = element_blank()) + 
  facet_wrap(~model) + 
  ggsci::scale_fill_d3(alpha = 0.4)  + 
  ggtitle("Cells with G=5 best fit ")
ggsave(p, filename = paste0(figpath, "cells_with_G5_best_g5vg2.pdf"), width = 4.5, height = 2.5)

## best fit distribution 
ctmd = b1met %>% select(barcode_check, celltype_label_3, p3_dist_3)
merge = full_join(dmean_test, ctmd)
merge$G_BEST = as.character(merge$G_BEST)

p = ggplot(merge, aes(x = p3_dist_3, fill = G_BEST )) + 
  geom_bar(position = "fill") + ggsci::scale_fill_d3() + 
  theme_minimal()
ggsave(p ,filename = paste0(figpath, "maxbic_celtype.pdf"), width = 5, height = 2.1)

merge$dv = "dv"
p = ggplot(merge, aes(x = dv, fill = G_BEST )) + 
  geom_bar(position = "fill", show.legend = FALSE) + 
  ggsci::scale_fill_d3() + 
  ylab("proportion with optimal BIC") + 
  theme_minimal() + 
  scale_y_continuous(position = "right") + 
  theme(axis.text.x =  element_blank(), axis.title.x = element_blank())
ggsave(p ,filename = paste0(figpath, "maxbic_GLOBAL.pdf"), width = 1, height = 3) 

# plot example cell with BIC best G = 3 
cellsub = ls[[2]][1,1]
cell = MclustDR(cmd[[3]][[cellsub]])
pdf(file = paste0(figpath,cell, "G3.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()
# plot the same cells G = 2 distribution 
#cellsub = ls[[2]][1,1]
cell = MclustDR(cmd[[2]][[cellsub]])
pdf(file = paste0(figpath,cell, "G2.pdf"), width = 3, height = 4)
plot(cell, what = "contour")
dev.off()

## plot global fits showing G=2 tends to be optimal 
# Mixture model components vs BIC plot 
mrdf = readRDS("V2/dsb_process_plots/generated_data/multi_component_model_comparison.rds")
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3", "4", "5", "6"))
p = ggplot(mrdf , aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,  show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  facet_wrap(~ p3_dist_3, nrow = 2) + 
  xlab(" number of mixing components in Gaussian mixture model ") + 
  ylab(" BIC ") + 
  theme(axis.title.x =element_text(size = 18)) + 
  theme(axis.title.y =element_text(size = 18)) + 
  theme(strip.text = element_text(size = 14))
ggsave(p, filename = paste0(figpath, "mixture_component_vs_bic.pdf"), width = 15 ,height = 5)  

# same plot, show only global distribution 
p = ggplot(mrdf , aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,  show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 10)) + 
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 14)) 
p
ggsave(p, filename = paste0(figpath, "GLOBALmixture_component_vs_bic.pdf"), width = 2 ,height = 3)  


## mu 1 for G = 2 vs mu1 for G = 3 
figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)

# load normalized B1 data and meta data
b1met = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")@meta.data %>% filter(batch == "1")
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))

##### mixture model parameter evaluation 
# Fit a gaussian mixture with g =  1 : 6 mixing components to evaluate G that globally maximizes BIC
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))
cmd = list()
for (i in 2:3) {
  cmd[[i]]  = apply(norm_adt1, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = TRUE)  
    return(g) ; gc()
  })
}
mu1_g2 = lapply(cmd[[2]], function(x){ x$parameters$mean[1] }) %>% unlist()
mu1_g3 = lapply(cmd[[3]], function(x){ x$parameters$mean[1] }) %>% unlist()
plot(mu1_g2, mu1_g3)
d = data.frame( mu1_g2, mu1_g3)
saveRDS(d, file = paste0(figpath,"d_g3_vs_g3_mu1.rds"))
d = readRDS(file = here('V2/dsb_process_plots/figures/d_g3_vs_g3_mu1.rds'))
p =ggplot(d, aes(x = mu1_g2, mu1_g3)) +
  theme_bw() + 
  geom_bin2d(bins = 400) + 
  scale_fill_viridis_c(option = "B") +
  ggpubr::stat_cor() + 
  theme(legend.position = c(0.8, 0.4))   + 
  xlab("µ1 2-component mixture") + 
  ylab("µ1 3-component mixture")
ggsave(p,filename = paste0(figpath, '3_vs_2_mu1.pdf'), width = 2.1, height = 3.3)


####################
# Normalization of g3 BIC best cells with µ1 from g=2 and µ1 from g=3 models 
d = readRDS(file = here('V2/dsb_process_plots/figures/d_g3_vs_g3_mu1.rds'))
rownames(d) = str_sub(rownames(d), 1,-3)

# load raw data 
h1 = readRDS(here("data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds"))
neg = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))

# dsb norm the g3 bic optimal cells 
dlog = log(h1@assay$CITE@raw.data + 10)[ ,rownames(d)]
nlog = log(neg[[1]] + 10)[ ,1:1000]

# process with dsb 
sd_nlog = apply(nlog, 1 , sd)
mean_nlog = apply(nlog, 1 , mean)
norm_adt = apply(dlog, 2, function(x) (x  - mean_nlog) / sd_nlog) 

# define isotype controls 
isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
                             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT" )
# construct noise matrix from both µ1 parameters 
noise_matrix_g2 = t(rbind(norm_adt[isotype.control.name.vec, ], mu1g2 = d$mu1_g2))
noise_matrix_g3 = t(rbind(norm_adt[isotype.control.name.vec, ], mu1g3 = d$mu1_g3))

# pc1 regression 
g2pc1 = prcomp(noise_matrix_g2, scale = TRUE)$x[ ,1]
g3pc1 = prcomp(noise_matrix_g3, scale = TRUE)$x[ ,1]
denoised_adt1 = limma::removeBatchEffect(norm_adt, covariates = g2pc1)
denoised_adt2 = limma::removeBatchEffect(norm_adt, covariates = g3pc1)

cormat = cor(t(denoised_adt1), t(denoised_adt2), method = 'pearson')
pdf(file = paste0(figpath,'dsb_g3_vs_g2_denoising.pdf'), width = 4, height = 4)
plot(diag(cormat),ylim = c(0.50, 1), main = 'pearson correlation')
dev.off()

range(diag(cormat))
# 0.9536114 0.9999648
mean(diag(cormat))
# 0.9918949
median(diag(cormat))
# 0.9949126

d2 = data.frame(g2pc1, g3pc1)

p = 
  ggplot(d2, aes(x = g2pc1, y = g3pc1)) + 
  theme_bw() +
  geom_bin2d(bins = 100) +
  scale_fill_viridis_c(option = "B") + 
  geom_abline(slope = 1, linetype = 'dashed') + 
  xlab('isotype controls combined with \n 2 component mixture µ1') +  
  ylab('isotype controls combined with \n 3 component mixture µ1') + 
  theme(axis.title = element_text(size =14)) + 
  ggpubr::stat_cor(aes(label = ..r.label..)) + 
  theme(legend.position = c(0.8,0.25))
p
ggsave(p, filename = paste0(figpath,'tech_component_g2_vs_g3.pdf'), width = 3.5, height = 4.3)





