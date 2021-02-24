# Evaluate per cell Gaussian mixture model 
set.seed(1)
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(mclust))
# define a custom infix for 'not in' 
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
  xlim(c(0,45)) + 
  geom_histogram(bins = 200) + 
  ggtitle("cells with G = 3 best fit \n BIC G3 - BIC G2")
p2 = ggplot(g3best, aes(x = g2_g1_delta )) + 
  theme_bw() + 
  geom_histogram(bins = 200) + 
  ggtitle("cells with G = 3 best fit \n BIC G2 - BIC G1")
p3 = cowplot::plot_grid(p1,p2)
ggsave(p3, filename = paste0(figpath, 'BIC_DELTA_cells_g3best.pdf'), width = 7, height = 3.5)

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
  geom_bar(position = "fill", show.legend = FALSE) + ggsci::scale_fill_d3() + 
  theme_minimal() + 
  theme(axis.text.x =  element_blank(), axis.title.x = element_blank())
ggsave(p ,filename = paste0(figpath, "maxbic_GLOBAL.pdf"), width = 2, height = 3.8) 

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
  geom_boxplot(outlier.size = 0, outlier.alpha = 0.1, lwd = 0.2,  outlier.color = "grey38", show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  facet_wrap(~ p3_dist_3, nrow = 2) + 
  xlab(" number of mixing components in Gaussian mixture model ") + 
  ylab(" BIC ")
ggsave(p, filename = paste0(figpath, "mixture_component_vs_bic.pdf"), width = 14 ,height = 5)  

# same plot, show only global distribution 
p = ggplot(mrdf , aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.size = 0, outlier.alpha = 0.1, lwd = 0.2,  outlier.color = "grey38", show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 16)) + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, "GLOBALmixture_component_vs_bic.pdf"), width = 3 ,height = 3)  

