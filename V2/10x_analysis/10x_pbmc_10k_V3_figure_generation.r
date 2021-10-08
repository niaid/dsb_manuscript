# 10x PBMC 10k data 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
"%ni%" = Negate("%in%")

# set params 
project_title = "10X PBMC 10K V3"
prot_plot = c("CD19_TotalSeqB", "CD3_TotalSeqB", "CD14_TotalSeqB",  "CD4_TotalSeqB" , "CD56_TotalSeqB", "CD8a_TotalSeqB")
neg_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC 10K V3neg_prot2.rds"))
pos_cells = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC 10K V3pos_prot.rds"))
df_dsb = readRDS(file = here("V2/10x_analysis/generated_data/10X PBMC 10K V3dsb_merged_result.RDS"))

# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/")
datapath = here("V2/10x_analysis/generated_data/")

## make tidy dataframe for plot
index1 = colnames(df_dsb)[14]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# cluster umap plots ; calc centers in umap space 
centers = df_dsb %>% 
  dplyr::group_by(clusters) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# labeled umap 
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank()) +  
  geom_point(mapping = aes(color = clusters), size = 0.5, shape = 16, alpha = 0.3, show.legend = FALSE) + 
  ggsci::scale_color_d3(palette = "category20") + 
  ggrepel::geom_text_repel(data = centers, box.padding = 0.5,
                           mapping = aes(label = clusters, size = 2.3, fontface = "bold"),
                           show.legend = FALSE) 
ggsave(p, filename = paste0(figpath, project_title, "clusters.png"), width = 3.3, height = 3.3)

# protein distributions 
plotsub = dsb %>% filter(DSB > -5) %>% filter(DSB < 40) 
plotsub_spread = plotsub %>% spread(prot, DSB, drop = TRUE)
plotsub_spread =  na.omit(plotsub_spread)
p = ggplot(plotsub %>% filter(prot %in% prot_plot), aes(x = UMAP_1, y = UMAP_2, color = DSB)) +
  geom_point(shape = 16) + 
  theme(legend.key.size = unit(0.8, units = "cm"), legend.title = element_text(size = 18, face = "bold"), legend.text = element_text(size = 20, face = "bold")) + 
  scale_color_viridis_c(option = "B") +
  facet_wrap(~ prot,nrow = 2) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15, face = "bold")) + 
  theme(legend.position = "right")
ggsave(p, filename = paste0(figpath, project_title, "clusters_DSBdist.png"), width = 12, height = 7.4)

# biaxial plots
mg_theme = list(
  theme_bw(),
  theme(axis.title.x =element_text(size = 18),axis.title.y = element_text(size = 18)), 
  geom_bin2d(bins = 200, show.legend = FALSE),
  viridis::scale_fill_viridis(option = "B") 
)

# plot manual gate distributions
p = ggplot(plotsub_spread, aes(x = plotsub_spread[[prot_plot[1]]], y = plotsub_spread[[prot_plot[2]]]))  + xlab("CD19") + ylab("CD3") + mg_theme 
ggsave(p, filename = paste0(figpath, project_title, "3_19_mg.pdf"), width = 3, height = 3)  
p = ggplot(plotsub_spread, aes(x = plotsub_spread[[prot_plot[4]]], y = plotsub_spread[[prot_plot[3]]]))  + xlab("CD4") + ylab("CD14") + mg_theme 
ggsave(p, filename = paste0(figpath, project_title, "4_14_mg.pdf"), width = 3, height = 3)  
p = ggplot(plotsub_spread, aes(x = plotsub_spread[[prot_plot[4]]], y = plotsub_spread[[prot_plot[6]]]))  + xlab("CD4") + ylab("CD8") + mg_theme 
ggsave(p, filename = paste0(figpath, project_title, "4_8_mg.pdf"), width = 3, height = 3)  

# protein heatmap by clusters 
prots = dsb$prot %>% unique 
adt_plot = df_dsb %>% 
  group_by(clusters) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("clusters") %>% 
  t %>% 
  as.data.frame

# average heatmap
pheatmap::pheatmap(adt_plot, color = viridis::viridis(18, option = "B"), 
                   fontsize_row = 9, border_color = NA,
                   treeheight_row = 10, treeheight_col = 10,
                   filename = paste0(figpath, project_title, "average_cluster_dsb_heatmap.pdf"),
                   width = 4.5, height = 5)


##########################
##DSB process plots 
pseudocount.use = 10
adtu_log = log(neg_cells + pseudocount.use) 
adt_log = log(pos_cells + pseudocount.use)

# background rescale counts 
mu_u = apply(adtu_log, 1 , mean)
sd_u = apply(adtu_log, 1 , sd)
norm_adt2 = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u) 

# gaussian mixture on background rescaled matrix 
library(mclust)
cellwise_model = apply(norm_adt2, 2, function(x) {
  g = Mclust(x, G=2, warn = TRUE , verbose = FALSE)  
  return(g) 
})

# tidy Gaussian mixture model data 
tm  = lapply(cellwise_model, function(x){broom::tidy(x)[ 2, 5:6]})
tm = do.call(rbind, tm)
tm$barcode_check = colnames(adt_log)
mr  = lapply(cellwise_model, broom::glance)
mr = do.call(rbind, mr)
mr$barcode_check = colnames(adt_log)

# merge model results with metadata 
df_dsb$barcode_check = rownames(df_dsb)
md = full_join(df_dsb, mr, by = "barcode_check")
md = full_join(md, tm, by = "barcode_check")

################
# calculate latent component noise vector 
cellwise_background_mean = lapply(cellwise_model, function(x) {x$parameters$mean[1] }) %>% unlist(use.names = FALSE)

# define latent noise component 
isotypes = rownames(adt_log)[15:17]; isotypes 
noise_matrix = rbind(norm_adt2[isotypes, ], cellwise_background_mean)
get_noise_vector = function(noise_matrix) { 
  g = prcomp(t(noise_matrix), scale = TRUE)
  return(g$x[ ,1]) 
} 
noise_vector = get_noise_vector(noise_matrix)

# add noise vector to cellwise mixture model 
PC = as.data.frame(noise_vector) %>% rownames_to_column("barcode_check")
md = full_join(md, PC, by = "barcode_check")

# noise vector vs protein library size clusters 
p = ggplot(md, aes(x = log10umiprot, y = noise_vector)) + 
  theme_bw() + 
  ggpubr::stat_cor(method = 'pearson', aes(label = ..r.label..)) + 
  geom_bin2d(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(strip.background = element_blank()) + 
  geom_smooth(color = "#3e8ede", method = 'lm') + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) + 
  facet_wrap(~clusters, scales = 'free')
ggsave(p, filename = paste0(figpath, project_title, "ALTBINnoise_vector_vs_libsize.pdf") , width = 4.5, height = 4, )
# not faceted by celltype 
p = ggplot(md, aes(x = log10umiprot, y = noise_vector)) + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede") + 
  theme(axis.text.x = element_text(size = 5)) + 
  theme(strip.background = element_blank()) + 
  xlab("log10 prot library size") + 
  ylab("dsb Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title, "BINnoise_vector_vs_libsize.pdf") , width = 3.5, height = 3.5)

# isotype control vs mixture mean 1 
iso = norm_adt2[ isotypes,  ]
iso_mean = apply(iso, 2, mean)

# add isotype mean to metadata.  
md = cbind(md, iso_mean)

# plot 
p = ggplot(data = md, 
           aes(x = iso_mean, y = mean.1)) + 
  geom_bin2d(bins = 180, show.legend = FALSE) + 
  ggpubr::stat_cor(method = "pearson") + 
  scale_fill_viridis_c(option = "B") + 
  theme_bw() + 
  xlab("isotype means") + 
  ylab(" µ1 background mean ") + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title,"meanisovs_meanbackground.pdf"), width = 3, height = 3)

# correlation across latent component. 
isotype_value = iso %>% t %>% as.data.frame() %>% rownames_to_column("bc")
mean1df  = md %>% select(mean.1, bc = barcode_check)
dcor = full_join(isotype_value, mean1df, by = "bc") %>% column_to_rownames("bc")
dcor = dcor[complete.cases(dcor), ]
cplot = cor(dcor, method = "spearman")
pcor = Hmisc::rcorr(as.matrix(dcor), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "µ1")
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
pdf(file = paste0(figpath,project_title,"background_correlation_plot_spearman.pdf"), width = 4,height = 4)
corrplot::corrplot(cplot,method="color", col=col(200),  
                   cl.lim = c(0,1),
                   type="upper",
                   order = "original",
                   addCoef.col = "black", 
                   cl.pos = "b", 
                   tl.col="black", tl.srt=45, #
                   diag=FALSE, addgrid.col = "ghostwhite"
)
dev.off()
# Save cell metadata created 
data.table::fwrite(md, file = paste0(datapath, project_title, '_cellmetadata.txt'),sep = "\t")


#####################
# mixture model parameter evaluation 
cmd = list()
for (i in 1:3) {
  cmd[[i]]  = apply(norm_adt2, 2, function(x) {
    g = Mclust(x, G=i, warn = TRUE , verbose = FALSE)  
    return(g) ; gc()
  })
}
mr2 = list()
for (i in 1:length(cmd)) {
  mr2[[i]]  = lapply(cmd[[i]], broom::glance)
  mr2[[i]] = do.call(rbind, mr2[[i]])
  mr2[[i]]$barcode_check = colnames(adt_log)
}
mrdf = do.call(rbind, mr2)
cmd = df_dsb %>% select(barcode_check, clusters)
mrdf = full_join(mrdf, cmd, by = "barcode_check")

# save mixture fit dataframe 
write_delim(mrdf, path = paste0(datapath, project_title, "_mixture_model_comparison.txt"),delim = "\t")

# calculate percent of cells with G = 2 best fit 
test = mrdf %>% 
  select(BIC, G,barcode_check) %>% 
  as.data.frame() %>% 
  arrange(barcode_check) %>% 
  group_by(barcode_check) %>% 
  mutate(bic_ = max(BIC)) %>% 
  filter(bic_ == BIC)
test_best = test %>% select(barcode_check, G_BEST=  G)


# plot g vs BIC 
mrdf$G = factor(mrdf$G, levels = c("1", "2" ,"3"))
#table(is.na(mrdf$BIC)) # n = 1 cell 
#mrdf = mrdf[,3:4]
#mrdf = mrdf[complete.cases(mrdf), ]
p = ggplot(mrdf, aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.shape = NA, lwd = 0.2,   show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  xlab("Mixing Components") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 16)) + 
  theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title, "GLOBALmixture_model_BICvsG.pdf"), width = 2, height = 3)


###############################
# write stats table for s table.  
d = data.frame(
  dataset = project_title,
  nprot = nrow(pos_cells), 
  ncells = ncol(pos_cells), 
  n_background = ncol(neg_cells), 
  cor_tech_size = cor(md$log10umiprot, md$noise_vector), 
  cor_isotype_background = cor(md$iso_mean, md$mean.1), 
  cells_removed_for_visualization = dim(df_dsb)[1] - dim(plotsub_spread)[1]
)
data.table::fwrite(d, file = paste0(datapath, project_title, "TABLE_STATS.txt"), sep = "\t")

## 
mu1_g2 = lapply(cmd[[2]], function(x){ x$parameters$mean[1] }) %>% unlist()
mu1_g3 = lapply(cmd[[3]], function(x){ x$parameters$mean[1] }) %>% unlist()
plot(mu1_g2, mu1_g3)
d = data.frame( mu1_g2, mu1_g3)
saveRDS(d, file = paste0(datapath,"d_g3_vs_g3_mu1.rds"))

p =ggplot(d, aes(x = mu1_g2, mu1_g3)) +
  theme_bw() + 
  geom_bin2d(bins = 400) + 
  scale_fill_viridis_c(option = "B") +
  ggpubr::stat_cor() + 
  theme(legend.position = c(0.8, 0.4))   + 
  xlab("µ1 2-component mixture") + 
  ylab("µ1 3-component mixture")
p
ggsave(p,filename = paste0(figpath, '3_vs_2_mu1.pdf'), width = 2.1, height = 3.3)

d = readRDS(file = here('V2/dsb_process_plots/figures/d_g3_vs_g3_mu1.rds'))
rownames(d) = str_sub(rownames(d), 1,-3)

# load raw data 
pseudocount.use = 10
adtu_log = log(neg_cells + pseudocount.use) 
adt_log = log(pos_cells + pseudocount.use)

# background rescale counts 
mu_u = apply(adtu_log, 1 , mean)
sd_u = apply(adtu_log, 1 , sd)
norm_adt2 = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u) 

# define isotype controls 
isotype.control.name.vec = isotypes

# construct noise matrix from both µ1 parameters 
noise_matrix_g2 = t(rbind(norm_adt2[isotype.control.name.vec, ], mu1g2 = d$mu1_g2))
noise_matrix_g3 = t(rbind(norm_adt2[isotype.control.name.vec, ], mu1g3 = d$mu1_g3))

# pc1 regression 
g2pc1 = prcomp(noise_matrix_g2, scale = TRUE)$x[ ,1]
g3pc1 = prcomp(noise_matrix_g3, scale = TRUE)$x[ ,1]
denoised_adt1 = limma::removeBatchEffect(norm_adt2, covariates = g2pc1)
denoised_adt2 = limma::removeBatchEffect(norm_adt2, covariates = g3pc1)

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


test = mrdf %>% 
  select(BIC, G,barcode_check) %>% 
  as.data.frame() %>% 
  arrange(barcode_check) %>% 
  group_by(barcode_check) %>% 
  mutate(bic_ = max(BIC)) %>% 
  filter(bic_ == BIC)
test_best = test %>% select(barcode_check, G_BEST=  G)
test_best$G_BEST %>% table 






