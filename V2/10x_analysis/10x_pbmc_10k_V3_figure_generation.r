# 10x PBMC 10k data 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
library(here)
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
index1 = colnames(df_dsb)[10]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# cluster umap plots ; calc centers in umap space 
centers = df_dsb %>% 
  dplyr::group_by(clusters) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# labeled umap 
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2)) + 
  theme_minimal() + 
  geom_point(mapping = aes(color = clusters),  shape = 16, alpha = 0.3, show.legend = FALSE) + 
  ggsci::scale_color_d3(palette = "category20") + 
  ggrepel::geom_text_repel(data = centers, box.padding = 0.5,
                           mapping = aes(label = clusters, size = 2.3, fontface = "bold"),
                           show.legend = FALSE) 
ggsave(p, filename = paste0(figpath, project_title, "clusters.png"), width = 3.3, height = 3.4)


# protein distributions 
plotsub = dsb %>% filter(DSB > -5) %>% filter(DSB < 40) 
plotsub_spread = plotsub %>% spread(prot, DSB, drop = TRUE)
plotsub_spread =  na.omit(plotsub_spread)
outliers_removed = dim(df_dsb)[1] - dim(plotsub_spread)[1]

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
  theme(axis.title.x =element_text(size = 18),
        axis.title.y = element_text(size = 18)), 
  geom_bin2d(bins = 250, show.legend = FALSE),
  viridis::scale_fill_viridis()
)

# plot manual gate distributions
# prot_plot = c("CD19_TotalSeqB", "CD3_TotalSeqB", "CD14_TotalSeqB",  "CD4_TotalSeqB" , "CD56_TotalSeqB", "CD8a_TotalSeqB")
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
pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"), fontsize_row = 8,border_color = NA,
                   filename = paste0(figpath, project_title, "average_cluster_dsb_heatmap.pdf"), width = 5, height = 5)

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
cellwise_background_mean = lapply(cellwise_model, function(x) {x$parameters$mean[1] })
cellwise_background_mean = unlist(cellwise_background_mean, use.names = FALSE)

# define latent noise component 
isotypes = rownames(adt_log)[15:17]
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
prot_lib_size = colSums(adt_log) %>% as.data.frame() 
colnames(prot_lib_size) = "nUMI_Prot"
prot_lib_size = prot_lib_size %>% rownames_to_column("barcode_check")
md = full_join(md, prot_lib_size, by = "barcode_check")


p = ggplot(md, aes(x = log10(nUMI_Prot), y = noise_vector)) + 
  theme_bw() + 
  geom_bin2d(bins = 200, show.legend = FALSE) + 
  scale_fill_viridis_c() + 
  geom_smooth(color = "#3e8ede") + 
  xlab("log10 prot library size") + 
  ylab("Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title, "BINnoise_vector_vs_libsize.pdf") , width = 3, height = 3)

# isotype control vs mixture mean 1 
iso = norm_adt2[ isotypes,  ]
iso_mean = apply(iso, 2, mean)

# mean comparison 
ndf = rbind(iso_mean, cellwise_background_mean) %>% t %>% as.data.frame()
p = ggplot(data = ndf %>% filter(iso_mean < 20 & cellwise_background_mean<12), 
           aes(x = iso_mean, y = cellwise_background_mean)) + 
  geom_bin2d(bins = 180, show.legend = FALSE) + 
  ggpubr::stat_cor(method = "pearson") + 
  scale_fill_viridis_c() + 
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
table(is.na(mrdf$BIC)) # n = 1 cell 
mrdf = mrdf[,3:4]
mrdf = mrdf[complete.cases(mrdf), ]
p = ggplot(mrdf, aes(x = G, y = BIC, fill = G)) +
  theme_bw() + 
  geom_boxplot(outlier.size = 0, outlier.alpha = 0.1, lwd = 0.2,  outlier.color = "grey38", show.legend = FALSE) +
  ggsci::scale_fill_d3() + 
  theme(strip.background = element_blank()) + 
  xlab(" Mixing Components ") + 
  ylab("Model  BIC ") +
  theme(axis.text.x = element_text(size = 16)) + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, project_title, "GLOBALmixture_model_BICvsG.pdf"), width = 3, height = 3)

###############################
# write stats table 
neg_c = neg_cells %>% t %>% as.data.frame() 
negmean = apply(neg_c, 2, mean) %>% mean

pos_c = pos_cells %>% t %>% as.data.frame() 
posmean = apply(pos_c, 2, mean) %>% mean

nprot = nrow(pos_cells)
ncells = ncol(pos_cells)
n_background = ncol(neg_cells)

stats10x = round(c(nprot, ncells, posmean, n_background, negmean), 0) 
names(stats10x) = c("n_prot",	"n_cells",	"prot. Avg.",	"n_drops",	"neg prot. Avg.")

stats10x = t(as.data.frame(stats10x)) %>% as.data.frame()

# outliers removed for umap visualization 
stats10x$outliers_removed_for_Vis = outliers_removed
nrm_isoplot = dim(ndf)[1] - (dim(ndf %>% filter(iso_mean < 20 & cellwise_background_mean<12))[1])

# outliers removed for correlation visualization
stats10x$outlier_rm_statcor = nrm_isoplot

# pct of cells with g = 2 best fit. 
pct2best = (test_best %>% filter(G_BEST == 2) %>% nrow) / nrow(test_best)
stats10x$percentg2best = pct2best

# save stats table 
write_delim(stats10x, path = paste0(datapath, project_title, "_stats.txt"), delim = "\t")
