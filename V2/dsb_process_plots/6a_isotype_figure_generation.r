suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
library(here)
setwd(here())

figpath = here("V2/dsb_process_plots/figures/"); dir.create(figpath)
datapath = here("V2/dsb_process_plots/generated_data/"); dir.create(datapath)

#######################################
# figure generation 
# read processed data generated above
df = readRDS(file = here("V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds"))
norm_adt = readRDS(file = here('V2/dsb_process_plots/generated_data/dsb_norm_adt_mtx.rds'))
iso  = readRDS(file = here("V2/dsb_process_plots/generated_data/isotype_values_dsb.rds"))


# plot noise vector vs lib size 
p = ggplot(df, aes(x = log10(nUMI_Prot), y = -1*(noise_vector))) + 
  theme_bw() + 
  geom_bin2d(bins = 200, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede") + 
  xlab("log10 prot library size") + 
  ylab("Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, "noise_vector_vs_libsize.pdf"), width = 8 ,height = 8.5)  


# facet by main lineage 
celltypes = df$celltype_label_1 %>% unique() %>% sort()
tcell = celltypes[c(2,3,4,5,10)]
myeloid = celltypes[c(6,8,9)]
bcell = celltypes[c(1)]
nk = celltypes[c(7)]
plot_sub = df %>%
  mutate(lineage = 
  if_else(celltype_label_1 %in% tcell, "T Cells",
  if_else(celltype_label_1 %in% myeloid, "Myeloid",
  if_else(celltype_label_1 %in% bcell, "B Cells", false = "NK")))) 
plot_sub$lineage = factor(plot_sub$lineage, levels = c("T Cells","Myeloid","B Cells","NK" ))
p = ggplot(plot_sub, aes(x = log10(nUMI_Prot), y = -1*(noise_vector))) + 
  facet_wrap(~lineage, scales = "free", ncol = 2) + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  theme(legend.position = c(0.1, 0.85), legend.key.size = unit(0.3, units = "cm")) + 
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede") + 
  theme(strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 20)) + 
  xlab("log10 protein library size") + 
  ylab("Technical Component")  + 
  theme(axis.title.x = element_text(size = 19)) +
  theme(axis.title.y = element_text(size = 19)) 
ggsave(p, filename = paste0(figpath, "lineage_noise_vector_vs_libsize.pdf"), width = 6, height = 5.3)  


# for a supplemental figure plot the noise vector by each p3 dist cell type 
p = ggplot(df, aes(x = log10(nUMI_Prot), y = -1*(noise_vector))) + 
  facet_wrap(~p3_dist_3, scales = "free", nrow = 2) + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  theme(legend.position = "bottom") +
  scale_fill_viridis_c(option = "B") + 
  geom_smooth(color = "#3e8ede") + 
  theme(strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 10)) + 
  xlab("log10 prot library size") + 
  ylab("Denoising Covariate")  + 
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4)) 
ggsave(p, filename = paste0(figpath, "celltype_noise_vector_vs_libsize.pdf"), width = 15.5, height = 5.5)  


# protin library size vs µ1
p1 = ggplot(df, aes(x = log10(nUMI_Prot), y = mean.1)) + 
  theme_bw() + 
  geom_bin2d(bins = 300, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.3, 0.7), legend.key.size = unit(0.3, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson",  aes(label = ..r.label..))  + 
  geom_smooth(color = "black", method = "lm") + 
  theme(strip.background = element_blank()) + 
  theme(strip.text = element_text(size = 10)) + 
  xlab("log10 prot library size") + 
  ylab(" background mean µ1 ")  + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15))
ggsave(p1, filename = paste0(figpath,"_background_mean_vs_libsize.pdf"), width = 3, height = 4)


# isotype control vs µ1 
p2 = ggplot(data = df, aes(x = iso_mean, y = mean.1))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.25, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson") + 
  geom_smooth(color = "black", method = "lm") + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" mean of isotype controls ") + 
  ylab(" background mean µ1  ")
ggsave(p2, filename = paste0(figpath,"meanisovs_meanbackground.pdf"), width = 3.3, height = 3.3)

# mean 1 mean 2 
# isotype control vs µ1 
p2 = ggplot(data = df, aes(x = mean.1, y = mean.2))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.25, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson") + 
  geom_smooth(color = "black", method = "lm") + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" background mean µ1 ") + 
  ylab(" positive mean µ2  ")
ggsave(p2, filename = paste0(figpath,"mean1_mean2.pdf"), width = 3.3, height = 3.3)


# isotype control vs µ2
p2 = ggplot(data = df, aes(x = iso_mean, y = mean.2))  + 
  theme_bw() + 
  geom_bin2d(bins = 100, show.legend = TRUE) + 
  scale_fill_viridis_c(option = "B") + 
  theme(legend.position = c(0.15, 0.7), legend.key.size = unit(0.25, units = "cm")) + 
  ggpubr::stat_cor( method = "pearson") + 
  geom_smooth(color = "black", method = "lm") + 
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  xlab(" mean of isotype controls ") + 
  ylab(" positive mean µ2  ")
ggsave(p2, filename = paste0(figpath,"meanisovs_mean2.pdf"), width = 3.3, height = 3.3)


### plot result distributions and correlations 
df3 = df %>% gather(protein_population, value, mean.1:mean.2)
p = ggplot(df3, aes(x = value, fill = protein_population)) + 
  geom_histogram(position =  "identity", bins = 100, show.legend = FALSE) + 
  theme_bw() +
  scale_fill_manual(values = c("#3e8ede", "red")) + 
  theme(legend.position = "bottom") + 
  theme(legend.title = element_blank()) + 
  xlab(label = "2 component mixture model means") + 
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 15)) +
  ylab(label = "") +
  xlim(c(-3,15)) 
ggsave(p, filename = paste0(figpath, "FULL_mixture_model_param_dist.pdf"), width = 5, height = 5)


# correlation across latent component.
isotype_value = iso %>% t %>% as.data.frame()
iso_mu = cbind(isotype_value, df$mean.1)
cplot = cor(iso_mu, method = "pearson")
pcor = Hmisc::rcorr(as.matrix(iso_mu), type = "spearman")
colnames(cplot) = rownames(cplot) = c("Isotype 1", "Isotype 2", "Isotype 3", "Isotype 4", "µ1")

# plot correlation matrix 
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
pdf(file = paste0(figpath,"background_correlation_plot_spearman.pdf"), width = 4,height = 4)
corrplot::corrplot(cplot,
                   method="color", 
                   col=col(200),  
         type="upper",
         addCoef.col = "black", 
         cl.pos = "b", 
         tl.col="black", tl.srt=45, 
         cl.lim = c(0,1),
         diag=FALSE, addgrid.col = "ghostwhite"
)
dev.off()


############### additional correlations in text 

# correlation across latent component.
isotype_value = iso %>% t %>% as.data.frame()
isomu = Matrix::rowMeans(isotype_value)
iso_mu = cbind(isotype_value,
               'isotypemean' = isomu, 
               'mu1' = df$mean.1,
               'plibsize' = df$nUMI_Prot)
cplot = cor(iso_mu, method = "pearson")
pdf(file = paste0(figpath,"background_correlation_plot_spearman_extra.pdf"), width = 4,height = 4)
corrplot::corrplot(cplot,
                   method="color", 
                   col=col(200),  
                   type="upper",
                   addCoef.col = "black", 
                   cl.pos = "b", 
                   tl.col="black", tl.srt=45, 
                   cl.lim = c(0,1),
                   diag=F, addgrid.col = "ghostwhite"
)
dev.off()







