# makes visualizations of dsb, umap results used in Figure 2 and  norm comparison in in Fig 1
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
'%ni%' = Negate('%in%')

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# read combined dataframe for visualization
df = readRDS("V2/dsb_normalize_cluster_pipeline/generated_data/h1_sng_metadata_umapdim_prots_dataframe.rds")
cu = pals::kelly(n = 22); cu[1] = "dodgerblue" 

# plot umap 
centers = df %>% dplyr::group_by(p3_dist_3) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p = ggplot(df, aes(x = UMAP1, y = UMAP2)) + 
  theme_minimal() + 
  geom_point(mapping = aes(color = p3_dist_3), size = 0.2, alpha = 0.4, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data = centers, size = 5.5, fontface = "bold",
                           mapping = aes(label = p3_dist_3), show.legend = FALSE) +
  scale_color_manual(values = cu) + 
  theme(legend.title =  element_blank()) + 
  theme(legend.text = element_text(colour="black", size=8, face="bold")) 
ggsave(p, filename = paste0(figpath,"h1_umap_DSB.png"),width = 8, height = 7)
ggsave(p, filename = paste0(figpath,"h1_umap_DSB.pdf"),width = 8, height = 7)

# plot distribution of lineage proteins on umap
prot_plot = c("CD3_PROT", "CD4_PROT", "CD8_PROT", "CD56_PROT", "CD16_PROT",
              "CD14_PROT", "CD19_PROT", "IgD_PROT", "CD303_PROT")
index1 = prot_plot[1]
index2 = prot_plot[length(prot_plot)]
df2 = df %>% select(UMAP1, UMAP2, prot_plot) %>% gather(prot, dsb, index1:index2)
df2$prot = factor(df2$prot, levels = prot_plot)
p = ggplot(df2 %>% filter(prot %in% prot_plot)) + aes(x = UMAP1,y =UMAP2, color = dsb) + 
  geom_point(show.legend = TRUE, size = 1, shape = 16, alpha = 0.7) + 
  scale_color_viridis_c(option = "B", limits = c(-5,30)) + 
  theme(strip.background = element_blank()) + 
  theme(legend.key.size  = unit(1, units = "cm")) + 
  theme(strip.text = element_text(color = "black", size = 20, face = "bold"), 
        legend.text = element_text(color = "black", size = 18, face = "bold"),
        legend.title = element_text(color = "black", size = 20, face = "bold")) + 
  facet_wrap(~prot, nrow = 3)
ggsave(p, filename = paste0(figpath,"all_lineage.png"), height = 9, width = 12)
  
## Figure 2 comparison of NK cell noise vs true staining proteins based on DSB in normalized data and neg drops. 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
raw = h1@assay$CITE@raw.data
dsb = h1@assay$CITE@data
md = h1@meta.data

# get NK cell barcode vectors 
dsb_cells = md %>% filter(celltype_label_3 == "CD16++ NK") %$% barcode_check

# subset raw NK data 
raw = as.data.frame(as.matrix(t(raw[ , dsb_cells ])))
dsb = as.data.frame(as.matrix(t(dsb[ , dsb_cells ])))

# load negative drops 
neg_adt = readRDS("data/V2_Data/background_data/adt_neg_dmx_list.rds")
neg_adt = do.call(cbind, neg_adt)
neg = as.data.frame(as.matrix(t(neg_adt)))
noise_plot = c("IgA_PROT", "IgM_PROT", "CD57_PROT","AnnexinV_PROT")

# calc median log + 1 in cells and negative drops add var for nooise proteins to highlight 
df1 = cbind(apply(log1p(raw), 2, median), apply(log1p(neg), 2,median )) %>%
  as.data.frame() %>% 
  rownames_to_column("protein")
names(df1) = c("protein", "median_cluster1_cells", "median_raw_negative_drop")
df1 = df1 %>% mutate(plot_noise = ifelse( protein %in% noise_plot, yes = '1', no = '0'))

# calculate median DSB 
df2 = cbind(apply(dsb, 2, median), apply(log1p(neg), 2,median )) %>% as.data.frame() %>% rownames_to_column("protein")
names(df2) = c("protein", "median_DSB_cluster1_cells", "median_raw_negative_drop") 
prot_plot = df1 %>% filter(median_cluster1_cells > 2.4) %$% protein
df2 = df2 %>% mutate(plot_noise  = ifelse(protein %in% noise_plot, yes = '1', no = '0'))

# log prot count 
p1 = ggplot(df1, aes( x = median_raw_negative_drop, y = median_cluster1_cells, label = protein, color = plot_noise)) + 
  geom_point(show.legend = FALSE) + 
  geom_abline() + 
  theme_bw() + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df1 %>% filter(median_cluster1_cells > 2.4) , size = 2.7, show.legend = FALSE) + 
  labs(x = "negative drops median log+1",  y = "highlighted cluster median log+1")
ggsave(p1,filename = paste0(figpath,"nk_log_vs_empty.pdf"), width = 3, height = 3.7)
  
# DSB
p2 = ggplot(df2, aes(x = median_raw_negative_drop, y = median_DSB_cluster1_cells, label = protein, color = plot_noise)) + 
  geom_point(show.legend = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "negative drops median log+1", y = "highlighted cluster  DSB Normalized") +
  theme_bw() + 
  geom_hline(yintercept = 3, color = "red3") + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df2 %>% filter(protein %in% prot_plot), 
                           size = 2.7, nudge_y = 0.5, show.legend = FALSE) 
ggsave(p2,filename = paste0(figpath,"nk_dsb_vs_empty.pdf"), width = 3, height = 3.7)


### CLR vs DSB density histograms on nk cluster 
plot_subset =  df1 %>% filter(median_cluster1_cells > 2.4) %>% arrange(desc(median_cluster1_cells)) %$% protein
index1 = plot_subset[1] ; index2 = plot_subset[length(plot_subset)]

# CLR normalize data 
h1 = SetAssayData(h1,assay.type = "CLR", slot = "raw.data",new.data = h1@assay$CITE@raw.data)
h1 = NormalizeData(h1, normalization.method = "genesCLR", assay.type = "CLR")

# subset NK cells
nk = h1 %>% 
  SetAllIdent(id = "celltype_label_3") %>%
  SubsetData(ident.use = "CD16++ NK", subset.raw = TRUE) 

# merge dsb and clr df
clrdf = cbind(as.data.frame(t(nk@assay$CLR@data)), nk@meta.data) %>% 
  select(barcode_check, plot_subset) %>% 
  gather(prot, CLR, index1:index2)

dsbdf = cbind(as.data.frame(t(nk@assay$CITE@data)), nk@meta.data) %>%
  select(barcode_check, plot_subset) %>% 
  gather(prot, dsb, index1:index2)

# tidy by norm method. 
mdf = full_join(dsbdf, clrdf) %>% gather(norm_method, value, dsb:CLR)
mdf$prot = factor(mdf$prot, levels = plot_subset)

# visualize distributions
p = ggplot(mdf, aes(x = value, y = prot , fill = norm_method, color = norm_method)) + 
  theme_bw() + 
  facet_wrap(~norm_method, scales = "free_x") + 
  ggridges::geom_density_ridges2(show.legend = FALSE, alpha = 0.3) + 
  ggsci::scale_fill_d3() + 
  ggsci::scale_color_d3() + 
  theme(strip.background = element_blank()) + 
  theme(axis.text.y =  element_text(size = 8, color = "black")) + 
  xlab("") + ylab("") + 
  geom_vline(xintercept = 0, color = "black", size = 0.7) 
ggsave(p, filename = paste0(figpath, "dsb_clr_nkcell_distribution_topmedian.pdf"), width = 3.4, height = 4.1)


######## across clusters log vs dsb counts 
prots = rownames(dsb)
dsb_df = cbind(as.data.frame(t(h1@assay$CITE@data)), md) %>% 
  group_by(celltype_label_3) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  group_by(prot, celltype_label_3) %>% 
  summarize(median_dsb = median(count)
  )
dsb_df$median_neg  = apply(log1p(neg), 2,median )

p = ggplot(dsb_df, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ celltype_label_3, nrow = 3) + 
  theme(strip.background = element_blank(), axis.title = element_text(face = "bold", size = 20)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab(" median log + 1 counts in empty droplets") + 
  ylab(" DSB normalized counts within cell clusters ") + 
  ggrepel::geom_text_repel(data = dsb_df %>% filter(median_dsb > 4), segment.size = 0.5,size = 2.4)
ggsave(p, filename = paste0(figpath,"neg_vs_dsb.pdf"), width = 18, height = 9)


## manual gate distributions of normalization comparison for Figure 1
plot_param = list(
  theme_bw(),
  theme(axis.title.x = element_text(size = 24, face = "bold" ),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18)) 
  )

# rm outlier cells 
h1 = SubsetData(h1, subset.name = "CD11c_PROT", accept.low = -10)
h1 = SubsetData(h1, subset.name = "CD14_PROT", accept.low = -4, accept.high = 12)

# outlier cells removed 
dim(md)[1] - dim(h1@meta.data)[1]
#  58 cells 

# other normalization methods including dsb without denoising  created in script 1
log10 = log10(as.matrix(h1@assay$CITE@raw.data + 1))
lognorm = NormalizeData(h1, assay.type = "CITE",normalization.method = "LogNormalize", scale.factor = 1e4)
arcsin =  asinh(as.matrix(h1@assay$CITE@raw.data))
sqr = sqrt(as.matrix(h1@assay$CITE@raw.data))
clr = h1@assay$CLR@data
dsb_nodenoise = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/NonDenoised_dsb_Mtx.rds"))
dsb_mu1denoise = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/MU1Denoised_dsb_Mtx.rds"))

# manual gate plot layers 
mg_layer = list(theme_bw(),  
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                geom_point(size = 0.2, shape = 16, alpha = 0.5) ,  
                geom_vline(xintercept = 0, size = 0.8, color = "red3") ,
                geom_hline(yintercept = 0, size = 0.8, color = "red3") , 
                theme(plot.title = element_text(size = 12) )
                )

# arcsin
p1 = ggplot(as.data.frame(t(arcsin)), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + ggtitle("Arcsin transformation \nasinh(x)") + xlim(-0.5, 8) + ylim(-0.5, 6)

# log norm global library size scaling (as in RNA) for comparison
p2 = ggplot(as.data.frame(t(as.matrix(lognorm@assay$CITE@data))), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +  ggtitle("library size scaling\nlog(1 + x/sum(x) ) * 1e4") 

# log 10 
p3 = ggplot(as.data.frame(t(log10)),  aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + ggtitle("Log transformation\nlog10(1 + x)") + xlim(-0.5, 3) + ylim(-0.5, 2.5)
  
# CLR 
p4 = ggplot(as.data.frame(t(h1@assay$CLR@data)), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + ggtitle("CLR \nNormalization") + xlim(-0.5, 3) + ylim(-0.5, 3)

# DSB
p5 = ggplot(as.data.frame(t(h1@assay$CITE@data)), aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer + ggtitle("DSB \nstep I and II (default) ") + xlim(c(-3, 13)) +  ylim(c(-3,13))

# DSB without denoising
p6 = ggplot(as.data.frame(t(dsb_nodenoise)), aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer  + ggtitle("DSB \nstep I only") + xlim(c(-3, 13)) +  ylim(c(-3,13))

# combine 
p = cowplot::plot_grid(p5,p6, p4, p2, p1,p3,nrow = 1)
ggsave(p, filename = paste0(figpath, "norm_comparison.png"), width = 15, height = 2.8)
