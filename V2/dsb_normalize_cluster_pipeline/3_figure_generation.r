# makes visualizations of dsb, umap results used in Figure 2 and  norm comparison in in Fig 1
set.seed(1)
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
'%ni%' = Negate('%in%')

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# read combined dataframe for visualization
df = readRDS(here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_sng_metadata_umapdim_prots_dataframe.rds"))
cu = pals::kelly(n = 22); cu[1] = "dodgerblue" 

# plot umap 
centers = df %>% dplyr::group_by(p3_dist_3) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p = ggplot(df, aes(x = -1*UMAP1, y = UMAP2)) + 
  theme_void() + 
  geom_point(mapping = aes(fill = p3_dist_3), size = 1, shape = 21, stroke = 0, alpha = 0.6, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data = centers, size = 8.5, mapping = aes(label = p3_dist_3), show.legend = FALSE) +
  scale_fill_manual(values = cu) + 
  theme(legend.title =  element_blank()) + 
  theme(legend.text = element_text(colour="black", size=8, face="bold")) 
ggsave(p, filename = paste0(figpath,"h1_umap_DSB.png"),width = 8, height = 7)
ggsave(p, filename = paste0(figpath,"h1_umap_DSB.pdf"),width = 8, height = 7)

# plot distribution of lineage proteins on umap
prot_plot = c("CD3_PROT", "CD4_PROT", "CD8_PROT", "CD56_PROT", "CD16_PROT",
              "CD14_PROT", "CD19_PROT", "IgD_PROT", "CD303_PROT")
index1 = prot_plot[1]
index2 = prot_plot[length(prot_plot)]
df2 = df %>% select(UMAP1, UMAP2, prot_plot) %>% gather(prot, dsb, index1:index2) %>% filter(dsb > -10)
df2$prot = factor(df2$prot, levels = prot_plot)
p = ggplot(df2 %>% filter(prot %in% prot_plot)) + aes(x = -1*UMAP1,y =UMAP2, color = dsb) + 
  geom_point(show.legend = TRUE, size = 0.2, shape = 16, alpha = 0.7) + 
  scale_color_viridis_c(option = "B", limits = c(-5,30)) + 
  theme_void() + 
  theme(strip.background = element_blank()) + 
  theme(legend.key.size  = unit(1, units = "cm")) + 
  theme(strip.text = element_text(color = "black", size = 18, face = "bold"), 
        legend.text = element_text(color = "black", size = 18, face = "bold"),
        legend.title = element_text(color = "black", size = 18, face = "bold")) + 
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
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))
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

# remove prot string  
noise_plot = str_sub(noise_plot, 1, -6)
df2$protein = str_sub(df2$protein, 1, -6)
df1$protein = str_sub(df1$protein, 1, -6)
pplt = df1 %>% filter(median_cluster1_cells > 2.4) %$% protein

# log prot count 
p1 = ggplot(df1, aes( x = median_raw_negative_drop, y = median_cluster1_cells, label = protein, color = plot_noise)) + 
  geom_point(shape = 16 , alpha  = 0.7, show.legend = FALSE) + 
  geom_abline() + 
  theme_bw() + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df1 %>% filter(median_cluster1_cells > 2.4) , size = 4, segment.color = "grey",  show.legend = FALSE) + 
  labs(x = "negative drops median log+1",  y = "highlighted cluster median log+1")
ggsave(p1,filename = paste0(figpath,"nk_log_vs_empty.pdf"), width = 3, height = 3.5)
  
# DSB
p2 = ggplot(df2, aes(x = median_raw_negative_drop, y = median_DSB_cluster1_cells, label = protein, color = plot_noise)) + 
  geom_point(shape = 16 , alpha  = 0.7, show.legend = FALSE) + 
  # geom_point(data = df2 %>% filter(median_DSB_cluster1_cells > 3.5), shape = 21, size = 2.5, fill = '#FF7F0EFF', show.legend = FALSE)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "negative drops median log+1", y = "highlighted cluster median dsb normalized") +
  theme_bw() + 
  geom_hline(yintercept = 3.5, color = "red3") + 
  scale_color_manual(values = c("black", "red")) +
  ggrepel::geom_text_repel(data = df2 %>% filter(protein %in% pplt), segment.color = "grey", size = 4, nudge_y = 0.5, show.legend = FALSE) 
ggsave(p2,filename = paste0(figpath,"nk_dsb_vs_empty.pdf"), width = 3, height = 3.5)


### CLR vs DSB density histograms on nk cluster 
plot_subset =  df1 %>% filter(median_cluster1_cells > 2.4) %>% arrange(desc(median_cluster1_cells)) %$% protein
plot_subset = paste(plot_subset, "_PROT", sep = "")
index1 = plot_subset[1] ; index2 = plot_subset[length(plot_subset)]

# CLR normalize data; all cells
h1 = SetAssayData(h1,assay.type = "CLR", slot = "raw.data",new.data = h1@assay$CITE@raw.data)
h1 = NormalizeData(h1, normalization.method = "genesCLR", assay.type = "CLR")

# subset NK cells
nk = h1 %>% 
  SetAllIdent(id = "celltype_label_3") %>%
  SubsetData(ident.use = "CD16++ NK", subset.raw = TRUE) 

# merge dsb and clr df; redefine plot subset with _prot string
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
  ggridges::geom_density_ridges2(show.legend = FALSE, alpha = 0.8, size = 0.3) + 
  ggsci::scale_fill_d3() + 
  ggsci::scale_color_d3() + 
  theme(strip.background = element_blank()) + 
  theme(axis.text.y =  element_text(size = 8, color = "black")) + 
  xlab("") + ylab("") + 
  geom_vline(xintercept = 0, color = "red2", size = 0.7) 
ggsave(p, filename = paste0(figpath, "dsb_clr_nkcell_distribution_topmedian.pdf"), width = 3.4, height = 4.1)


######## across clusters log vs dsb counts 
dsb_df = cbind(as.data.frame(t(h1@assay$CITE@data)), md) %>% 
  group_by(celltype_label_3, p3_dist_3) %>% 
  gather(prot, count, AnnexinV_PROT:CD20_PROT) %>% 
  group_by(prot, p3_dist_3, celltype_label_3) %>% 
  summarize(median_dsb = median(count), 
  )

# add neg median for each protein  to dsb_df
neg_med = data.frame(median_neg = apply(log1p(raw), 2, median) ) %>% rownames_to_column("prot")
dsb_df = full_join(dsb_df, neg_med, by = "prot")

# cat vars cluster and annotated celltype 
dsb_df$cluster_celltype = paste(dsb_df$p3_dist_3, dsb_df$celltype_label_3, sep = ": ")

# all 
p = ggplot(dsb_df, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ cluster_celltype, nrow = 3) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 10)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab(" median log + 1 counts in empty droplets") + 
  ylab(" dsb normalized median within cluster ") + 
  ggrepel::geom_text_repel(data = dsb_df %>% filter(median_dsb > 4), segment.size = 0.5, size = 2.4, force = 2)
ggsave(p, filename = paste0(figpath,"neg_vs_dsb2.pdf"), width = 18, height = 9)


#### version 2 (subset for main panel)
marker_highlight = c("CD20", "CD19", "IgD", "IgM", "CD123", "CD303", "CD1d", "CD1c", "CD14","CD103", 
                    "CD16", "CD3", "CD4", "CD8", "CD28", "CD161", "CD45RO", "CD45RA", "CD33", "CD86")
highl = c("13: Unswitched B cells", "5: Transitional B cells", "12: Switched B cells",  "16: pDC", "14: mDC",
          "10: Non-classical Monocytes", "2: Classical Monocytes", "3: CD8+ Naive T", "19: CD8+ CD103+ T", "8: Unconventional T CD161+/CD3+")
subs = dsb_df %>% filter(cluster_celltype %in% highl)
subs$cluster_celltype = str_replace(string = subs$cluster_celltype, pattern = "Unconventional T", replacement = "T cell")
highl = str_replace(string = highl, pattern = "Unconventional T", replacement = "T cell")
subs$cluster_celltype = factor(subs$cluster_celltype, levels = highl)
subs$prot = str_sub(subs$prot, 1, -6)
p = ggplot(subs, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point(shape = 16 , alpha  = 0.5) + 
  geom_point(data = subs %>% filter(prot %in% marker_highlight), shape = 21, size = 2.5, fill = '#FF7F0EFF')+
  ggrepel::geom_text_repel(data = subs %>% filter(prot %in% marker_highlight & median_dsb > 3.5), segment.size = 0.5, segment.color = "grey", size = 4.7, force = 3) +
  theme_bw() + 
  geom_segment(data = subs %>% filter(prot %in% marker_highlight ), aes(yend=median_dsb, xend=median_neg)) +
  facet_wrap(~ cluster_celltype, nrow = 2) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 16)) + 
  theme(strip.text = element_text(size = 14)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab(" median log + 1 protein in empty droplets") + 
  ylab(" median dsb normalized  protein ") 
ggsave(p, filename = paste0(figpath,"neg_vs_dsb_highlight.pdf"), width = 14, height = 6)


### remaining clusters in supplement
marker_highlight2 = c("CD56", "CD71", "CD27", "CD244", "KLRG1", "CD195", "CD38", "CD127",  "CD16",
                      "CD3", "CD4", "CD8", "CD28", "CD161", "CD45RO", "CD45RA", "CD34", "CD33", "CD86")
subs2 = dsb_df %>% filter(cluster_celltype %ni% highl)
subs2 = subs2 %>% filter(cluster_celltype %ni% "8: Unconventional T CD161+/CD3+")
subs2$cluster_celltype = str_replace(string = subs2$cluster_celltype, pattern = "Unconventional T", replacement = "T cell")
highl = subs2$cluster_celltype %>% unique 
subs2$prot = str_sub(subs2$prot, 1, -6)
p = ggplot(subs2, aes(x = median_neg, y = median_dsb, label = prot)) + 
  geom_point(shape = 16 , alpha  = 0.5) + 
  geom_point(data = subs2 %>% filter(prot %in% marker_highlight2), shape = 21, size = 2.5, fill = '#FF7F0EFF')+
  ggrepel::geom_text_repel(data = subs2 %>% filter(prot %in% marker_highlight2 & median_dsb > 3.5), segment.size = 0.5, segment.color = "grey", size = 4.7, force = 5) +
  theme_bw() + 
  geom_segment(data = subs2 %>% filter(prot %in% marker_highlight2 ), aes(yend=median_dsb, xend=median_neg)) +
  facet_wrap(~ cluster_celltype, nrow = 2) + 
  theme(strip.background = element_blank(), axis.title = element_text(size = 16)) + 
  theme(strip.text = element_text(size = 10)) + 
  geom_hline(yintercept = 3.5, color = "red") + 
  xlab(" median log + 1 protein in empty droplets") + 
  ylab(" median dsb normalized  protein ") 
ggsave(p, filename = paste0(figpath,"supplneg_vs_dsb_highlight_2.pdf"), width = 16, height = 6)

