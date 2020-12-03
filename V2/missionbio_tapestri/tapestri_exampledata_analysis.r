suppressMessages(library(tidyverse))
suppressMessages(library(dsb))
suppressMessages(library(Seurat))
suppressMessages(library(here))

# proj title for fig paths 
project_title = "missionbio tapestri"

## data paths 
figpath = here("V2/missionbio_tapestri/figures/"); dir.create(figpath)
datapath = here("V2/missionbio_tapestri/generated_data/"); dir.create(datapath)
test = read.delim(file = here("data/mission_bio_data/AML-4-cell-line-multiomics-adt-counts.tsv"), sep = "\t",header = T)

# make dataframe
test = as.data.frame(test)
prot = test %>% spread(ab_description, raw) 
prot[is.na(prot)] =  0 

# transpose into cells x prot matrix 
prot = prot %>% 
  column_to_rownames("cell_barcode") %>% 
  t %>% 
  as.data.frame()

# calculate library size of droplets to make rough thresholds for cell containing and ambient droplets 
prot_size = log10(Matrix::colSums(prot, na.rm = TRUE))
md = as.data.frame(prot_size)
md$bc = colnames(prot)
hist(md$prot_size, breaks = 100)

# define a vector of background / empty droplet barcodes based on protein library size
background_drops = md[md$prot_size < 2.5 & md$prot_size > 1.4, ]$bc
negative_mtx_rawprot = prot[ , background_drops] %>%  as.matrix()

# define a vector of cell-containing droplet barcodes based on protein library size 
positive_cells = md[md$prot_size > 2.7, ]$bc
cells_mtx_rawprot = prot[ , positive_cells] %>% as.matrix()

# no isotype data available 
# normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot,
  empty_drop_matrix = negative_mtx_rawprot,
  denoise.counts = FALSE,
  use.isotype.control = FALSE
)


##########################
# visualization of drop distribution 
plot_layer = list(theme_bw() , 
                  ggsci::scale_fill_d3(), ggsci::scale_color_d3() ,
                  geom_histogram(aes(y=..count..), alpha=0.5, bins = 50,position="identity"),
                  geom_density(alpha = 0.5), 
                  ylab("Number of Drops"),  xlab("log10 protein library size"), 
                  theme(axis.title.x = element_text(size = 14)),
                  theme(plot.title = element_text(face = "bold",size = 14)),
                  theme(legend.position = c(0.8, 0.7), legend.margin = margin(0,0,0,0))
)
pv = md  %>% filter(bc %in% colnames(cells_mtx_rawprot)) %>% mutate(class = "cell_containing")
nv = md  %>% filter(bc %in% colnames(negative_mtx_rawprot)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = prot_size, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", nrow(test), "\n", 
    "estimated cell containing drops  = ", ncol(cells_mtx_rawprot), "\n",
    "negative droplets = ", ncol(negative_mtx_rawprot)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = prot_size, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)


##########################
# clustering analysis 

# calculate distance matrix for clustering 
p_dist = dist(t(dsb_norm_prot))
p_dist = as.matrix(p_dist)

# Graph based clustering 
s = CreateSeuratObject(raw.data =  cells_mtx_rawprot, min.cells = 0, min.genes = 0)
s = FindClusters(s, resolution = 0.5, distance.matrix = p_dist)

# add CLR norm 
s = SetAssayData(object = s, assay.type = "CLR", slot = "raw.data",new.data = cells_mtx_rawprot)
s = NormalizeData(s, assay.type = "CLR", normalization.method = "genesCLR")

# heatmap of protein values 
prots = rownames(s@raw.data)
adt_data = cbind(s@meta.data, as.data.frame(t(dsb_norm_prot)))
adt_plot = adt_data %>% 
  group_by(res.0.5) %>% 
  summarize_at(.vars = prots, .funs = mean) %>% 
  column_to_rownames("res.0.5") %>% 
  t %>% 
  as.data.frame

x = pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"),
                   filename = paste0(figpath,project_title, "_tapestri_heatmap.pdf"),
                   fontsize_row = 8, border_color = NA, width = 4, height = 4)


# umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(dsb_norm_prot), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(dsb_norm_prot)))
df_clr = cbind(s@meta.data, umap_res, as.data.frame(t(s@assay$CLR@data)))
saveRDS(df_dsb, file = paste0(datapath, project_title,"df_dsb.rds"))
saveRDS(df_clr, file = paste0(datapath, project_title,"df_clr.rds"))

#######################
### Visualizations 

# clusters by umap dims 
p = ggplot(df_dsb, aes(x = UMAP_1, y= UMAP_2, color = res.0.5 )) + theme_bw() + 
  geom_point(shape = 16) + 
  ggsci::scale_color_d3(alpha = 0.7)
ggsave(p, filename = paste0(figpath,"tapestri_umap_dsb.pdf"), width = 4, height = 4)

# long format for DSB heatmap 
index1 = colnames(df_dsb)[7]; index2 = colnames(df_dsb)[ncol(df_dsb)]
dsb = df_dsb %>% gather(prot, DSB, index1:index2)

# plots 
p = ggplot(dsb, aes(x = UMAP_1, y = UMAP_2, color = DSB)) +
  geom_point(shape = 16) + 
  theme(legend.key.size = unit(0.8, units = "cm"), 
        legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 20, face = "bold")) + 
  scale_color_viridis_c(option = "B") +
  facet_wrap(~ prot,nrow = 2) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15, face = "bold")) + 
  theme(legend.position = "right")
ggsave(p, filename = paste0(figpath,  "clusters_DSBdist.png"), width = 12, height = 7.4)

# scatter plot flow 
mg = list(
  geom_point(shape = 16),
  geom_density_2d(), 
  geom_vline(xintercept = 0, color = "red3"),
  geom_hline(yintercept = 0, color = "red3") 
  )

p1 = ggplot(df_dsb, aes(x = CD19_0D, CD3_5A)) + mg
p2 = ggplot(df_dsb, aes(x = CD34_9C, CD30_6E)) + mg
p = plot_grid(p1, p2)
ggsave(p, filename = paste0(figpath, "scatter_prots_tapestri.pdf"),width = 6, height = 3 )


# vln plot 
dsb$prot = factor(dsb$prot, levels = x$tree_row$labels[x$tree_row$order])
dsb$res.0.5 = factor(dsb$res.0.5, levels = x$tree_col$labels[x$tree_col$order])
p = ggplot(dsb, aes(x = prot, y = DSB, fill = prot, color = prot)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 10)) + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) +
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  coord_flip() 
p
ggsave(p, filename = paste0(figpath, "vln_prots_tapestri.pdf"),width = 10, height = 4)


#### add CLR plot 
# save results dataframe 

# long format for CLR heatmap 
index1 = colnames(df_dsb)[7]; index2 = colnames(df_dsb)[ncol(df_dsb)]
clr = df_clr %>% gather(prot, CLR, index1:index2)
clr$prot = factor(dsb$prot, levels = x$tree_row$labels[x$tree_row$order])
clr$res.0.5 = factor(dsb$res.0.5, levels = x$tree_col$labels[x$tree_col$order])
p = ggplot(clr, aes(x = prot, y = CLR, fill = prot, color = prot)) + 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_text(size = 10)) + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) +
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1, scales = "free_x") +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  coord_flip() 
ggsave(p, filename = paste0(figpath, "CLR_vln_prots_tapestri.pdf"),width = 10, height = 4)
p = ggplot(clr, aes(x = prot, y = CLR, fill = prot)) + 
  theme_bw() +
  ggsci::scale_fill_d3() + 
  geom_violin(scale = "count", show.legend = F, draw_quantiles = TRUE) + 
  geom_hline(yintercept = 0, color = "red3") + 
  facet_wrap(~res.0.5, nrow = 1) +
  coord_flip()

ggsave(p, filename = paste0(figpath, "CLRvln_prots_tapestri.pdf"),width = 6, height = 3 )

