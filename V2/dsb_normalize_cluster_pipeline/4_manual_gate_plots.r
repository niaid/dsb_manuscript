##############
# manual gate plots for Fig 2 and S3
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
'%ni%' = Negate('%in%')

# file paths

figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# data  
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
h1 = SetAssayData(h1,assay.type = "CLR", slot = "raw.data",new.data = h1@assay$CITE@raw.data)
h1 = NormalizeData(h1, normalization.method = "genesCLR", assay.type = "CLR")
# rm ~ 50 outlier cells for visualization 
h1 = SubsetData(h1, subset.name = "CD11c_PROT", accept.low = -10)
h1 = SubsetData(h1, subset.name = "CD14_PROT", accept.low = -4, accept.high = 12)

# lineage 
df = h1@assay$CITE@data %>% t %>% as.data.frame()
p = ggplot(df, aes(x = CD19_PROT, y = CD3_PROT )) + 
  plot_param +
  geom_bin2d(bins = 300)+
  scale_fill_viridis_c(option = "B") + 
  geom_density_2d(color = "white")+
  geom_hline(yintercept = 5) +
  geom_vline(xintercept =8)
ggsave(p, filename = paste0(figpath, "dsb_lineage.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath, "dsb_lineage.pdf"), width = 6, height = 5)

# source DSB gate functions 
source(here("V2/dsb_normalize_cluster_pipeline/manual_gates/dsb_gates.r"))

# lineage negative and outlier cell removal for plots
linneg = GateLinneg(h1, return.seurat = T)
linneg = SubsetData(linneg, subset.name = "CD14_PROT",accept.high = 14)

# plot 
p = GenePlot4(linneg, gene1 = "CD14_PROT", gene2 = "CD16_PROT", pt.size = 0.5) +
  plot_param + 
  geom_hline(yintercept = 4) + geom_vline(xintercept = 2)
ggsave(p, filename = paste0(figpath, "dsb_linneg_cd14_cd16.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath, "dsb_linneg_cd14_cd16.pdf"), width = 6, height = 5)

# Classical monocytes 
cd14 = GateMono14(linneg, return.seurat = T)
subset1 = SubsetData(linneg, cells.use = linneg@cell.names[!linneg@cell.names %in% cd14@cell.names])

# plot 
p = GenePlot4(subset1, gene1 = "CD56_PROT", gene2 = "CD16_PROT", pt.size = 0.5) + plot_param 
ggsave(p, filename = paste0(figpath,"dsb_subset1_cd56_cd16.png"),  width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_subset1_cd56_cd16.pdf"),  width = 6, height = 5)

# dump gate 
neg16neg56 = Gateneg16neg56(subset1,return.seurat = T)
subset2 = SubsetData(subset1, cells.use = subset1@cell.names[!(subset1@cell.names %in% neg16neg56@cell.names)])

# non classical and intermediate
p = GenePlot4(subset2, gene1 = "CD11c_PROT", gene2 = "CD14_PROT", pt.size = 0.5) +
  plot_param +
  geom_hline(yintercept = 2.5) + geom_vline(xintercept = 9)
ggsave(p, filename = paste0(figpath,"dsb_subset2_cd11c_cd14.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_subset2_cd11c_cd14.pdf"), width = 6, height = 5)


#### T cell gating 
linneg = GateLinnegT(h1, return.seurat = T)
df = linneg@assay$CITE@data %>% t %>% as.data.frame()
p = ggplot(df, aes(x = CD4_PROT, y = CD8_PROT )) + 
  plot_param +
  geom_bin2d(bins = 300)+
  scale_fill_viridis_c(option = "B") + 
  geom_density_2d(color = "white")+
  geom_hline(yintercept = 5) +
  geom_vline(xintercept =6.5)
ggsave(p, filename = paste0(figpath, "dsb_Tcells.pdf"), width = 6, height = 5)

# cd4 and CD8 T cells 
cd4 = GateT4(linneg, return.seurat = T)
cd8 = GateT8(linneg, return.seurat = T)
p =GenePlot4(cd4, gene1 = "CD45RO_PROT", gene2 = "CD62L_PROT", pt.size = 0.5) +
  plot_param +
  geom_hline(yintercept = 6) + geom_vline(xintercept = 6)
ggsave(p, filename = paste0(figpath,"dsb_cd4.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_cd4.pdf"), width = 6, height = 5)

p =GenePlot4(cd8, gene1 = "CD45RO_PROT", gene2 = "CD62L_PROT", pt.size = 0.5) +
  plot_param +
  geom_hline(yintercept = 6) + geom_vline(xintercept = 6)
ggsave(p, filename = paste0(figpath,"dsb_cd8.png"), width = 6, height = 5)
ggsave(p, filename = paste0(figpath,"dsb_cd8.pdf"), width = 6, height = 5)

# B cell 
bc = GateBC(h1, return.seurat = T)
p =GenePlot4(bc, gene1 = "IgD_PROT", gene2 = "CD27_PROT", pt.size = 0.5) +
  plot_param +
  geom_hline(yintercept = 4) + geom_vline(xintercept = 10)
ggsave(p, filename = paste0(figpath,"Bcell.pdf"), width = 6, height = 5)


## plot same distributions with CLR 
source(here("V2/dsb_normalize_cluster_pipeline/manual_gates/clr_gates.r"))

# define separate object  
cl = SubsetData(h1,cells.use = h1@cell.names)
cl@assay$CITE = cl@assay$CLR
cl@assay$CLR = NULL
GenePlot4(cl,gene1 = "CD3_PROT", gene2 = "CD19_PROT")

linneg = GateLinneg(cl, return.seurat = T)
# p = GenePlot4(h1, gene1 = "CD19_PROT", gene2 = "CD3_PROT", x.axis.size = 18, y.axis.size = 18, pt.size = 0.1)
# ggsave(p, filename = "figures/cd19cells.png",width = 4, height = 3)
p = GenePlot4(linneg, gene1 = "CD14_PROT", gene2 = "CD16_PROT",
              x.axis.size = 30, y.axis.size = 30, pt.size = 0.1) + theme_classic()
ggsave(p, filename = paste0(figpath, "clr_linneg_cd14_cd16.png"),width = 7, height = 5)
ggsave(p, filename = paste0(figpath, "clr_linneg_cd14_cd16.pdf"),width = 7, height = 5)


# classical mono
cd14 = GateMono14(linneg, return.seurat = T)
subset1 = SubsetData(linneg, cells.use = linneg@cell.names[!linneg@cell.names %in% cd14@cell.names])
p = GenePlot4(subset1, gene1 = "CD56_PROT", gene2 = "CD16_PROT", pt.size = 0.1) + theme_classic()
ggsave(p, filename = paste0(figpath, "clr_subset1_cd56_cd16.png"),width = 7, height = 5)
ggsave(p, filename = paste0(figpath, "clr_subset1_cd56_cd16.pdf"),width = 7, height = 5)

# non classical and intermediate 
neg16neg56 = Gateneg16neg56(subset1,return.seurat = T)
subset2 = SubsetData(subset1, cells.use = subset1@cell.names[!(subset1@cell.names %in% neg16neg56@cell.names)])
p =GenePlot4(subset2, gene1 = "CD11c_PROT", gene2 = "CD14_PROT", pt.size = 0.1) + theme_classic()
p
ggsave(p, filename = paste0(figpath,"clr_subset2_cd11c_cd14.png"),width = 7, height = 5)
ggsave(p, filename = paste0(figpath,"clr_subset2_cd11c_cd14.pdf"),width = 7, height = 5)
