# makes visualizations of dsb, umap results used in Figure 2 and  norm comparison in in Fig 1
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
suppressMessages(library(dsb))
set.seed(1)
'%ni%' = Negate('%in%')

# file paths
figpath = here("V2/dsb_normalize_cluster_pipeline/figures/")
datapath = here("V2/dsb_normalize_cluster_pipeline/generated_data/")

# unstained control dsb norm with empty drops 
un = readRDS(file = here("data/V2_Data/unstained_control_singlets.rds"))
neg_adt = readRDS(here("data/V2_Data/background_data/adt_neg_dmx_list.rds"))
b1cells = un@meta.data[un@meta.data$batch == 1, ]$barcode_check
b2cells = un@meta.data[un@meta.data$batch == 2, ]$barcode_check
pos_adt1 = un@assay$CITE@raw.data[ ,b1cells]
pos_adt2 = un@assay$CITE@raw.data[ ,b2cells]
pos_adt = c(pos_adt1, pos_adt2)
# apply denoised scaled by background protein normalization per batch 
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT", 
             "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")
un_dsb = list()
for (i in 1:length(neg_adt)) {
  un_dsb[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                       empty_drop_matrix = neg_adt[[i]], 
                                       denoise.counts = TRUE, 
                                       use.isotype.control = TRUE, 
                                       isotype.control.name.vec = isotypes)
}
un_dsb = do.call(cbind, un_dsb)
un = SetAssayData(un, assay.type = "CITE",slot = "data",new.data = un_dsb)


###### nondenoised 
un_dsb_nd = list()
for (i in 1:length(neg_adt)) {
  un_dsb_nd[[i]] =  DSBNormalizeProtein(cell_protein_matrix = pos_adt[[i]], 
                                     empty_drop_matrix = neg_adt[[i]], 
                                     denoise.counts = TRUE, 
                                     use.isotype.control = TRUE, 
                                     isotype.control.name.vec = isotypes)
}
un_dsb_nd = do.call(cbind, un_dsb_nd)
un = SetAssayData(un, assay.type = "nond",slot = "data",new.data = un_dsb_nd)



# add visualization 
h1 = readRDS(file = here("V2/dsb_normalize_cluster_pipeline/generated_data/h1_d0_singlets_ann_Seurat2.4_dsbnorm.rds"))
raw = h1@assay$CITE@raw.data
dsb = h1@assay$CITE@data
md = h1@meta.data

# CLR normalize data 
h1 = SetAssayData(h1,assay.type = "CLR", slot = "raw.data",new.data = h1@assay$CITE@raw.data)
h1 = NormalizeData(h1, normalization.method = "genesCLR", assay.type = "CLR")

# repeat for unstained 
un = SetAssayData(un,assay.type = "CLR", slot = "raw.data",new.data = un@assay$CITE@raw.data)
un = NormalizeData(un, normalization.method = "genesCLR", assay.type = "CLR")

# rm outlier cells 
h1 = SubsetData(h1, subset.name = "CD11c_PROT", accept.low = -10)
h1 = SubsetData(h1, subset.name = "CD14_PROT", accept.low = -4, accept.high = 12)


## manual gate distributions of normalization comparison for Figure 1
plot_param = list(
  theme_bw(),
  theme(axis.title.x = element_text(size = 24, face = "bold" ),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18)) 
)


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

## repeat for unstained
log10u = log10(as.matrix(un@assay$CITE@raw.data + 1))
lognormu = NormalizeData(un, assay.type = "CITE",normalization.method = "LogNormalize", scale.factor = 1e4)
arcsinu =  asinh(as.matrix(un@assay$CITE@raw.data))
sqru = sqrt(as.matrix(un@assay$CITE@raw.data))
clru = un@assay$CLR@data


# manual gate plot layers 
mg_layer = list(theme_bw(),  
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                geom_point(size = 0.4, shape = 16, alpha = 0.4) ,  
                geom_vline(xintercept = 0, size = 1.5, color = "red3") ,
                geom_hline(yintercept = 0, size = 1.5, color = "red3") , 
                theme(plot.title = element_text(size = 18) )
)

# arcsin
p1 = ggplot(as.data.frame(t(arcsin)), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +
  ggtitle("Arcsin transformation \nasinh(x)") + xlim(-0.7, 8) + ylim(-0.7, 6) + 
  geom_density_2d(data = as.data.frame(t(arcsinu)), color = '#3e8ede', alpha = 0.9, size = 1)
#geom_point(data = as.data.frame(t(arcsinu)),shape = 21, stroke = 0, fill = '#3e8ede', alpha = 0.5, size = 1.2)

# log norm global library size scaling (as in RNA) for comparison
p2 = ggplot(as.data.frame(t(as.matrix(lognorm@assay$CITE@data))), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer +  ggtitle("library size scaling\nlog(1 + x/sum(x) ) * 1e4") + 
  geom_density_2d(data = as.data.frame(t(as.matrix(lognormu@assay$CITE@data))),color = '#3e8ede', alpha = 0.9, size = 1)

# log 10 
p3 = ggplot(as.data.frame(t(log10)),  aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + ggtitle("Log transformation\nlog10(1 + x)") + xlim(-0.5, 3) + ylim(-0.5, 2.5) + 
  geom_density_2d(data = as.data.frame(t(log10u)),color = '#3e8ede', alpha = 0.9, size = 1)

# CLR 
p4 = ggplot(as.data.frame(t(h1@assay$CLR@data)), aes(x = CD4_PROT, y = CD14_PROT)) + 
  mg_layer + ggtitle("CLR \nNormalization") + xlim(-0.5, 3) + ylim(-0.5, 3) + 
  geom_density_2d(data = as.data.frame(t(un@assay$CLR@data)),color = '#3e8ede', alpha = 0.9, size = 1)

# DSB
p5 = ggplot(as.data.frame(t(h1@assay$CITE@data)), aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer + ggtitle("dsb: both protein-specific ambient noise & \ncell-intrinsic technical component removed") +
  theme(plot.title = element_text(size = 14)) + 
  geom_density_2d(data = as.data.frame(t(un@assay$CITE@data)), color = '#3e8ede', alpha = 0.9, size = 1)
p5
# DSB without denoising
p6 = ggplot(as.data.frame(t(dsb_nodenoise)), aes(x = CD4_PROT, y = CD14_PROT)) +
  mg_layer  + ggtitle("dsb: \nonly protein specific ambient noise removed") + xlim(c(-3, 13)) +  ylim(c(-3,13)) + 
  theme(plot.title = element_text(size = 14)) + 
  geom_density_2d(data = as.data.frame(t(un_dsb_nd)), color = '#3e8ede', alpha = 0.9, size = 1)

p6
# combine 
p = cowplot::plot_grid(p5,p6, p4, p2, p1,p3,nrow = 2)
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained.png"), width = 12.5, height = 8.5)
ggsave(p, filename = paste0(figpath, "norm_comparison_stained_unstained.png"), width = 13, height = 9)
