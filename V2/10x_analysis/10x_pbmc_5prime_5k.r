# 10x analysis 
suppressMessages(library(Seurat)) 
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
library(dsb)
library(parallelDist)
library(here)

# set params 
project_title = "10X PBMC 5` 5K V2"
expected_cells = "~5,000"
res = 0.8
kparam = 40

# thresholds 
max_neg_logprotumi = 2.8
min_cell_logprotumi = 3
genemax = 3000
genemin = 80
mtthresh = 0.1


# savepaths 
figpath = paste0(here("V2/10x_analysis/figures/"), project_title, "/")
datapath = here("V2/10x_analysis/generated_data/")
dir.create(figpath); dir.create(datapath)
source("V2/functions/preprocessing_functions.R")
source("V2/functions/analysis_functions.R")

# read data 
raw = readRDS(file = here("data/10x_rds/10x_pbmc_5prime_5k.rds"))

# Separate RNA and protein matrix 
prot = raw[grep(rownames(raw), pattern = "Total"), ]
rna = raw[rownames(raw)[rownames(raw) %ni% rownames(prot)], ]

# calculate metadata 
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE)
pctmt = colSums(rna[mtgene, ])/colSums(rna)
umi  = colSums(rna)
log10umi = log10(umi)
umiprot = colSums(prot)
log10umiprot = log10(umiprot)
nGene = colSums(rna > 0)

# check to see if there are protein detected in drops with no RNA 
pdrop = names(umiprot)[names(umiprot) %ni% names(umi)]
stopifnot(length(pdrop) == 0)

# Confirm RNA assay reads in additional barcodes with no data 
rnadrop = names(umi)[names(umi) %ni% names(umiprot)]
qrna = quantile(umi[rnadrop])

# subset matrices 
int = intersect(names(umiprot), names(umi))
bcmax = length(int)
pctmt = pctmt[int]; umi = umi[int]; log10umi = log10umi[int]; nGene = nGene[int]

# combine into metadata 
md = as.data.frame(cbind(pctmt, umi, log10umi, nGene, umiprot, log10umiprot))

hist(md$log10umiprot[md$log10umiprot < 5], breaks = 100)

# define negative background subset protein matrix to empty drop subset
neg_drops2 = md %>% rownames_to_column("bc") %>% 
  filter(log10umiprot < max_neg_logprotumi & log10umiprot > 2)  %>% 
  filter(nGene < genemin) %$% bc

neg_prot2 = prot[ , neg_drops2] %>%  as.matrix()

# subset out outlier drops from positive protein matrix RNA based; add absolute ceiling for conservative cell estimate. 
positive_cells = md %>% rownames_to_column("bc") %>% 
  filter(log10umiprot > min_cell_logprotumi) %>% 
  filter(nGene < genemax & nGene > 200) %>% 
  filter(pctmt < mtthresh) %$% 
  bc

# subset protein matrix for cells 
pos_prot = prot[ , positive_cells] %>% as.matrix()



###########
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
pv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(pos_prot)) %>% mutate(class = "cell_containing")
nv = md %>% rownames_to_column("bc") %>% filter(bc %in% colnames(neg_prot2)) %>% mutate(class = "background")
ddf = rbind(pv, nv)
p = ggplot(ddf, aes(x = log10umiprot, fill = class, color = class )) +
  ggtitle(paste0(
    project_title, " Threshold 2 \n", 
    "theoretical max barcodes = ", bcmax, "\n", 
    "cell containing drops after QC = ", ncol(pos_prot), "\n",
    "negative droplets = ", ncol(neg_prot2)
  )) + plot_layer
xtop = axis_canvas(p, axis = "x") + geom_density(data = ddf, aes(x = log10umiprot, fill = class)) + ggsci::scale_fill_d3(alpha = 0.5)
p2 = insert_xaxis_grob(p, xtop, grid::unit(.4, "null"), position = "top")
p3 = ggdraw(p2)
ggsave(p3, filename = paste0(figpath,project_title, "protein_joint_lib_distribution.pdf"), width = 4.5, height = 3.5)

# save raw 
saveRDS(pos_prot,file = paste0(datapath, project_title, "pos_prot.rds"))
saveRDS(neg_prot2,file = paste0(datapath, project_title, "neg_prot2.rds"))


#############
# DSB normalize 
# define isotypes 
isotypes = rownames(pos_prot)[13:15]
# threshold 1 normalization
mtx2 = DSBNormalizeProtein(cell_protein_matrix = pos_prot,
                           empty_drop_matrix = neg_prot2,
                           denoise.counts = TRUE,
                           use.isotype.control = TRUE,
                           isotype.control.name.vec = isotypes)

# inspect distributions (plot in next script)
mgl = list(theme_bw(),  geom_point(size  = 0.3, alpha = 0.4) , geom_density_2d(color = 'red', size = 0.3))
dfplot = as.data.frame(t(mtx2))

############# 
# cluster 
rna_cells = rna[ ,positive_cells]
md_cells = md %>% rownames_to_column("bc") %>% filter(bc %in% positive_cells) %>% column_to_rownames("bc")
s = CreateSeuratObject(raw.data = rna_cells, min.cells = 40, min.genes = genemin, meta.data = md_cells)
s = SetAssayData(s, assay.type = "CITE", slot = "data",new.data = mtx2)

##### cluster 
prot = rownames(s@assay$CITE@data)
prot_subset = setdiff(prot, isotypes)

# Subset sd background normalized denoised protein 
s2_adt = GetAssayData(s, assay.type = "CITE", slot = "data")
s2_adt3 = s2_adt[prot_subset, ]
p3_dist = parDist(t(s2_adt3))
p3_dist = as.matrix(p3_dist)

# Protein clustering 
s = FindClusters(s, 
                 distance.matrix = p3_dist,
                 k.param = kparam,
                 print.output = F, 
                 resolution = res,
                 random.seed = 1,
                 algorithm = 3,
                 modularity.fxn = 1)
s = StashIdent(s, save.name = "clusters")

# run umap 
library(reticulate); use_virtualenv("r-reticulate")
library(umap)

# set umap config
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap
ump = umap(t(s2_adt3), config = config)
umap_res = ump$layout %>% as.data.frame() 
colnames(umap_res) = c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb = cbind(s@meta.data, umap_res, as.data.frame(t(s@assay$CITE@data)))
saveRDS(df_dsb, file = paste0(datapath,project_title, "dsb_merged_result.RDS"))

