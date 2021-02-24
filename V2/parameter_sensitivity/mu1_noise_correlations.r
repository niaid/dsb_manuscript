suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here)) 
suppressMessages(library(Seurat)) 
suppressMessages(library(mclust)) 
set.seed(1)
####### 

figpath = here("V2/parameter_sensitivity/figures_mu1noise/"); dir.create(figpath)
datapath = here("V2/parameter_sensitivity/generateddata_mu1noise/"); dir.create(datapath)

# load dsb step 1 ambient corrected data from b1 
norm_adt1 = readRDS(file = here("V2/dsb_process_plots/generated_data/b1_dsb_norm_adt_mtx.rds"))
h = readRDS(file = "data/V2_Data/CITEseq_raw_PMID32094927_seurat2.4.rds")
isotypes = c("MouseIgG1kappaisotype_PROT", "MouseIgG2akappaisotype_PROT",
             "Mouse IgG2bkIsotype_PROT",  "RatIgG2bkIsotype_PROT")

# downsample uniformly across celltypes 
cell_subset = h %>% 
  SetAllIdent(id = "batch") %>% 
  SubsetData(ident.use = '1') %>% 
  SetAllIdent(id = 'celltype_label_3') %>% 
  SubsetData(random.seed = 1, max.cells.per.ident = 30)
cell_subset = cell_subset@meta.data$barcode_check

##### apply per cell mixture model  
cmd = apply(norm_adt1, 2, function(x){ g = Mclust(x, G=2, warn = TRUE , verbose = TRUE) }) 

# extract vector of all mu1 and mu2 proteins for each cell
m2p = lapply(cmd, function(x) x$classification[x$classification == 2] %>% names())
m1p = lapply(cmd, function(x) x$classification[x$classification == 1] %>% names())
names(m1p) = names(cmd)

# get i = 100 samples of k=4 proteins for each cell
psamples = list()
for (i in 1:length(m1p)) {
  protvec = m1p[[i]]
  protvec = setdiff(protvec, isotypes)
  psamples[[i]] = replicate(100, expr = sample(x = protvec,size = 4, replace = FALSE), simplify = FALSE)
}

# get the mean of the k protein samples from mu 1 for each i for each n  
mx = list()
for (i in 1:length(psamples)) {
  pv = psamples[[i]]
  cell_n = names(cmd)[i]
  mx[[i]] = lapply(pv, function(x){mean(norm_adt1[x,cell_n])})
}

# vector of the resampled mu1 from the k=4 subsamples
mu1_ls = lapply(mx, unlist)
mu1_df = do.call(cbind, mu1_ls) %>% as.matrix() %>% as.data.frame()
colnames(mu1_df) = names(cmd)
saveRDS(mu1_df, file = paste0(datapath, "mu1_df.rds"))

# load cell metadata containing true values mu 2 
md = readRDS(file = "V2/dsb_process_plots/generated_data/mixture_model_metadata_merged.rds")
mdsub = md %>% filter(barcode_check %in% names(cmd)) 


# correlate k=4 mu1 protein subsamples with µ2 
rho2 = list()
yy2 = mdsub$mean.2
for (i in 1:100) {
  xx = mu1_df[i,] %>% as.matrix %>% as.vector
  rho2[[i]] = cor.test(x = xx,  y = yy2,method = "pearson")$estimate
}
rv2 = unlist(rho2, use.names = FALSE) %>% as.vector()
saveRDS(rv2, file = paste0(datapath, "rv2.rds"))

# correlate isotype control mean with full mu 1
yy = Matrix::colMeans(norm_adt1[isotypes, ]) %>% as.matrix %>% as.vector
m2 = mdsub$mean.2
ct2 =  cor.test(m2, yy, method = "pearson")$estimate

# vlsualization and comparison with mu2 and isotype controls
pdf(file = paste0(figpath, "pearson_mu1_cor_m2.pdf"),height = 4, width = 4.5)
rethinking::dens(rv2, 
                 col = '#3e8ede', 
                 adj = 0.5,
                 lwd = 3, 
                 show.HPDI = 0.50, 
                 font.main= 1, cex.main = 0.8,
                 xlim = c(0.24, 0.48), xlab ='Pearson Correlation',
                 main = '100 random samples of k = 4 µ1 proteins from n = 28229 cells: \n Distribution of µ1 k-sample means Pearson correlation w/ µ2 (blue)', 
)
grid(lty = 1, lwd = 1, col = "grey")
abline(v = ct2, col="red", lwd=3, lty=2)
legend('bottom', col = 'red',lty = 2, lwd = 2, cex = 0.7,
       legend = 'Pearson correlation \n isotype controls and µ2')
dev.off()


# show the distribution of mu1 k samples from a single cell  
xx = mu1_df[,1] %>% as.matrix %>% as.vector
true = mdsub$mean.1[1]
pdf(file = paste0(figpath, "singlecell_mu1samplesvstrue.pdf"),height = 4, width = 4.5)
rethinking::dens(xx,
                 font.main= 1, cex.main = 0.8,
                 main = '100 random samples of k = 4 µ1 proteins from a single cell', 
                 lwd = 3,show.HPDI = 0.50)
abline(v = true, col = 'red')
legend('topright', col = 'red',lty = 1, lwd = 2, cex = 0.7,
       legend = 'true µ1 from all proteins \n N1(µ1,sigmasq) in a single cell ')
dev.off()










