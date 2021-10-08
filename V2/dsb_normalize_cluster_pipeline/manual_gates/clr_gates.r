# CLR gates for monocyte gating comparison in S3  
#Myeloid  lineage 
GateLinneg = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 2 &
                         adt["CD3_PROT", ] < 0.75))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

# classical monocytes 
GateMono14 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 2 &
                         adt["CD3_PROT", ] < 0.75 & 
                         adt["CD14_PROT", ] > 1 & 
                         adt["CD16_PROT", ] < 2 ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

# subset 
Gateneg16neg56 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 2 &
                         adt["CD3_PROT", ] < 0.75 & 
                         adt["CD14_PROT", ] > 1 & 
                         adt["CD16_PROT", ] < 2  & 
                         adt["CD56_PROT", ] < 1 ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}