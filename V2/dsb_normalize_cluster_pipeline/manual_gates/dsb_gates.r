
GateLinneg = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD3_PROT", ] < 5))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}


# classical mono
GateMono14 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD3_PROT", ] < 5 & 
                         adt["CD14_PROT", ] > 2 & 
                         adt["CD16_PROT", ] < 4 ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

Gateneg16neg56 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD3_PROT", ] < 5 & 
                         adt["CD14_PROT", ] > 2 & 
                         adt["CD14_PROT", ] < 15 & 
                         adt["CD16_PROT", ] < 3  & 
                         adt["CD56_PROT", ] < 5 ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

#### TCELL 
GateLinnegT = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD14_PROT", ] < 3 &
                         adt["CD3_PROT", ] > 5))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

GateT4 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD14_PROT", ] < 3 &
                         adt["CD3_PROT", ] > 5 & 
                         adt["CD4_PROT", ] > 6.5 & 
                         adt["CD8_PROT", ] < 5 
  ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}
GateT8 = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] < 8 &
                         adt["CD14_PROT", ] < 3 &
                         adt["CD3_PROT", ] > 5 & 
                         adt["CD4_PROT", ] < 6.5 & 
                         adt["CD8_PROT", ] > 5 
  ))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}

GateBC = function(SeuratObject, return.seurat = T) {
  adt = as.matrix(SeuratObject@assay$CITE@data)
  cells <- names(which(adt["CD19_PROT", ] > 8 &
                         adt["CD14_PROT", ] < 3 &
                         adt["CD3_PROT", ] < 5))
  if(return.seurat == TRUE) {
    sub = SubsetData(SeuratObject, cells.use = cells)
    return(sub)
  } else { return(cells) }
}