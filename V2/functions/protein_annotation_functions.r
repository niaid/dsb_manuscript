# convenience wrapper for avarage heatmap 
AggregateProteinDataPlot = function(protein_assay_data, rownamed_metadata, mdname){
  
  adt = protein_assay_data
  md = rownamed_metadata
  
  adt = adt %>% t %>% as.data.frame 
  proteins = colnames(adt)
  
  adt = adt %>% rownames_to_column("cell")
  adt$mdname = md[[mdname]]
  
  md = md[match(x = rownames(md), table = adt$cell), ]
  
  adt_plot = adt %>% 
    group_by(mdname) %>% 
    summarize_at(.vars = proteins, .funs = base::mean) %>% 
    #return(adt_plot) 
    column_to_rownames("mdname") %>% 
    t %>% 
    as.data.frame
  rownames(adt_plot) = str_sub(rownames(adt_plot), start = 1, end = -6)
  return(adt_plot)
}

AggregateProteinDataMedian = function(protein_assay_data, rownamed_metadata, mdname){
  
  adt = protein_assay_data
  md = rownamed_metadata
  
  adt = adt %>% t %>% as.data.frame 
  proteins = colnames(adt)
  
  adt = adt %>% rownames_to_column("cell")
  adt$mdname = md[[mdname]]
  
  md = md[match(x = rownames(md), table = adt$cell), ]
  
  adt_plot = adt %>% 
    group_by(mdname) %>% 
    summarize_at(.vars = proteins, .funs = stats::median) %>% 
    #return(adt_plot) 
    column_to_rownames("mdname") %>% 
    t %>% 
    as.data.frame
  rownames(adt_plot) = str_sub(rownames(adt_plot), start = 1, end = -6)
  return(adt_plot)
}
