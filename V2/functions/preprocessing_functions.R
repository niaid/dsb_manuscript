# Preprocessing Functions 

# These are modified Seurat functions to accomodate Cellrange V3 output within the Seurat V2.3.4 framework
# Changed the names of genes loc to features.tsv to accomodate new naming convention in Cellranger v 3.0 
# Changed the file extensions to .tsv.gz to accomodate cell ranger v 3 compressed output. 

Read10X_MPM = function (data.dir) {
  full.data <- list()
  for (i in seq_along(data.dir)) {
    if (!dir.exists(data.dir)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- paste0(data.dir, "barcodes.tsv.gz")
    gene.loc <- paste0(data.dir, "features.tsv.gz")
    matrix.loc <- paste0(data.dir, "matrix.mtx.gz")
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)) {
      stop("Gene name file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, 
                                                                      FUN = ExtractField, field = 2, delim = "\\t")))
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}




# read10X V3 - direct copy from seurat to accomodate cellrange output. 
Read10X_V3 = function (data.dir = NULL, gene.column = 2) 
{
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(run)) {
      stop("Directory provided does not exist")
    }
    if (!grepl("\\/$", run)) {
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    features.loc <- paste0(run, "features.tsv.gz")
    matrix.loc <- paste0(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing")
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                              yes = gene.loc, no = features.loc), header = FALSE, 
                                stringsAsFactors = FALSE)
    rownames(x = data) <- make.unique(names = feature.names[, 
                                                            gene.column])
    if (ncol(x = feature.names) > 2) {
      if (length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, ])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}
