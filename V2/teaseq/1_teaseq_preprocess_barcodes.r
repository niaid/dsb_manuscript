# This script was provided by Lucas Graybuck

# extract the sameuniquw barcode IDs for the ~30k cells used in the TEAseq manuscript

# " To get these cell_uuids, you would need to pull in the QC metrics from the eLife paper,
# and do some matching: Figure 4—source data 1  Single cell quality metrics for TEA-seq samples.
# I’m attaching a script that will do this join and filtering (swanson_adt_data.R).
# I’m also including the passing cell_uuids (Swanson_filtered_uuids.csv.gz " 


# the script below is a modified version of swanson_adt_data.R that matches directory structure of this proj

# teaeq barcode preprocess 
library(GEOquery)
library(data.table)
library(purrr)
library(dplyr)
suppressMessages(library(here))

elife_dir <- here("data/revision_data/tea_seq_preprocess/Swanson_eLife_data")
dir.create(elife_dir)

geo_dir <-  here("data/revision_data/tea_seq_preprocess/Swanson_GEO_data")
dir.create(geo_dir)

# (MPM) set save path for new data 
datapath = here("data/revision_data/tea_seq_preprocess/generated_data/"); dir.create(datapath)

# Download eLife metrics
elife_qc_url <- "https://cdn.elifesciences.org/articles/63632/elife-63632-fig4-data1-v1.zip"
elife_qc_zip <- file.path(elife_dir, basename(elife_qc_url))
download.file(elife_qc_url, elife_qc_zip)
unzip(elife_qc_zip, exdir = elife_dir)

elife_file <- list.files(elife_dir, pattern = "Figure4.+csv$", full.names = TRUE)

# TEA-seq GEO Sample IDs:
sample_ids <- c("GSM5123951", "GSM5123952", "GSM5123953", "GSM5123954")
# Note: GSM4949911 is also TEA-seq, but from an early run with fewer ADTs.
# These 4 are from the same experiment.

# Get file info from GEO
GSM_info <- lapply(
  sample_ids,
  getGEO,
  GSEMatrix = FALSE
)

# Download cellranger-arc metadata, which include cell_uuids
arc_meta_files <- lapply(
  GSM_info,
  function(GSM) {
    meta <- Meta(GSM)
    supp_files <- meta[grepl("supplementary_file", names(meta))]
    supp_files <- unlist(supp_files)
    supp_files[grepl("arc_per_barcode_metrics", supp_files)]
  }
)

walk(arc_meta_files,
     function(filename) {
       out_file <- basename(filename)
       download.file(filename, file.path(geo_dir, out_file))
     })

# Download ADT counts
adt_files <- lapply(
  GSM_info,
  function(GSM) {
    meta <- Meta(GSM)
    supp_files <- meta[grepl("supplementary_file", names(meta))]
    supp_files <- unlist(supp_files)
    supp_files[grepl("adt_counts", supp_files)]
  }
)

walk(adt_files,
     function(filename) {
       out_file <- basename(filename)
       download.file(filename, file.path(geo_dir, out_file))
     })

### Read and filter metadata

# eLife metrics use sequence-based barcodes suffixed with a well number.
# We can generate these from the barcodes and filenames stored in the arc data for matching

arc_files <- list.files(geo_dir, pattern = "arc_per_barcode_metrics")

arc_meta <- map_dfr(
  arc_files,
  function(arc_file) {
    meta <- fread(file.path(geo_dir, arc_file))
    # Here's where we pull out the well number
    well_num <- sub(".+W([0-9])_.+","\\1",arc_file)
    # and swap this in for the barcode
    meta$barcode <- sub("1$", well_num, meta$barcode)
    meta
  }
)

# Now, let's get the pass/fail calls in from the eLife data:
elife_qc <- read.csv(elife_file)

# filter for passing cells
elife_pass <- elife_qc %>%
  filter(pass_fail == "Pass")

# and filter the metadata to get those UUIDs:
filtered_meta <- arc_meta %>%
  filter(barcode %in% elife_pass$barcode)

nrow(filtered_meta)
# save filtered metadata with new barcode column with well appended 
saveRDS(filtered_meta, file = paste0(datapath, 'filteredmetadata_uniqueIDS_barcode.rds'))
saveRDS(arc_meta, file = paste0(datapath, 'UNfilteredmetadata_uniqueIDS_barcode.rds'))

filtered_uuids <- filtered_meta[ ,"cell_uuid"]
#fwrite(filtered_uuids,file = paste0(datapath,"Swanson_filtered_uuids.csv.gz"))
write_csv(x = data.frame(cell_uuid = filtered_uuids$cell_uuid),
          col_names = TRUE,
          file = paste0(datapath,"Swanson_filtered_uuids.csv"))

### Read and filter ADT data

# We can use the same trick to match up the ADT barcodes
adt_files <- list.files(geo_dir, pattern = "adt_counts")

adt_data <- map_dfr(
  adt_files,
  function(adt_file) {
    adt_dt <- fread(file.path(geo_dir, adt_file))
    # Here's where we pull out the well number
    well_num <- sub(".+W([0-9])_.+","\\1",adt_file)
    # and swap this in for the barcode
    adt_dt$barcode <- paste0(adt_dt$cell_barcode, "-", well_num)
    adt_dt
  }
)
# MPM 
adt_data = adt_data %>% select(barcode,everything())
filtered_adt_data <- adt_data %>%
  filter(barcode %in% elife_pass$barcode)
nrow(filtered_adt_data)

# save filtered and unfiltered adt data with well-appended barcode information 
fwrite(filtered_adt_data,file = paste0(datapath, "Swanson_filtered_adt_counts.csv.gz"))
fwrite(adt_data,file = paste0(datapath, "Swanson_UNfiltered_adt_counts.csv.gz"))
