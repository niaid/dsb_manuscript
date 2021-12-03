suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(here))

# figpath 
datapath = here('V2/table/generated_data/'); dir.create(datapath)

# files with summary data 
tenx_summary_files = list.files(path = here('V2/10x_analysis/generated_data'), pattern = "STATS", full.names = TRUE)
asap_summary_file = list.files(path = here('V2/asapseq/generated_data'),pattern = "STATS", full.names = TRUE)
teaseq_summary_file = list.files(path = here('V2/teaseq/generated_data_v2/TEA-seq'),pattern = "STATS", full.names = TRUE)
h1_summary_file = list.files(path = here('V2/dsb_process_plots/generated_data'),pattern = "STATS", full.names = TRUE)
fil = as.list(c(tenx_summary_files, asap_summary_file, teaseq_summary_file, h1_summary_file))

# save table and format as pdf for fig 3i 
tab = lapply(fil, data.table::fread) %>% bind_rows()
data.table::fwrite(tab,file = paste0(datapath,'combined_table.txt'), sep = "\t")


