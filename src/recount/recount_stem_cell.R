## source("https://bioconductor.org/biocLite.R")
## biocLite('recount')
library('recount')

setwd('gpfs/home/jbrugg/map/data/recount/stem_cells/')

#### Find projects of interest ####
## download all sra metadata
metadata <- all_metadata(subset = "sra", verbose = TRUE)
metadata <- subset(metadata, !is.na(sharq_beta_tissue) & !is.na(sharq_beta_cell_type))
metadata <- subset(metadata, 
                   sharq_beta_tissue == "stem cell" | 
                     sharq_beta_cell_type == "esc" | 
                     sharq_beta_cell_type == "ips")

metadata <- subset(metadata, sharq_beta_cell_type != "other")
# project SRP042161 is some weird single cell glioblastoma/gliomasphere experiment
# tissue == "intestine" & cell type == "esc" for all 875 samples

metadata <- subset(metadata, project != "SRP042161")

projects <- unique(metadata$project)

## download projects
project_data <- list()

for(project_id in projects){
  download_study(project_id)
  ## load(file.path(project_id, 'rse_gene.Rdata'))
  ## project_data[project_id] <- rse_gene
}


#### Format project data ####
## format count data
## project_counts_list <- lapply(project_data, function(x) assays(x)$counts)
## project_counts <- do.call(cbind, project_counts_list)
## project_counts <- project_counts[, colnames(project_counts) %in% metadata$run]

## ## format metadata
## # some 54 runs did not have downloads...
## metadata <- subset(metadata, metadata$run %in% colnames(project_counts))
