# Set the working directory
getwd()
setwd('/home/data/jiaguangshuai/01_brain_atlas')

# Import the Seurat package
library(Seurat)

# Load sample metadata from a CSV file
sample_metadata <- read.csv("./02_sample_info_Singleron.csv")
head(sample_metadata)
tail(sample_metadata)

# Function to create a Seurat object from sample data and attach metadata
create_seurat_from_sample <- function(sample_id, run_id, sample_metadata) {
  data_directory <- file.path("./01_Singleron_matrix/", sample_id, run_id, 'Matrix')
  expression_data <- Read10X(data_directory, gene.column = 2)
  project_identifier <- paste(sample_id, run_id, sep = "_")
  
  seurat_object <- CreateSeuratObject(expression_data, min.cells = 3, min.features = 100, 
                                      project = project_identifier)
  
  # Append metadata to the Seurat object
  metadata_fields <- c('Sample_ID', 'Run_ID', 'Animal', 'Sex', 'Region_1', 'Region_2', 
                       'Region_LR', 'Platform', 'Age')
  for (field in metadata_fields) {
    seurat_object[[field]] <- sample_metadata[[field]]
  }
  seurat_object$sample_id <- sample_id
  
  return(seurat_object)
}

# Process each sample and store the results in a list
processed_samples <- list()
for (i in 1:nrow(sample_metadata)) {
  current_sample_id <- sample_metadata$Sample_ID[i]
  current_run_id <- sample_metadata$Run_ID[i]
  current_metadata <- sample_metadata[i, c('Sample_ID', 'Run_ID', 'Animal', 'Sex', 'Region_1', 'Region_2', 
                                           'Region_LR', 'Platform', 'Age')]
  
  seurat_object <- create_seurat_from_sample(current_sample_id, current_run_id, current_metadata)
  processed_samples[[paste(current_sample_id, current_run_id, sep = "_")]] <- seurat_object
}

# Load lists of mitochondrial and ribosomal genes
mitochondrial_genes <- read.table("./03_mitochondrial_gene_list.txt", header = FALSE, quote = "")$V1
ribosomal_genes <- read.table("./03_ribosome_gene_list.txt", header = FALSE, quote = "")$V1

# Function to calculate and append the percentage of a specific gene set in a Seurat object
append_gene_set_percentage <- function(seurat_object, gene_set, set_label) {
  matched_genes <- intersect(rownames(seurat_object), gene_set)
  seurat_object[[set_label]] <- PercentageFeatureSet(seurat_object, features = matched_genes)
  return(seurat_object)
}

# Append mitochondrial and ribosomal gene percentages to all sample objects
for (sample_name in names(processed_samples)) {
  processed_samples[[sample_name]] <- append_gene_set_percentage(processed_samples[[sample_name]], mitochondrial_genes, "percent.mito")
  processed_samples[[sample_name]] <- append_gene_set_percentage(processed_samples[[sample_name]], ribosomal_genes, "percent.ribo")
}

# Initialize a merged Seurat object for budgerigar brain single-nucleus transcriptomics analysis
Singleron_budgerigar_brain_obj <- processed_samples[[1]]

# Merge all individual Seurat objects into the final combined object
for (i in 2:length(processed_samples)) {
  Singleron_budgerigar_brain_obj <- merge(Singleron_budgerigar_brain_obj, y = processed_samples[[i]])
}

saveRDS(Singleron_budgerigar_brain_obj, file = '03_Singleron_budgerigar_brain_obj.rds')
