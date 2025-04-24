# Set the working directory
setwd('~/01_brain_atlas/')

# Load the necessary R packages
library(Seurat)
library(DoubletFinder)
library(dplyr)

# Load the budgerigar brain data
BGI_budgerigar_brain_obj <- readRDS('03_BGI_budgerigar_brain_obj.rds')
Singleron_budgerigar_brain_obj <- readRDS('03_Singleron_budgerigar_brain_obj.rds')

budgerigar_brain_obj <- merge(BGI_budgerigar_brain_obj, y = Singleron_budgerigar_brain_obj)

# Split the `budgerigar_brain_obj` object by sequencing run
splitSeurat <- SplitObject(budgerigar_brain_obj, split.by = "Run_ID")

# Initialize a list to store processed data
processed_list <- list()

# Loop through each split Seurat object
for (i in 1:length(splitSeurat)) {
  # Pre-process Seurat object with standard methods
  splitSeurat[[i]] <- NormalizeData(splitSeurat[[i]])
  splitSeurat[[i]] <- FindVariableFeatures(splitSeurat[[i]], selection.method = "vst", nfeatures = 2000)
  splitSeurat[[i]] <- ScaleData(splitSeurat[[i]])
  splitSeurat[[i]] <- RunPCA(splitSeurat[[i]])
  splitSeurat[[i]] <- RunUMAP(splitSeurat[[i]], dims = 1:10)
  splitSeurat[[i]] <- FindNeighbors(splitSeurat[[i]], dims = 1:10)
  splitSeurat[[i]] <- FindClusters(splitSeurat[[i]])

  # Perform parametric sweep to find doublets
  sweep_res_list <- paramSweep(splitSeurat[[i]], PCs = 1:10, sct = F)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats) # Get the best parameter point
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() # Extract the best pk value
  
  # Model homotypic cell doublets
  annotations <- splitSeurat[[i]]@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  
  # Calculate expected number of doublets
  n_exp_poi <- round(0.05 * nrow(splitSeurat[[i]]@meta.data))
  n_exp_poi_adj <- round(n_exp_poi * (1 - homotypic_prop))
  
  # Perform doublet detection
  splitSeurat[[i]] <- doubletFinder(splitSeurat[[i]], PCs = 1:10, pN = 0.25, pK = pK_bcmvn, 
                                     nExp = n_exp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

  # Add doublet detection results to metadata
  doub_find_res <- splitSeurat[[i]]@meta.data %>% select(contains('DF.classifications'))
  doub_find_score <- splitSeurat[[i]]@meta.data %>% select(contains('pANN'))
  splitSeurat[[i]] <- AddMetaData(splitSeurat[[i]], metadata = doub_find_res, col.name = "doubFind_res")
  splitSeurat[[i]] <- AddMetaData(splitSeurat[[i]], metadata = doub_find_score, col.name = "doubFind_score")

  # Add the processed object to the list
  processed_list[[i]] <- splitSeurat[[i]]
}


# Merge all processed Seurat objects into a single object
budgerigar_brain_obj <- merge(x = processed_list[[1]], y = processed_list[2:length(processed_list)])

# Save the processed data list
saveRDS(budgerigar_brain_obj, "04_doubletFinder_processed_data.rds")


