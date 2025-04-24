library(Seurat)
budgerigar_brain_obj <- readRDS('/data/input/Files/scRNA-seq/04_doubletFinder_processed_data.rds')
# Calculate log10GenesPerUMI
budgerigar_brain_obj$log10GenesPerUMI <- log10(budgerigar_brain_obj$nFeature_RNA) / log10(budgerigar_brain_obj$nCount_RNA)

# Subset samples based on criteria
budgerigar_brain_obj <- subset(budgerigar_brain_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & 
                                 percent.mito < 5 &
                                 nCount_RNA > 300 & nCount_RNA < 10000 & log10GenesPerUMI > 0.8 &
                                 doubFind_res == 'Singlet')
# Data preprocessing and dimensionality reduction
budgerigar_brain_obj_list <- SplitObject(budgerigar_brain_obj, split.by = "Sample_ID")

# NormalizeDataFindVariableFeatures
for (name in names(budgerigar_brain_obj_list)) {
   
    budgerigar_brain_obj_list[[name]] <- NormalizeData(budgerigar_brain_obj_list[[name]], normalization.method = "LogNormalize", assay = "RNA")

    
    budgerigar_brain_obj_list[[name]] <- FindVariableFeatures(budgerigar_brain_obj_list[[name]], assay = "RNA", nfeatures = 2000, selection.method = "vst")
}