library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
budgerigar_brain_obj <- readRDS('/data/work/seurat_integrated/06_budgerigar_brain_obj_rpca.rds')
options(repr.plot.width = 9, repr.plot.height = 9)
DefaultAssay(budgerigar_brain_obj) <- "integrated"
# PCA
budgerigar_brain_obj <- RunPCA(budgerigar_brain_obj, assay = "integrated", verbose = T) 
budgerigar_brain_obj <- FindNeighbors(budgerigar_brain_obj, dims = 1:30, reduction = "pca")
budgerigar_brain_obj <- FindClusters(budgerigar_brain_obj, resolution = 0.4)
# UMAP, TSNE, FindNeighbors  
budgerigar_brain_obj <- RunUMAP(budgerigar_brain_obj, dims = 1:30, verbose = T )
saveRDS(budgerigar_brain_obj, file = '/data/work/seurat_integrated/07_budgerigar_brain_obj_unannotated.rds')