library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
budgerigar_5months <- readRDS('/data/work/seurat_integrated/09_budgerigar_brain_5months.rds')
budgerigar_5months
Idents(budgerigar_5months) <- budgerigar_5months$`Region_1`
budgerigar_5months_Cere <- subset(budgerigar_5months, idents = c('Cerebellum'))
Idents(budgerigar_5months_Cere) <- budgerigar_5months_Cere$`Sample_ID`
# Data preprocessing and dimensionality reduction
budgerigar_brain_obj_list <- SplitObject(budgerigar_5months_Cere, split.by = "Sample_ID")

# NormalizeData
for (name in names(budgerigar_brain_obj_list)) {

    budgerigar_brain_obj_list[[name]] <- NormalizeData(budgerigar_brain_obj_list[[name]], normalization.method = "LogNormalize", assay = "RNA")

    budgerigar_brain_obj_list[[name]] <- FindVariableFeatures(budgerigar_brain_obj_list[[name]], assay = "RNA", nfeatures = 2000, selection.method = "vst")
}

integrated_features <- SelectIntegrationFeatures(object.list = budgerigar_brain_obj_list)

budgerigar_brain_obj_list <- lapply(X = budgerigar_brain_obj_list, FUN = function(x) {
    x <- ScaleData(x, features = integrated_features, verbose = FALSE)
    x <- RunPCA(x, features = integrated_features, verbose = FALSE)
})

####rpca
anchors <- FindIntegrationAnchors(object.list = budgerigar_brain_obj_list, anchor.features = integrated_features, reduction = "rpca",k.anchor = 20)

budgerigar_5months_Cere <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")

# Scale data (default: only for variable features)
variable_feature <- rownames(budgerigar_5months_Cere)
DefaultAssay(budgerigar_5months_Cere) <- "integrated"

# PCA
budgerigar_5months_Cere <- RunPCA(budgerigar_5months_Cere, assay = "integrated", verbose = T) 
budgerigar_5months_Cere <- FindNeighbors(budgerigar_5months_Cere, dims = 1:30, reduction = "pca")
# FindClusters
budgerigar_5months_Cere <- RunUMAP(budgerigar_5months_Cere, dims = 1:30, verbose = T )
for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  budgerigar_5months_Cere <- FindClusters(budgerigar_5months_Cere, resolution = i)
  print(DimPlot(budgerigar_5months_Cere, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(budgerigar_5months_Cere) <- budgerigar_5months_Cere$`integrated_snn_res.0.7`

main_type_anno <- c( 
'0' = 'Glutamatergic',
'1' = 'Glutamatergic',
'2' = 'Oligodendrocyte',
'3' = 'Glutamatergic',
'4' = 'Oligodendrocyte',
'5' = 'Bergmann',
'6' = 'Astrocyte',
'7' = 'GABAergic',
'8' = 'Glutamatergic',
'9' = 'GABAergic',
'10' = 'Glutamatergic',
'11' = 'Glutamatergic',
'12' = 'OPC',
'13' = 'Microglia',
'14' = 'Glutamatergic',
'15' = 'Blood',
'16' = 'OPC',
'17' = 'GABAergic',
'18' = 'Vascular cell',
'19' = 'Glutamatergic',
'20' = 'GABAergic',
'21' = 'Glutamatergic',
'22' = 'Mural cell',
'23' = 'GABAergic')
                   
#names(main_type_anno) 
budgerigar_5months_Cere <- RenameIdents(budgerigar_5months_Cere, main_type_anno)
budgerigar_5months_Cere$Cere_cell_type_res_0.7 <- Idents(budgerigar_5months_Cere)

cell_colors <- c(
"Bergmann" ="#c7deef",
'Glutamatergic' = '#E84e40', 
"GABAergic" = "#A4D38E", 
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF",
"Astroependymal" = "#33A1FF"
)
DimPlot(budgerigar_5months_Cere, reduction = "umap", label = F,group.by ="Cere_cell_type_res_0.7",alpha = 0.3,cols = cell_colors)

#load subclass data
Cere_Glu <- readRDS('/data/work/5month_Cere/Glu/10_Cere_Glu_annotated_!Ex1_Ex14.rds')
Cere_GABA <- readRDS('/data/work/5month_Cere/GABA/10_Cere_GABA_annotated_!Inh2.rds')

budgerigar_5months_Cere$Cere_cell_type_res_0.7 = ifelse(budgerigar_5months_Cere$Cere_cell_type_res_0.7=='Glutamatergic',
                          as.character(Cere_Glu$Glu_level_subclass_annotation[match(colnames(budgerigar_5months_Cere),colnames(Cere_Glu))]),
                          as.character(budgerigar_5months_Cere$Cere_cell_type_res_0.7))
budgerigar_5months_Cere$Cere_cell_type_res_0.7 = ifelse(budgerigar_5months_Cere$Cere_cell_type_res_0.7=='GABAergic',
                          as.character(Cere_GABA$GABA_level_subclass_annotation[match(colnames(budgerigar_5months_Cere),colnames(Cere_GABA))]),
                          as.character(budgerigar_5months_Cere$Cere_cell_type_res_0.7))

# Check for NA values in subclass_annotation column
na_indices <- is.na(budgerigar_5months_Cere$Cere_cell_type_res_0.7)

# Print the number of NA values found
cat("Number of NA values in Cere_subclass_annotation:", sum(na_indices), "\n")

# Remove cells with NA in subclass_annotation
budgerigar_5months_Cere <- subset(budgerigar_5months_Cere, cells = colnames(budgerigar_5months_Cere)[!na_indices])

# Verify that the NAs have been removed
cat("Number of NA values after filtering:", sum(is.na(budgerigar_5months_Cere$Cere_cell_type_res_0.7)), "\n")

budgerigar_5months_Cere$subclass_annotation <- budgerigar_5months_Cere$Cere_cell_type_res_0.7

cell_colors <- c(
"Bergmann" ="#c7deef",
"Granule" ="#CE381F", #FF6B95
"Interneuron_1" ="#dc8e97",
"Interneuron_2" ="#e3d1db",
"Purkinje" ="#916ba6",
"SNCG_Ex" ="#d2e0ac",
"TRPC4_UBC" ="#e97371",
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF"
)

DimPlot(budgerigar_5months_Cere, reduction = "umap", label = T,group.by ="subclass_annotation",alpha =0.9,cols = cell_colors)
