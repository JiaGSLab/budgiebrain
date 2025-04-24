library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
budgerigar_5months <- readRDS('/data/work/seurat_integrated/09_budgerigar_brain_5months.rds')
budgerigar_5months
Idents(budgerigar_5months) <- budgerigar_5months$`Region_1`
budgerigar_5months_Th <- subset(budgerigar_5months, idents = c('Thalamus'))
Idents(budgerigar_5months_Th) <- budgerigar_5months_Th$`Sample_ID`

# Data preprocessing and dimensionality reduction
budgerigar_brain_obj_list <- SplitObject(budgerigar_5months_Th, split.by = "Sample_ID")

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
anchors <- FindIntegrationAnchors(object.list = budgerigar_brain_obj_list, anchor.features = integrated_features, reduction = "rpca")

budgerigar_5months_Th <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")

variable_feature <- rownames(budgerigar_5months_Th)
budgerigar_5months_Th <- ScaleData(budgerigar_5months_Th, features = variable_feature, vars.to.regress = c("nCount_RNA"))

options(repr.plot.width = 9, repr.plot.height = 9)
DefaultAssay(budgerigar_5months_Th) <- "integrated"

# PCA
budgerigar_5months_Th <- RunPCA(budgerigar_5months_Th, assay = "integrated", verbose = T) 
budgerigar_5months_Th <- FindNeighbors(budgerigar_5months_Th, dims = 1:40, reduction = "pca")
# FindClusters
budgerigar_5months_Th <- RunUMAP(budgerigar_5months_Th, dims = 1:40, verbose = T )
for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  budgerigar_5months_Th <- FindClusters(budgerigar_5months_Th, resolution = i)
  print(DimPlot(budgerigar_5months_Th, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(budgerigar_5months_Th) <- budgerigar_5months_Th$`integrated_snn_res.0.7`

main_type_anno <- c( 
'0' = 'Glutamatergic',
'1' = 'Glutamatergic',
'2' = 'Oligodendrocyte',
'3' = 'Glutamatergic',
'4' = 'Astrocyte',
'5' = 'Oligodendrocyte',
'6' = 'Glutamatergic',
'7' = 'Glutamatergic',
'8' = 'GABAergic',
'9' = 'GABAergic',
'10' = 'Oligodendrocyte',
'11' = 'OPC',
'12' = 'OPC',
'13' = 'Glutamatergic',
'14' = 'Microglia',
'15' = 'Glutamatergic',
'16' = 'Glutamatergic',
'17' = 'GABAergic',
'18' = 'Dopaminergic',
'19' = 'Mural cell',
'20' = 'Blood',
'21' = 'GABAergic',
'22' = 'Glutamatergic',
'23' = 'Astroependymal',
'24' = 'Vascular cell'    
)
                   

#names(main_type_anno) 
budgerigar_5months_Th <- RenameIdents(budgerigar_5months_Th, main_type_anno)
budgerigar_5months_Th$Th_cell_type_res_0.7 <- Idents(budgerigar_5months_Th)
cell_colors <- c(
'Glutamatergic' = '#E84e40', 
"GABAergic" = "#A4D38E", 
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF",
"Astroependymal" = "#33A1FF",
"Dopaminergic" = "#F29B8D"
)
DimPlot(budgerigar_5months_Th, reduction = "umap", label = F ,group.by ='Th_cell_type_res_0.7',cols = cell_colors,alpha = 0.3)



Th_Glu <- readRDS('/data/work/5month_Th/Glu/10_budgerigar_brain_5months_Th_Glu_annotated!GABA_Ex3.rds')
Th_Glu
Th_GABA <- readRDS('/data/work/5month_Th/GABA/10_budgerigar_brain_5months_Th_GABA_annotated!Glu.rds')
Th_GABA

budgerigar_5months_Th$Th_cell_type_res_0.7 = ifelse(budgerigar_5months_Th$Th_cell_type_res_0.7=='Glutamatergic',
                          as.character(Th_Glu$Glu_level_subclass_annotation[match(colnames(budgerigar_5months_Th),colnames(Th_Glu))]),
                          as.character(budgerigar_5months_Th$Th_cell_type_res_0.7))
budgerigar_5months_Th$Th_cell_type_res_0.7 = ifelse(budgerigar_5months_Th$Th_cell_type_res_0.7=='GABAergic',
                          as.character(Th_GABA$GABA_level_subclass_annotation[match(colnames(budgerigar_5months_Th),colnames(Th_GABA))]),
                          as.character(budgerigar_5months_Th$Th_cell_type_res_0.7))

# Check for NA values in subclass_annotation column
na_indices <- is.na(budgerigar_5months_Th$Th_cell_type_res_0.7)

# Print the number of NA values found
cat("Number of NA values in Th_subclass_annotation:", sum(na_indices), "\n")

# Remove cells with NA in subclass_annotation
budgerigar_5months_Th <- subset(budgerigar_5months_Th, cells = colnames(budgerigar_5months_Th)[!na_indices])

# Verify that the NAs have been removed
cat("Number of NA values after filtering:", sum(is.na(budgerigar_5months_Th$Th_cell_type_res_0.7)), "\n")
budgerigar_5months_Th$subclass_annotation <- budgerigar_5months_Th$Th_cell_type_res_0.7
cell_colors <- c(
"ISL1_MEIS2_Inh" ="#AF7FD3",
"KIAA1217_ITPR1_Ex" ="#909FC1",
"CEMIP_MEIS2_Inh" = "#D697C4",   
"CCK_FOXP2_Ex" ="#A2A2A2",
"DLX1_MEIS2_Inh" ="#98B3A2" ,
"Dopaminergic" = "#C8CEDA",
"ENSMUNG00000008668_Inh"= "#A9C4E1",
"GRIA3_GRIK1_Ex" ="#5E83A0",
"NHP2_HMOX2_Inh" ="#E4D3DD",
"SLIT1_PTPRF_Ex" ="#FE6D96",
"SNCG_GATA3_Inh" ="#D8BC97",
"SNCG_NEFL_Ex" ="#D37639",  
"VIP_NTS_Ex" = "#81AF9D",
    
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF",
"Astroependymal" = "#33A1FF"
)
DimPlot(budgerigar_5months_Th, reduction = "umap", label = T ,group.by ='subclass_annotation',alpha =0.9,cols = cell_colors)