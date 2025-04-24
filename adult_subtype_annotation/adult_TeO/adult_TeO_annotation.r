library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")

budgerigar_5months <- readRDS('/data/work/seurat_integrated/09_budgerigar_brain_5months.rds')
budgerigar_5months
Idents(budgerigar_5months) <- budgerigar_5months$`Region_1`
budgerigar_5months_TeO <- subset(budgerigar_5months, idents = c('Midbrain'))
Idents(budgerigar_5months_TeO) <- budgerigar_5months_TeO$`Sample_ID`
# Data preprocessing and dimensionality reduction
budgerigar_brain_obj_list <- SplitObject(budgerigar_5months_TeO, split.by = "Sample_ID")

#NormalizeData
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

budgerigar_5months_TeO <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")
variable_feature <- rownames(budgerigar_5months_TeO)
budgerigar_5months_TeO <- ScaleData(budgerigar_5months_TeO, features = variable_feature, vars.to.regress = c("nCount_RNA"))
options(repr.plot.width = 9, repr.plot.height = 9)
DefaultAssay(budgerigar_5months_TeO) <- "integrated"

# PCA
budgerigar_5months_TeO <- RunPCA(budgerigar_5months_TeO, assay = "integrated", verbose = T) 
budgerigar_5months_TeO <- FindNeighbors(budgerigar_5months_TeO, dims = 1:30, reduction = "pca")
# FindClusters
budgerigar_5months_TeO <- RunUMAP(budgerigar_5months_TeO, dims = 1:30, verbose = T )
for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  budgerigar_5months_TeO <- FindClusters(budgerigar_5months_TeO, resolution = i)
  print(DimPlot(budgerigar_5months_TeO, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(budgerigar_5months_TeO) <- budgerigar_5months_TeO$`integrated_snn_res.0.7`

main_type_anno <- c( 
'0' = 'GABAergic',
'1' = 'Glutamatergic',
'2' = 'Oligodendrocyte',
'3' = 'Glutamatergic',
'4' = 'GABAergic',
'5' = 'GABAergic',
'6' = 'Glutamatergic',
'7' = 'Astrocyte',
'8' = 'Glutamatergic',
'9' = 'GABAergic',
'10' = 'GABAergic',
'11' = 'Glutamatergic',
'12' = 'GABAergic',
'13' = 'GABAergic',
'14' = 'GABAergic',
'15' = 'OPC',
'16' = 'Glutamatergic',
'17' = 'Microglia',
'18' = 'OPC',
'19' = 'GABAergic',
'20' = 'Vascular cell',
'21' = 'Astroependymal',
'22' = 'Blood',
'23' = 'Mural cell')
#names(main_type_anno) 
budgerigar_5months_TeO <- RenameIdents(budgerigar_5months_TeO, main_type_anno)
budgerigar_5months_TeO$TeO_cell_type_res_0.7 <- Idents(budgerigar_5months_TeO)

#Te0_cell_type_res_0.7
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
"Astroependymal" = "#33A1FF"
)
DimPlot(budgerigar_5months_TeO, reduction = "umap", label = T,group.by ="TeO_cell_type_res_0.7",alpha = 0.3,cols = cell_colors)

TeO_Glu <- subset(budgerigar_5months_TeO, idents = c('Glutamatergic'))
TeO_GABA <- subset(budgerigar_5months_TeO, idents = c('GABAergic'))
saveRDS(TeO_Glu, file = '/data/work/5month_TeO/10_TeO_Glu_unannotation.rds')
saveRDS(TeO_GABA, file = '/data/work/5month_TeO/10_TeO_GABA_unannotation.rds')

TeO_Glu <- readRDS('/data/work/5month_TeO/Glu/10_budgerigar_brain_5months_TeO_Glu_annotated_!Ex3_6.rds')
TeO_Glu
TeO_GABA <- readRDS('/data/work/5month_TeO/GABA/10_TeO_GABA_annotated_!Inh3_5_19.rds')
TeO_GABA

##subclass_annotation
budgerigar_5months_TeO$TeO_cell_type_res_0.7 = ifelse(budgerigar_5months_TeO$TeO_cell_type_res_0.7=='Glutamatergic',
                          as.character(TeO_Glu$Glu_level_subclass_annotation[match(colnames(budgerigar_5months_TeO),colnames(TeO_Glu))]),
                          as.character(budgerigar_5months_TeO$TeO_cell_type_res_0.7))
budgerigar_5months_TeO$TeO_cell_type_res_0.7 = ifelse(budgerigar_5months_TeO$TeO_cell_type_res_0.7=='GABAergic',
                          as.character(TeO_GABA$GABA_level_subclass_annotation[match(colnames(budgerigar_5months_TeO),colnames(TeO_GABA))]),
                          as.character(budgerigar_5months_TeO$TeO_cell_type_res_0.7))

# Check for NA values in subclass_annotation column
na_indices <- is.na(budgerigar_5months_TeO$TeO_cell_type_res_0.7)

# Print the number of NA values found
cat("Number of NA values in TeO_subclass_annotation:", sum(na_indices), "\n")

# Remove cells with NA in subclass_annotation
budgerigar_5months_TeO <- subset(budgerigar_5months_TeO, cells = colnames(budgerigar_5months_TeO)[!na_indices])

# Verify that the NAs have been removed
cat("Number of NA values after filtering:", sum(is.na(budgerigar_5months_TeO$TeO_cell_type_res_0.7)), "\n")
budgerigar_5months_TeO$subclass_annotation <- budgerigar_5months_TeO$TeO_cell_type_res_0.7
cell_colors <- c(
"CA8_KIAA1217_Inh" ="#9AABBF",
"CDH9_TLL1_Ex" ="#F099CE",
"CEMIP_RELN_Inh" = "#7B6E9C",   
"EBF1_EBF3_Inh" ="#702783",
"GATA3_DACH1_Inh" ="#D2D4D1" ,
"IRX1_IRX2_Inh" = "#B796B1",
    
"KIAA1217_CACNA1G_Ex"= "#CCCFD4",
"MEIS2_MOXD1_Inh" ="#6AA77E",
"MEIS2_SVIL_Ex" ="#C2661D",
"POSTN_NXPH1_Ex" ="#F4F2A5",
"SCN4B_RGS7_Ex" ="#E4AA60",
"SHOX2_SATB1_Inh" ="#F1D0E5",
"ZFHX4_CHL1_Inh" ="#E7F2E1",      
"ZFHX4_SLIT1_Ex" = "#B2CF9F",   
    
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF",
"Astroependymal" = "#33A1FF"
)

DimPlot(budgerigar_5months_TeO, reduction = "umap", label = F,group.by ="subclass_annotation",alpha =0.9,cols = cell_colors)

saveRDS(budgerigar_5months_TeO, file = '/data/output/11_budgerigar_brain_5months_TeO_annotation_subclass.rds')