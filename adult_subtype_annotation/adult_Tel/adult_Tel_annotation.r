library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
budgerigar_brain_obj <- readRDS('/data/work/seurat_integrated/08_budgerigar_brain_obj_annotated.rds')
Idents(budgerigar_brain_obj) <- budgerigar_brain_obj$`Age`
budgerigar_5months <- subset(budgerigar_brain_obj, idents = c('5months'))
saveRDS(budgerigar_5months, file = '/data/work/seurat_integrated/09_budgerigar_brain_5months.rds')
Idents(budgerigar_5months) <- budgerigar_5months$`Region_1`
budgerigar_5months_Tel <- subset(budgerigar_5months, idents = c('Telencephalon'))
#5month_Tel
Idents(budgerigar_5months_Tel) <- budgerigar_5months_Tel$`Sample_ID`
# Data preprocessing and dimensionality reduction
budgerigar_brain_obj_list <- SplitObject(budgerigar_5months_Tel, split.by = "Sample_ID")

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

budgerigar_5months_Tel <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")

variable_feature <- rownames(budgerigar_5months_Tel)
budgerigar_5months_Tel <- ScaleData(budgerigar_5months_Tel, features = variable_feature, vars.to.regress = c("nCount_RNA"))
options(repr.plot.width = 9, repr.plot.height = 9)
DefaultAssay(budgerigar_5months_Tel) <- "integrated"

# PCA

budgerigar_5months_Tel <- RunPCA(budgerigar_5months_Tel, assay = "integrated", verbose = T) 
budgerigar_5months_Tel <- FindNeighbors(budgerigar_5months_Tel, dims = 1:30, reduction = "pca")
# Run UMAP for visualization
budgerigar_5months_Tel <- RunUMAP(budgerigar_5months_Tel, dims = 1:30)

# FindClusters
budgerigar_5months_Tel <- RunUMAP(budgerigar_5months_Tel, dims = 1:30, verbose = T )
for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  budgerigar_5months_Tel <- FindClusters(budgerigar_5months_Tel, resolution = i)
  print(DimPlot(budgerigar_5months_Tel, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(budgerigar_5months_Tel) <- budgerigar_5months_Tel$`integrated_snn_res.1.2`

main_type_anno <- c( '0' = 'Glutamatergic',
'1' = 'GABAergic',
'2' = 'Glutamatergic',
'3' = 'Glutamatergic',
'4' = 'Astrocyte',
'5' = 'Glutamatergic',
'6' = 'GABAergic',
'7' = 'GABAergic',
'8' = 'Glutamatergic',
'9' = 'Glutamatergic',
'10' = 'GABAergic',
'11' = 'GABAergic',
'12' = 'Glutamatergic',
'13' = 'Oligodendrocyte',
'14' = 'Glutamatergic',
'15' = 'Glutamatergic',
'16' = 'Glutamatergic',
'17' = 'Glutamatergic',
'18' = 'Glutamatergic',
'19' = 'Glutamatergic',
'20' = 'GABAergic',
'21' = 'Glutamatergic',
'22' = 'OPC',
'23' = 'Glutamatergic',
'24' = 'GABAergic',
'25' = 'GABAergic',
'26' = 'Astroependymal',
'27' = 'Mural cell',
'28' = 'Microglia',
'29' = 'Oligodendrocyte',
'30' = 'Blood',
'31' = 'OPC',
'32' = 'GABAergic',
'33' = 'Astrocyte',
'34' = 'GABAergic',
'35' = 'Astrocyte')
                   

#names(main_type_anno) 
budgerigar_5months_Tel <- RenameIdents(budgerigar_5months_Tel, main_type_anno)
budgerigar_5months_Tel$Tel_cell_type_res_1.2 <- Idents(budgerigar_5months_Tel)

#re-anno_vasc
Idents(mural) <- mural$`main_cell_type_res_0.4`
mural <- subset(mural, idents = c('Mural cell','Vascular cell'))
budgerigar_5months_Tel$Tel_cell_type_res_1.2 = ifelse(budgerigar_5months_Tel$Tel_cell_type_res_1.2=='Mural cell',
                          as.character(mural$main_cell_type_res_0.4[match(colnames(budgerigar_5months_Tel),colnames(mural))]),
                          as.character(budgerigar_5months_Tel$Tel_cell_type_res_1.2))

#Tel_cell_type_res_1.2
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
"Purkinje" = "#F3BAA5",
"Astroependymal" = "#33A1FF"
)
DimPlot(budgerigar_5months_Tel, label = FALSE, reduction = "umap",raster=FALSE,group.by ="Tel_cell_type_res_1.2", alpha = 0.3,cols = cell_colors)


##loading GABA GLU re-clustered data
Tel_Glu <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_Glu_annotated_!Ex5.rds')
Tel_GABA <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_GABA_annotated.rds')

##annotaion_subclass
budgerigar_5months_Tel$Tel_cell_type_res_1.2 = ifelse(budgerigar_5months_Tel$Tel_cell_type_res_1.2=='Glutamatergic',
                          as.character(Tel_Glu$Glu_level_subclass_annotation[match(colnames(budgerigar_5months_Tel),colnames(Tel_Glu))]),
                          as.character(budgerigar_5months_Tel$Tel_cell_type_res_1.2))
budgerigar_5months_Tel$Tel_cell_type_res_1.2 = ifelse(budgerigar_5months_Tel$Tel_cell_type_res_1.2=='GABAergic',
                          as.character(Tel_GABA$GABA_level_subclass_annotation[match(colnames(budgerigar_5months_Tel),colnames(Tel_GABA))]),
                          as.character(budgerigar_5months_Tel$Tel_cell_type_res_1.2))

# Check for NA values in Tel_subclass_annotation column
na_indices <- is.na(budgerigar_5months_Tel$Tel_cell_type_res_1.2)

# Print the number of NA values found
cat("Number of NA values in Tel_subclass_annotation:", sum(na_indices), "\n")

# Remove cells with NA in Tel_subclass_annotation
budgerigar_5months_Tel <- subset(budgerigar_5months_Tel, cells = colnames(budgerigar_5months_Tel)[!na_indices])

# Verify that the NAs have been removed
cat("Number of NA values after filtering:", sum(is.na(budgerigar_5months_Tel$Tel_cell_type_res_1.2)), "\n")
budgerigar_5months_Tel$subclass_annotation <- budgerigar_5months_Tel$Tel_cell_type_res_1.2

##annotation_supertype(subcelltype) from Tel_cell_type_res_1.2
budgerigar_5months_Tel$Tel_cell_type_res_1.2 = ifelse(budgerigar_5months_Tel$Tel_cell_type_res_1.2=='Glutamatergic',
                          as.character(Tel_Glu$Glu_level_subcelltype_annotation[match(colnames(budgerigar_5months_Tel),colnames(Tel_Glu))]),
                          as.character(budgerigar_5months_Tel$Tel_cell_type_res_1.2))
budgerigar_5months_Tel$Tel_cell_type_res_1.2 = ifelse(budgerigar_5months_Tel$Tel_cell_type_res_1.2=='GABAergic',
                          as.character(Tel_GABA$GABA_level_subcelltype_annotation[match(colnames(budgerigar_5months_Tel),colnames(Tel_GABA))]),
                          as.character(budgerigar_5months_Tel$Tel_cell_type_res_1.2))

# Check for NA values in Tel_subclass_annotation column
na_indices <- is.na(budgerigar_5months_Tel$Tel_cell_type_res_1.2)

# Print the number of NA values found
cat("Number of NA values in Tel_subcelltype_annotation:", sum(na_indices), "\n")

# Remove cells with NA in Tel_subclass_annotation
budgerigar_5months_Tel <- subset(budgerigar_5months_Tel, cells = colnames(budgerigar_5months_Tel)[!na_indices])

# Verify that the NAs have been removed
cat("Number of NA values after filtering:", sum(is.na(budgerigar_5months_Tel$Tel_cell_type_res_1.2)), "\n")
budgerigar_5months_Tel$subcelltype_annotation <- budgerigar_5months_Tel$Tel_cell_type_res_1.2


#plot
cell_colors <- c(
'DACH2_GRIA4_Ex' = '#FF6B95',
'DACH2_GRIK1_Ex' = '#e1a4c6',    
'DACH2_MEIS2_Ex' = '#B2E5FB',
'DACH2_CCN2_Ex' = '#e3d1db', 
'CACNA1H_CCN2_Ex' = '#c7deef',
'SATB2_SOX6_Ex' = '#E45A5F',  
'SATB2_KIAA1217_Ex' = '#edeaa4', 
'SATB2_SCUBE1_Ex' = '#F19294',
'SATB2_SATB1_Ex' = '#684797',
'SATB2_RARB_Ex' = '#BDA7CB',
'RELN_LHX8_Inh' = '#7587b1',  
'MGE_SST_Inh' = '#4A9D47',
'MGE_PVALB_Inh' = '#d2e0ac',    
'MGE_ST18_Inh' = '#74a893',
'CGE_ADARB2_Inh' = '#C0937E',
'CGE_LAMP5_Inh' = "#5EC0AD",    
'LGE_PCP4_Inh' = '#e0bc58',
'LGE_FOXP2_Inh' = '#dc8e97',   
"Astrocyte" = "#82D4FB",
"Astroependymal" = "#33A1FF",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF"
)
DimPlot(budgerigar_5months_Tel, label = T, reduction = "umap",raster=FALSE,group.by ="subclass_annotation", alpha = 0.6,cols = cell_colors)
saveRDS(budgerigar_5months_Tel, file = '/data/work/seurat_integrated/09_budgerigar_brain_5months_Tel_annotated.rds')