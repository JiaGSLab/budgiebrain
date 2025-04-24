library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")

Cere_Glu <- readRDS('/data/work/5month_Cere/Glu/10_budgerigar_brain_5months_Cere_Glu_unannotated.rds')
Idents(Cere_Glu) <- Cere_Glu$`Glu_integrated_snn_res.1.2`
Cere_Glu <- subset(Cere_Glu, subset = Glu_integrated_snn_res.1.2 != "0")
Cere_Glu <- subset(Cere_Glu, subset = Glu_integrated_snn_res.1.2 != "18")



#supertype
Idents(Cere_Glu) <- Cere_Glu$`Glu_integrated_snn_res.1.2`
main_type_anno <- c( 
'1' = 'Granule2', '2' = 'BDNF_Granule1', '3' = 'BDNF_Granule1', '4' = 'Granule2', 
    '5' = 'Granule2', '6' = 'BDNF_Granule1', '7' = 'TRPC4_UBC', '8' = 'SNCG_Ex', '9' = 'Granule2', 
    '10' = 'SNCG_Ex', '11' = 'Granule2', '12' = 'TRPC4_UBC', '13' = 'TRPC4_UBC', '14' = 'Granule2', 
    '15' = 'Granule2', '16' = 'SNCG_Ex', '17' = 'Granule2' )
#names(main_type_anno) 
Cere_Glu <- RenameIdents(Cere_Glu, main_type_anno)
Cere_Glu$Glu_level_subcelltype_annotation <- Idents(Cere_Glu)
cell_colors <- c(
'BDNF_Granule1' = '#807087', #FABDBB
'Granule2' = '#E45A5F',
'TRPC4_UBC' = '#e1a4c6',   
'SNCG_Ex' = '#B2E5FB')
DimPlot(Cere_Glu, label = F, reduction = "umap",group.by ="Glu_level_subcelltype_annotation", alpha = 0.5,cols = cell_colors)



#subclass_annotation
Idents(Cere_Glu) <- Cere_Glu$`Glu_integrated_snn_res.1.2`
main_type_anno <- c( 
'1' = 'Granule', '2' = 'Granule', '3' = 'Granule', '4' = 'Granule', 
    '5' = 'Granule', '6' = 'Granule', '7' = 'TRPC4_UBC', '8' = 'SNCG_Ex', '9' = 'Granule', 
    '10' = 'SNCG_Ex', '11' = 'Granule', '12' = 'TRPC4_UBC', '13' = 'TRPC4_UBC', '14' = 'Granule', 
    '15' = 'Granule', '16' = 'SNCG_Ex', '17' = 'Granule' )
#names(main_type_anno) 
Cere_Glu <- RenameIdents(Cere_Glu, main_type_anno)
Cere_Glu$Glu_level_subclass_annotation <- Idents(Cere_Glu)

saveRDS(Cere_Glu, file = '/data/work/5month_Cere/Glu/10_Cere_Glu_annotated_!Ex1_Ex14.rds')