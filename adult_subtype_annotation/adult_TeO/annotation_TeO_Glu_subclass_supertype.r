library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")

TeO_Glu <- readRDS('/data/work/5month_TeO/10_budgerigar_brain_5months_TeO_Glu_unannotated.rds')

Idents(TeO_Glu) <- TeO_Glu$`Glu_integrated_snn_res.4`

main_type_anno <- c( 
'0' = 'Ex1', '1' = 'Ex2', '2' = 'Ex3', '3' = 'Ex4', '4' = 'Ex5', 
    '5' = 'Ex6', '6' = 'Ex7', '7' = 'Ex8', '8' = 'Ex9', '9' = 'Ex10', 
    '10' = 'Ex11', '11' = 'Ex12', '12' = 'Ex13', '13' = 'Ex14', '14' = 'Ex14', 
    '15' = 'Ex13', '16' = 'Ex11', '17' = 'Ex13', '18' = 'Ex13', '19' = 'Ex15', 
    '20' = 'Ex9', '21' = 'Ex8', '22' = 'Ex16', '23' = 'Ex7', '24' = 'Ex16', 
    '25' = 'Ex5', '26' = 'Ex17', '27' = 'Ex18', '28' = 'Ex19', '29' = 'Ex20', 
    '30' = 'Ex21', '31' = 'Ex22', '32' = 'Ex23', '33' = 'Ex24', '34' = 'Ex25',
    '35' = 'Ex15','36' = 'Ex26', '37' = 'Ex14', '38' = 'Ex27', '39' = 'Ex19', 
    '40' = 'Ex25', '41' = 'Ex10', '42' = 'Ex28', '43' = 'Ex17', '44' = 'Ex17',
    '45' = 'Ex29')
                   

#names(main_type_anno) 
TeO_Glu <- RenameIdents(TeO_Glu, main_type_anno)
TeO_Glu$Glu_level_subcelltype <- Idents(TeO_Glu)

TeO_Glu <- subset(TeO_Glu, subset = Glu_level_subcelltype != "Ex3")
TeO_Glu <- subset(TeO_Glu, subset = Glu_level_subcelltype != "Ex6")

Idents(TeO_Glu) <- TeO_Glu$`Glu_level_subcelltype`

main_type_anno <- c( 
               'Ex1' = 'CCK_ZFHX4_Ex', 'Ex2' = 'SLIT2_ZFHX4_Ex', 'Ex4' = 'SLIT2_ZFHX4_Ex','Ex5' = 'KIAA1217_CACNA1G_Ex',
               'Ex7' = 'VIP_SOX6_Ex', 'Ex8' = 'NTS_ZFHX4_Ex', 'Ex9' = 'CABP7_CACNA1G_Ex',
               'Ex10' = 'Stressed_Ex', 'Ex11' = 'KIAA1217_GRIK1_Ex', 'Ex12' = 'Stressed_Ex', 'Ex13' = 'NXPH1_POSTN_Ex',
               'Ex14' = 'CABP7_MEIS2_Ex', 'Ex15' = 'SLIT2_ETV1_Ex', 'Ex16' = 'PENK_SCN4B_Ex', 'Ex17' = 'CCK_ZFHX4_Ex',
               'Ex18' = 'NTNG2_CACNA1G_Ex', 'Ex19' = 'LYPD1_SCN4B_Ex', 'Ex20' = 'SLIT2_FXYD6_Ex', 'Ex21' = 'SPON1_MEIS2_Ex',
               'Ex22' = 'SST_CACNA1G_Ex', 'Ex23' = 'GRIA4_ETV1_Ex', 'Ex24' = 'FOXP2_TCF7L2_Ex', 'Ex25' = 'ZIC1_ZIC2_Ex',
               'Ex26' = 'SPON1_TAC1_Ex', 'Ex27' = 'RSPO3_SPON1_Ex', 'Ex28' = 'NR4A2_TCF7L2_Ex', 'Ex29' = 'SLIT2_ZFHX4_Ex'
)
                   

#names(main_type_anno) 
TeO_Glu <- RenameIdents(TeO_Glu, main_type_anno)
TeO_Glu$Glu_level_subcelltype_annotation <- Idents(TeO_Glu)
cell_colors <- c(
'CCK_ZFHX4_Ex' = '#FF33A1',
'SLIT2_ZFHX4_Ex' = '#E45A5F',
'KIAA1217_CACNA1G_Ex' = '#e1a4c6',   
'VIP_SOX6_Ex' = '#B2E5FB',
'NTS_ZFHX4_Ex' = '#c7deef',
'CABP7_CACNA1G_Ex' = '#4A9D47',
'Stressed_Ex' = '#d2e0ac',    
'KIAA1217_GRIK1_Ex' = '#74a893',
'NXPH1_POSTN_Ex' = '#7587b1',   
'CABP7_MEIS2_Ex' = '#e3d1db',
'SLIT2_ETV1_Ex' = '#F19294',
'PENK_SCN4B_Ex' = '#C0937E',
'NTNG2_CACNA1G_Ex' = '#684797',
'LYPD1_SCN4B_Ex' = '#edeaa4',  
'SLIT2_FXYD6_Ex' = '#BDA7CB',
 'SPON1_MEIS2_Ex' = "#5EC0AD",   
'SST_CACNA1G_Ex' = '#e0bc58',
'GRIA4_ETV1_Ex' = '#dc8e97', 
"FOXP2_TCF7L2_Ex" = "#82D4FB",
"ZIC1_ZIC2_Ex" = "#807087",
"SPON1_TAC1_Ex" = "#FEE082",
"RSPO3_SPON1_Ex" = "#C7B889",
"NR4A2_TCF7L2_Ex" = "#8C33FF"
)
DimPlot(TeO_Glu, label = F, reduction = "umap",raster=F,group.by ="Glu_level_subcelltype_annotation",cols = cell_colors)


Idents(TeO_Glu) <- TeO_Glu$`Glu_level_subcelltype`
main_type_anno <- c( 
               'Ex1' = 'ZFHX4_SLIT1_Ex', 'Ex2' = 'ZFHX4_SLIT1_Ex', 'Ex4' = 'ZFHX4_SLIT1_Ex','Ex5' = 'KIAA1217_CACNA1G_Ex',
               'Ex7' = 'CDH9_TLL1_Ex', 'Ex8' = 'ZFHX4_SLIT1_Ex', 'Ex9' = 'SCN4B_RGS7_Ex',
               'Ex10' = 'ZFHX4_SLIT1_Ex', 'Ex11' = 'KIAA1217_CACNA1G_Ex', 'Ex12' = 'ZFHX4_SLIT1_Ex', 'Ex13' = 'POSTN_NXPH1_Ex',
               'Ex14' = 'MEIS2_SVIL_Ex', 'Ex15' = 'CDH9_TLL1_Ex', 'Ex16' = 'SCN4B_RGS7_Ex', 'Ex17' = 'ZFHX4_SLIT1_Ex',
               'Ex18' = 'SCN4B_RGS7_Ex', 'Ex19' = 'SCN4B_RGS7_Ex', 'Ex20' = 'ZFHX4_SLIT1_Ex', 'Ex21' = 'SCN4B_RGS7_Ex',
               'Ex22' = 'POSTN_NXPH1_Ex', 'Ex23' = 'SCN4B_RGS7_Ex', 'Ex24' = 'SCN4B_RGS7_Ex', 'Ex25' = 'ZFHX4_SLIT1_Ex',
               'Ex26' = 'ZFHX4_SLIT1_Ex', 'Ex27' = 'SCN4B_RGS7_Ex', 'Ex28' = 'ZFHX4_SLIT1_Ex', 'Ex29' = 'ZFHX4_SLIT1_Ex'
)
                   

#names(main_type_anno) 
TeO_Glu <- RenameIdents(TeO_Glu, main_type_anno)
TeO_Glu$Glu_level_subclass_annotation <- Idents(TeO_Glu)

saveRDS(TeO_Glu, file = '/data/work/5month_TeO/Glu/10_budgerigar_brain_5months_TeO_Glu_annotated_!Ex3_6.rds')