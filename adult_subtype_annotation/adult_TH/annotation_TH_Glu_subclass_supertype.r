library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
Th_Glu <- readRDS('/data/work/5month_Th/Glu/10_budgerigar_brain_5months_Th_Glu_unannotated.rds')
Th_Glu

#delete GABA
Idents(Th_Glu) <- Th_Glu$`Glu_integrated_snn_res.4.5`
Th_GABA <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 %in% c('7', '8', '15', '38'))
saveRDS(Th_GABA, file = '/data/work/5month_Th/Th_GABA_subset.rds')
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "7")
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "8")
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "15")
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "38")
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "35")
Th_Glu <- subset(Th_Glu, subset = Glu_integrated_snn_res.4.5 != "42")

Idents(Th_Glu) <- Th_Glu$`Glu_integrated_snn_res.4.5`
main_type_anno <- c( 
'0' = 'Ex1', '1' = 'Ex2', '2' = 'Ex3', '3' = 'Ex4', '4' = 'Ex2', 
    '5' = 'Ex4', '6' = 'Ex5', '9' = 'Ex5', 
    '10' = 'Ex2', '11' = 'Ex4', '12' = 'Ex4', '13' = 'Ex6', '14' = 'Ex1', 
    '16' = 'Ex7', '17' = 'Ex8', '18' = 'Ex1', '19' = 'Ex5', 
    '20' = 'Ex6', '21' = 'Ex9', 
    '22' = 'Ex7', '23' = 'Ex8', '24' = 'Ex8', 
    '25' = 'Ex10', '26' = 'Ex11', '27' = 'Ex12', '28' = 'Ex2', '29' = 'Ex10', 
    '30' = 'Ex4', '31' = 'Ex13', '32' = 'Ex11', '33' = 'Ex14', '34' = 'Ex15',
     '36' = 'Ex2', '37' = 'Ex16', '39' = 'Ex3', 
    '40' = 'Ex2', '41' = 'Ex11', '43' = 'Ex17'
)
#names(Glu_level_subcelltype) 
Th_Glu <- RenameIdents(Th_Glu, main_type_anno)
Th_Glu$Glu_level_subcelltype <- Idents(Th_Glu)

#delete low-quality cluster
Th_Glu <- subset(Th_Glu, subset = Glu_level_subcelltype != "Ex3")

#supertype_annotation
Idents(Th_Glu) <- Th_Glu$`Glu_level_subcelltype`
main_type_anno <- c( 
               'Ex1' = 'SLIT1_NR2F2_Ex', 'Ex2' = 'KIAA1217_ITPR1_Ex',  'Ex4' = 'VIP_NTS_Ex',
               'Ex5' = 'VIP_MEGF11_Ex', 'Ex6' = 'CCK_FOXP2_Ex','Ex7' = 'CACNA1G_SCN4B_Ex',
               'Ex8' = 'SLIT1_NR2F2_Ex',  'Ex9' = 'GRIA4_MEGF11_Ex',
               'Ex10' = 'CCK_FOXP2_Ex', 'Ex11' = 'SNCG_NEFL_Ex', 'Ex12' = 'ZNF385B_RORB_Ex', 'Ex13' = 'CALCB_ETV1_Ex',
               'Ex14' = 'SLIT1_NR2F2_Ex', 'Ex15' = 'CCN2_NR4A2_Ex', 'Ex16' = 'SLIT1_NR2F2_Ex', 'Ex17' = 'TAC1_SCUBE1_Ex'
)
                   

#names
Th_Glu <- RenameIdents(Th_Glu, main_type_anno)
Th_Glu$Glu_level_subcelltype_annotation <- Idents(Th_Glu)

cell_colors <- c(
'SLIT1_NR2F2_Ex' = '#E45A5F', 
'KIAA1217_ITPR1_Ex' = '#807087',
'VIP_NTS_Ex' = '#e1a4c6',   
'VIP_MEGF11_Ex' = '#B2E5FB',
'CCK_FOXP2_Ex' = '#c7deef',
'CACNA1G_SCN4B_Ex' = '#4A9D47',
'GRIA4_MEGF11_Ex' = '#d2e0ac',    
'SNCG_NEFL_Ex' = '#74a893',
'ZNF385B_RORB_Ex' = '#edeaa4',  
'CALCB_ETV1_Ex' = '#e3d1db',
'CCN2_NR4A2_Ex' = '#C0937E',
'TAC1_SCUBE1_Ex' = '#684797'
)
DimPlot(Th_Glu, label = F, reduction = "umap",raster=FALSE,group.by ="Glu_level_subcelltype_annotation",alpha = 0.9,cols = cell_colors)


Idents(Th_Glu) <- Th_Glu$`Glu_level_subcelltype`
main_type_anno <- c( 
               'Ex1' = 'SLIT1_PTPRF_Ex', 'Ex2' = 'KIAA1217_ITPR1_Ex',  'Ex4' = 'VIP_NTS_Ex',
               'Ex5' = 'VIP_NTS_Ex', 'Ex6' = 'CCK_FOXP2_Ex','Ex7' = 'VIP_NTS_Ex',
               'Ex8' = 'SLIT1_PTPRF_Ex',  'Ex9' = 'GRIA3_GRIK1_Ex',
               'Ex10' = 'CCK_FOXP2_Ex', 'Ex11' = 'SNCG_NEFL_Ex', 'Ex12' = 'GRIA3_GRIK1_Ex', 'Ex13' = 'CCK_FOXP2_Ex',
               'Ex14' = 'SLIT1_PTPRF_Ex', 'Ex15' = 'SLIT1_PTPRF_Ex', 'Ex16' = 'SLIT1_PTPRF_Ex', 'Ex17' = 'SLIT1_PTPRF_Ex'
)
                   

#names(main_type_anno) 
Th_Glu <- RenameIdents(Th_Glu, main_type_anno)
Th_Glu$Glu_level_subclass_annotation <- Idents(Th_Glu)
saveRDS(Th_Glu, file = '/data/work/5month_Th/Glu/10_budgerigar_brain_5months_Th_Glu_annotated!GABA_Ex3.rds')