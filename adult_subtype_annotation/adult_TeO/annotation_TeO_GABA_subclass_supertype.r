library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
TeO_GABA <- readRDS('/data/work/5month_TeO/10_budgerigar_brain_5months_TeO_GABA_unannotated.rds')
TeO_GABA

Idents(TeO_GABA) <- TeO_GABA$`GABA_integrated_snn_res.2`

main_type_anno <- c( 
'0' = 'Inh1', '1' = 'Inh2', '2' = 'Inh3', '3' = 'Inh4', '4' = 'Inh5', 
    '5' = 'Inh6', '6' = 'Inh7', '7' = 'Inh8', '8' = 'Inh9', '9' = 'Inh10', 
    '10' = 'Inh11', '11' = 'Inh12', '12' = 'Inh8', '13' = 'Inh13', '14' = 'Inh14', 
    '15' = 'Inh15', '16' = 'Inh16', '17' = 'Inh13', '18' = 'Inh17', '19' = 'Inh18', 
    '20' = 'Inh19', '21' = 'Inh10', '22' = 'Inh11', '23' = 'Inh4', '24' = 'Inh20', 
    '25' = 'Inh21', '26' = 'Inh7', '27' = 'Inh9', '28' = 'Inh6', '29' = 'Inh5', 
    '30' = 'Inh19', '31' = 'Inh22', '32' = 'Inh23', '33' = 'Inh7', '34' = 'Inh6',
    '35' = 'Inh24','36' = 'Inh25')
                   

#names(main_type_anno) 
TeO_GABA <- RenameIdents(TeO_GABA, main_type_anno)
TeO_GABA$GABA_level_subcelltype <- Idents(TeO_GABA)

TeO_GABA <- subset(TeO_GABA, subset = GABA_level_subcelltype != "Inh3")
TeO_GABA <- subset(TeO_GABA, subset = GABA_level_subcelltype != "Inh5")
TeO_GABA <- subset(TeO_GABA, subset = GABA_level_subcelltype != "Inh19")

#supertype_annotation
Idents(TeO_GABA) <- TeO_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh1' = 'ZFHX4_CHL1_Inh', 'Inh2' = 'MEIS1_ST18_Inh',  'Inh4' = 'MEIS2_CALB2_Inh',
               'Inh6' = 'IRX1_IRX2_Inh', 'Inh7' = 'CADPS2_GATA3_Inh', 'Inh8' = 'MEIS2_NDNF_Inh', 'Inh9' = 'MEIS2_RELN_Inh',
               'Inh10' = 'GATA3_DACH1_Inh', 'Inh11' = 'EMX2_PENK_Inh', 'Inh12' = 'ZFHX4_PAX7_Inh', 'Inh13' = 'PVALB_EBF1_Inh',
               'Inh14' = 'PVALB_LYPD1_Inh', 'Inh15' = 'PVALB_EBF1_Inh', 'Inh16' = 'MEIS2_PAX7_Inh','Inh17' = 'RELN_FOXP2_Inh',
               'Inh18' = 'ADARB2_GRIK1_Inh', 'Inh20' = 'PCP4_CADPS2_Inh', 'Inh21' = 'MEIS2_RELN_Inh',
               'Inh22' = 'ZFHX4_GATA3_Inh', 'Inh23' = 'IRX1_IRX2_Inh', 'Inh24' = 'GATA3_LYPD1_Inh', 'Inh25' = 'GATA3_NPY_Inh'
)
                   

#names(main_type_anno) 
TeO_GABA <- RenameIdents(TeO_GABA, main_type_anno)
TeO_GABA$GABA_level_subcelltype_annotation <- Idents(TeO_GABA)

cell_colors <- c(
'ZFHX4_CHL1_Inh' = '#FF33A1',
'MEIS1_ST18_Inh' = '#E45A5F',
'MEIS2_CALB2_Inh' = '#e1a4c6',   
'IRX1_IRX2_Inh' = '#B2E5FB',
'CADPS2_GATA3_Inh' = '#c7deef',
'MEIS2_NDNF_Inh' = '#4A9D47',   
'MEIS2_RELN_Inh' = '#d2e0ac',    
'GATA3_DACH1_Inh' = '#74a893',
'EMX2_PENK_Inh' = '#7587b1',   
'ZFHX4_PAX7_Inh' = '#e3d1db',
'PVALB_EBF1_Inh' = '#F19294',   
'PVALB_LYPD1_Inh' = '#C0937E',    
'MEIS2_PAX7_Inh' = '#684797',
'RELN_FOXP2_Inh' = '#edeaa4',      
'ADARB2_GRIK1_Inh' = '#BDA7CB', 
'PCP4_CADPS2_Inh' = '#e0bc58',  
'ZFHX4_GATA3_Inh' = '#dc8e97', 
"GATA3_LYPD1_Inh" = "#82D4FB",
"GATA3_NPY_Inh" = "#807087"
)
DimPlot(TeO_GABA, label = F, reduction = "umap",group.by ="GABA_level_subcelltype_annotation",cols = cell_colors)

#subclass_annotation
Idents(TeO_GABA) <- TeO_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh1' = 'ZFHX4_CHL1_Inh', 'Inh2' = 'SHOX2_SATB1_Inh',  'Inh4' = 'CEMIP_RELN_Inh',
               'Inh6' = 'IRX1_IRX2_Inh', 'Inh7' = 'GATA3_DACH1_Inh', 'Inh8' = 'MEIS2_MOXD1_Inh', 'Inh9' = 'CEMIP_RELN_Inh',
               'Inh10' = 'GATA3_DACH1_Inh', 'Inh11' = 'CA8_KIAA1217_Inh', 'Inh12' = 'ZFHX4_CHL1_Inh', 'Inh13' = 'EBF1_EBF3_Inh',
               'Inh14' = 'EBF1_EBF3_Inh', 'Inh15' = 'EBF1_EBF3_Inh', 'Inh16' = 'MEIS2_MOXD1_Inh','Inh17' = 'CEMIP_RELN_Inh',
               'Inh18' = 'ZFHX4_CHL1_Inh', 'Inh20' = 'SHOX2_SATB1_Inh', 'Inh21' = 'CEMIP_RELN_Inh',
               'Inh22' = 'ZFHX4_CHL1_Inh', 'Inh23' = 'IRX1_IRX2_Inh', 'Inh24' = 'GATA3_DACH1_Inh', 'Inh25' = 'ZFHX4_CHL1_Inh'
)
                   

#names(main_type_anno) 
TeO_GABA <- RenameIdents(TeO_GABA, main_type_anno)
TeO_GABA$GABA_level_subclass_annotation <- Idents(TeO_GABA)

saveRDS(TeO_GABA, file ='/data/work/5month_TeO/GABA/10_TeO_GABA_annotated_!Inh3_5_19.rds')