library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")

Cere_GABA <- readRDS('/data/work/5months_Cere/GABA/10_budgerigar_brain_5months_Cere_GABA_unannotated.rds')

Idents(Cere_GABA) <- Cere_GABA$`GABA_integrated_snn_res.1.4`

main_type_anno <- c( 
'0' = 'Inh1', '1' = 'Inh2', '2' = 'Inh3', '3' = 'Inh4', '4' = 'Inh5', 
    '5' = 'Inh6', '6' = 'Inh7', '7' = 'Inh7', '8' = 'Inh8', '9' = 'Inh9', 
    '10' = 'Inh10', '11' = 'Inh11', '12' = 'Inh6', '13' = 'Inh12', '14' = 'Inh8', 
    '15' = 'Inh6', '16' = 'Inh13', '17' = 'Inh8')
                   

#names(main_type_anno) 
Cere_GABA <- RenameIdents(Cere_GABA, main_type_anno)
Cere_GABA$GABA_level_subcelltype <- Idents(Cere_GABA)

#delete low-quality
Cere_GABA <- subset(Cere_GABA, subset = GABA_level_subcelltype != "Inh2")

#supertype_annotation
Idents(Cere_GABA) <- Cere_GABA$`GABA_level_subcelltype`

main_type_anno <- c( 
               'Inh1' = 'SNCG_THY1_Inh', 'Inh3' = 'PTPRK_PVALB_Inh', 'Inh4' = 'PTPRK_PVALB_Inh',
    'Inh5' = 'FOXP2_Purkinje',
               'Inh6' = 'GRIA1_THY1_Inh', 'Inh7' = 'NXPH2_PVALB_Inh', 'Inh8' = 'CACNA1B_PAX2_Inh', 'Inh9' = 'NXPH1_RELN_Inh',
               'Inh10' = 'SNCG_THY1_Inh', 'Inh11' = 'PTPRK_RELN_Inh', 'Inh12' = 'GRIA1_THY1_Inh', 'Inh13' = 'SNCG_THY1_Inh'
)
                   

#names(main_type_anno) 
Cere_GABA <- RenameIdents(Cere_GABA, main_type_anno)
Cere_GABA$GABA_level_subcelltype_annotation <- Idents(Cere_GABA)
cell_colors <- c(
'SNCG_THY1_Inh' = '#807087', #FABDBB
'PTPRK_PVALB_Inh' = '#E45A5F',
'FOXP2_Purkinje' = '#e1a4c6',   
'GRIA1_THY1_Inh' = '#B2E5FB',
'NXPH2_PVALB_Inh' = '#c7deef',
'CACNA1B_PAX2_Inh' = '#4A9D47',
'NXPH1_RELN_Inh' = '#d2e0ac',    
"PTPRK_RELN_Inh" = "#FF33A1"
)
DimPlot(Cere_GABA, label = F, reduction = "umap",group.by ="GABA_level_subcelltype_annotation", cols = cell_colors)

#subclass_annotation
Idents(Cere_GABA) <- Cere_GABA$`GABA_integrated_snn_res.1.4`
main_type_anno <- c( 
'0' = 'Interneuron_1',  '2' = 'Interneuron_2', '3' = 'Interneuron_2', '4' = 'Purkinje', 
    '5' = 'Interneuron_1', '6' = 'Interneuron_2', '7' = 'Interneuron_2', '8' = 'Interneuron_2', '9' = 'Interneuron_2', 
    '10' = 'Interneuron_1', '11' = 'Interneuron_1','12' = 'Interneuron_1', '13' = 'Interneuron_1', '14' = 'Interneuron_2', 
    '15' = 'Interneuron_1', '16' = 'Interneuron_1', '17' = 'Interneuron_2')
                   

#names(main_type_anno) 
Cere_GABA <- RenameIdents(Cere_GABA, main_type_anno)
Cere_GABA$GABA_level_subclass_annotation <- Idents(Cere_GABA)

saveRDS(Cere_GABA, file = '/data/work/5month_Cere/GABA/10_Cere_GABA_annotated_!Inh2.rds')