library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
library("ggplot2")
Th_GABA <- readRDS('/data/work/5month_Th/GABA/10_budgerigar_brain_5months_Th_GABA_unannotated.rds')
Th_GABA

Idents(Th_GABA) <- Th_GABA$`GABA_integrated_snn_res.3.5`
my4colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

VlnPlot(Th_GABA, features = c('SLC17A6','GAD1', 'GAD2','TH'), 
             stack = TRUE, 
             sort = TRUE, 
             cols = my4colors,
             split.by =  "GABA_integrated_snn_res.3.5" , #每种cluster 一个颜色
             flip = TRUE) +
  theme(legend.position = "none") + 
  ggtitle("Identity on x-axis")
ggsave("/data/work/5month_Th/GABA/GABA_integrated_snn_res.3.5_vlnplot.pdf", width = 15, height = 10)

#delete Glu
Th_GABA <- subset(Th_GABA, subset = GABA_integrated_snn_res.3.5 != "3")
Th_GABA <- subset(Th_GABA, subset = GABA_integrated_snn_res.3.5 != "14")

Idents(Th_GABA) <- Th_GABA$`GABA_integrated_snn_res.3.5`
main_type_anno <- c( 
   '0' = 'Inh1', '1' = 'Inh2', '2' = 'Inh3', '4' = 'Inh4', 
    '5' = 'Inh5', '6' = 'Inh6', '7' = 'Inh7', '8' = 'Inh8','9' = 'Inh9', 
    '10' = 'Inh10', '11' = 'Inh8', '12' = 'Inh2', '13' = 'Inh11', 
    '15' = 'Inh7','16' = 'Inh12', '17' = 'Inh13', '18' = 'Inh13', '19' = 'Inh14', 
    '20' = 'Inh1', '21' = 'Inh15', 
    '22' = 'Inh1', '23' = 'Inh1', '24' = 'Inh16', 
    '25' = 'Inh17', '26' = 'Inh18', '27' = 'Inh19', '28' = 'Inh8', '29' = 'Inh20', 
    '30' = 'Inh21', '31' = 'Inh18', '32' = 'Inh8'
)
                   

#names(main_type_anno) 
Th_GABA <- RenameIdents(Th_GABA, main_type_anno)
Th_GABA$GABA_level_subcelltype <- Idents(Th_GABA)

#supertype_annotation
Idents(Th_GABA) <- Th_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh1' = 'ENSMUNG00000008668_Inh', 'Inh2' = 'SNCG_NEFL_Inh', 'Inh3' = 'SNCG_NEFL_Inh', 'Inh4' = 'ENSMUNG00000008668_Inh','Inh5' = 'CHGB_GATA3_Inh',
               'Inh6' = 'NHP2_HMOX2_Inh', 'Inh7' = 'ISL1_MEIS2_Inh', 'Inh8' = 'NHP2_HMOX2_Inh', 'Inh9' = 'SCN4B_PVALB_Inh',
               'Inh10' = 'PAX6_MEIS2_Inh', 'Inh11' = 'PAX2_PAX5_Inh', 'Inh12' = 'DLX1_NDNF_Inh', 'Inh13' = 'CEMIP_MEIS2_Inh',
               'Inh14' = 'NHP2_HMOX2_Inh', 'Inh15' = 'PENK_GATA3_Inh','Inh16' = 'CHGB_GATA3_Inh', 'Inh17' = 'LHX8_Inh',
               'Inh18' = 'ENSMUNG00000008668_Inh','Inh19' = 'CADPS2_Inh', 'Inh20' = 'NPY_GATA3_Inh', 'Inh21' = 'ST18_Inh'
)
                   

#names(main_type_anno) 
Th_GABA <- RenameIdents(Th_GABA, main_type_anno)
Th_GABA$GABA_level_subcelltype_annotation <- Idents(Th_GABA)

cell_colors <- c(
'ENSMUNG00000008668_Inh' = '#807087', #FABDBB
'SNCG_NEFL_Inh' = '#E45A5F',
'CHGB_GATA3_Inh' = '#e1a4c6',   
'NHP2_HMOX2_Inh' = '#B2E5FB',
'ISL1_MEIS2_Inh' = '#c7deef',
'SCN4B_PVALB_Inh' = '#4A9D47',
'PAX6_MEIS2_Inh' = '#d2e0ac',    
'PAX2_PAX5_Inh' = '#74a893',
'DLX1_NDNF_Inh' = '#edeaa4',  
'CEMIP_MEIS2_Inh' = '#e3d1db',
'PENK_GATA3_Inh' = '#C0937E',
'LHX8_Inh' = '#684797',
'CADPS2_Inh' = '#BDA7CB',
 'NPY_GATA3_Inh' = "#5EC0AD",   
'ST18_Inh' = '#e0bc58'
)
DimPlot(Th_GABA, reduction = "umap", label = F,group.by= 'GABA_level_subcelltype_annotation',cols =cell_colors)

#subclass_annotation
Idents(Th_GABA) <- Th_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh1' = 'ENSMUNG00000008668_Inh', 'Inh2' = 'SNCG_GATA3_Inh', 'Inh3' = 'SNCG_GATA3_Inh', 'Inh4' = 'ENSMUNG00000008668_Inh','Inh5' = 'SNCG_GATA3_Inh',
               'Inh6' = 'NHP2_HMOX2_Inh', 'Inh7' = 'ISL1_MEIS2_Inh', 'Inh8' = 'NHP2_HMOX2_Inh', 'Inh9' = 'SNCG_GATA3_Inh',
               'Inh10' = 'DLX1_MEIS2_Inh', 'Inh11' = 'SNCG_GATA3_Inh', 'Inh12' = 'DLX1_MEIS2_Inh', 'Inh13' = 'CEMIP_MEIS2_Inh',
               'Inh14' = 'NHP2_HMOX2_Inh', 'Inh15' = 'SNCG_GATA3_Inh','Inh16' = 'SNCG_GATA3_Inh', 'Inh17' = 'ENSMUNG00000008668_Inh',
               'Inh18' = 'ENSMUNG00000008668_Inh','Inh19' = 'ENSMUNG00000008668_Inh', 'Inh20' = 'SNCG_GATA3_Inh', 'Inh21' = 'ENSMUNG00000008668_Inh'
)
                   

#names(main_type_anno) 
Th_GABA <- RenameIdents(Th_GABA, main_type_anno)
Th_GABA$GABA_level_subclass_annotation <- Idents(Th_GABA)