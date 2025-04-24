library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
library('ggplot2')

Tel_GABA <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_GABA_unannotated.rds')
Tel_GABA

Idents(Tel_GABA) <- Tel_GABA$`GABA_integrated_snn_res.1.5`
main_type_anno2 <- c( 
'0' = 'Inh_1', '1' = 'Inh_2', '2' = 'Inh_3', '3' = 'Inh_1', '4' = 'Inh_4', 
    '5' = 'Inh_3', '6' = 'Inh_5', '7' = 'Inh_6', '8' = 'Inh_7', '9' = 'Inh_8', 
    '10' = 'Inh_7', '11' = 'Inh_9', '12' = 'Inh_10', '13' = 'Inh_7', '14' = 'Inh_11', 
    '15' = 'Inh_12', '16' = 'Inh_13', '17' = 'Inh_14', '18' = 'Inh_11', '19' = 'Inh_15', 
    '20' = 'Inh_16', '21' = 'Inh_9', '22' = 'Inh_15', '23' = 'Inh_17', '24' = 'Inh_18', 
    '25' = 'Inh_19', '26' = 'Inh_13', '27' = 'Inh_20', '28' = 'Inh_21', '29' = 'Inh_22', 
    '30' = 'Inh_23', '31' = 'Inh_24', '32' = 'Inh_25', '33' = 'Inh_26', '34' = 'Inh_27',
    '35' = 'Inh_28')
                   
Tel_GABA <- RenameIdents(Tel_GABA, main_type_anno2)
Tel_GABA$GABA_level_subcelltype <- Idents(Tel_GABA)

Tel_GABA <- subset(Tel_GABA, subset = GABA_level_subcelltype != "Inh_16")
Tel_GABA <- subset(Tel_GABA, subset = GABA_level_subcelltype != "Inh_19")
Tel_GABA <- subset(Tel_GABA, subset = GABA_level_subcelltype != "Inh_28")

##GABA_level_subcelltype(supertype)
Idents(Tel_GABA) <- Tel_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh_1' = 'LGE_PCP4_DRD2_Inh', 'Inh_2' = 'LGE_PCP4_DRD1_Inh', 'Inh_3' = 'LGE_FOXP2_SLIT2_Inh', 'Inh_4' = 'SNCG_NRP1_Inh','Inh_5' = 'MGE_SST_NPY_Inh',
               'Inh_6' = 'LGE_PCP4_EBF1_Inh', 'Inh_7' = 'MGE_SST_PVALB_Inh', 'Inh_8' = 'LGE_FOXP2_SLIT2_Inh', 'Inh_9' = 'CGE_ADARB2_NR2F2_Inh',
               'Inh_10' = 'MGE_SST_NPY_Inh', 'Inh_11' = 'MGE_SST_PVALB_Inh', 'Inh_12' = 'RELN_LHX8_Inh', 'Inh_13' = 'MGE_SST_NPY_Inh',
               'Inh_14' = 'CGE_ADARB2_NR2F2_Inh', 'Inh_15' = 'CGE_RELN_LAMP5_Inh', 'Inh_17' = 'MGE_SST_ST18_Inh',
               'Inh_18' = 'MGE_SST_ST18_Inh', 'Inh_20' = 'LGE_PCP4_SLC1A3_Inh', 'Inh_21' = 'MGE_SST_SATB1_Inh',
               'Inh_22' = 'LGE_PCP4_IL17REL_Inh', 'Inh_23' = 'RELN_LHX8_Inh', 'Inh_24' = 'CGE_ADARB2_NPY_Inh', 'Inh_25' = 'MGE_SST_SCUBE1_Inh',
               'Inh_26' = 'CGE_ADARB2_VIPR1_Inh', 'Inh_27' = 'RELN_LHX8_Inh'
)
                   

#names(main_type_anno) 
Tel_GABA <- RenameIdents(Tel_GABA, main_type_anno)
Tel_GABA$GABA_level_subcelltype_annotation <- Idents(Tel_GABA)



##GABA_level_subclass
Idents(Tel_GABA) <- Tel_GABA$`GABA_level_subcelltype`
main_type_anno <- c( 
               'Inh_1' = 'LGE_PCP4_Inh', 'Inh_2' = 'LGE_PCP4_Inh', 'Inh_3' = 'LGE_FOXP2_Inh', 'Inh_4' = 'RELN_LHX8_Inh','Inh_5' = 'MGE_SST_Inh',
               'Inh_6' = 'LGE_PCP4_Inh', 'Inh_7' = 'MGE_PVALB_Inh', 'Inh_8' = 'LGE_FOXP2_Inh', 'Inh_9' = 'CGE_ADARB2_Inh',
               'Inh_10' = 'MGE_SST_Inh', 'Inh_11' = 'MGE_PVALB_Inh', 'Inh_12' = 'RELN_LHX8_Inh', 'Inh_13' = 'MGE_SST_Inh',
               'Inh_14' = 'CGE_ADARB2_Inh', 'Inh_15' = 'CGE_LAMP5_Inh', 'Inh_17' = 'MGE_ST18_Inh',
               'Inh_18' = 'MGE_ST18_Inh', 'Inh_20' = 'LGE_PCP4_Inh', 'Inh_21' = 'MGE_PVALB_Inh',
               'Inh_22' = 'LGE_PCP4_Inh', 'Inh_23' = 'RELN_LHX8_Inh', 'Inh_24' = 'CGE_ADARB2_Inh', 'Inh_25' = 'MGE_SST_Inh',
               'Inh_26' = 'CGE_ADARB2_Inh', 'Inh_27' = 'RELN_LHX8_Inh'
)
                   

Tel_GABA <- RenameIdents(Tel_GABA, main_type_anno)
Tel_GABA$GABA_level_subclass_annotation <- Idents(Tel_GABA)




cell_colors <- c(
'LGE_PCP4_DRD2_Inh' = '#E45A5F', #FABDBB
'LGE_PCP4_DRD1_Inh' = '#807087',
'LGE_FOXP2_SLIT2_Inh' = '#e1a4c6',   
'SNCG_NRP1_Inh' = '#B2E5FB',
'LGE_PCP4_EBF1_Inh' = '#c7deef',
'MGE_SST_NPY_Inh' = '#EBEBF3',
'MGE_SST_PVALB_Inh' = '#60386D',    
'MGE_SST_SATB1_Inh' = '#907DAB',
'MGE_SST_SCUBE1_Inh' = '#A7B4C7',  
'MGE_SST_ST18_Inh' = '#e3d1db',
'CGE_ADARB2_NR2F2_Inh' = '#C0937E',
'RELN_LHX8_Inh' = '#684797',
'CGE_RELN_LAMP5_Inh' = '#BDA7CB',
'LGE_PCP4_IL17REL_Inh' = "#5EC0AD",   
'CGE_ADARB2_NPY_Inh' = '#e0bc58',
'CGE_ADARB2_VIPR1_Inh' = '#dc8e97', 
"LGE_PCP4_SLC1A3_Inh" = "#FEE082"
)
DimPlot(Tel_GABA, label = F, reduction = "umap",raster=FALSE,group.by ="GABA_level_subcelltype_annotation",, alpha = 0.9,cols = cell_colors)

##GABA_level_cluster
Idents(Tel_GABA) <- Tel_GABA$`GABA_integrated_snn_res.1.5`
main_type_anno <- c( 
'0' = 'Tel_GABA_1', '1' = 'Tel_GABA_2', '2' = 'Tel_GABA_3', '3' = 'Tel_GABA_4', '4' = 'Tel_GABA_5', 
    '5' = 'Tel_GABA_6', '6' = 'Tel_GABA_7', '7' = 'Tel_GABA_8', '8' = 'Tel_GABA_9', '9' = 'Tel_GABA_10', 
    '10' = 'Tel_GABA_11', '11' = 'Tel_GABA_12', '12' = 'Tel_GABA_13', '13' = 'Tel_GABA_14', '14' = 'Tel_GABA_15', 
    '15' = 'Tel_GABA_16', '16' = 'Tel_GABA_17', '17' = 'Tel_GABA_18', '18' = 'Tel_GABA_19', '19' = 'Tel_GABA_20', 
    '21' = 'Tel_GABA_21',  '22' = 'Tel_GABA_22', '23' = 'Tel_GABA_23', '24' = 'Tel_GABA_24', 
    '26' = 'Tel_GABA_25', '27' = 'Tel_GABA_26', '28' = 'Tel_GABA_27', '29' = 'Tel_GABA_28', 
    '30' = 'Tel_GABA_29', '31' = 'Tel_GABA_30', '32' = 'Tel_GABA_31', '33' = 'Tel_GABA_32', '34' = 'Tel_GABA_33')
                   

Tel_GABA <- RenameIdents(Tel_GABA, main_type_anno)
Tel_GABA$GABA_level_cluster <- Idents(Tel_GABA)


DimPlot(Tel_GABA, label = F, reduction = "umap",raster=FALSE,group.by ="GABA_level_subcelltype_annotation",, alpha = 0.9,cols = cell_colors)
ggsave("/data/work/5month_Tel/GABA/5months_Tel_GABA_subcelltype2.pdf",bg = "transparent", width =10, height = 8)

saveRDS(Tel_GABA, file = '/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_GABA_annotated.rds')