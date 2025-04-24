library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
library('ggplot2')


Tel_Glu <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_Glu_unannotated.rds')

Idents(Tel_Glu) <- Tel_Glu$`Glu_integrated_snn_res.4`
main_type_anno <- c( 
'0' = 'Ex1', '1' = 'Ex1', '2' = 'Ex2', '3' = 'Ex3', '4' = 'Ex4', 
    '5' = 'Ex5', '6' = 'Ex4', '7' = 'Ex2', '8' = 'Ex6', '9' = 'Ex7', 
    '10' = 'Ex8', '11' = 'Ex5', '12' = 'Ex6', '13' = 'Ex9', '14' = 'Ex7', 
    '15' = 'Ex1', '16' = 'Ex10', '17' = 'Ex11', '18' = 'Ex12', '19' = 'Ex13', 
    '20' = 'Ex14', '21' = 'Ex1', '22' = 'Ex15', '23' = 'Ex16', '24' = 'Ex17', 
    '25' = 'Ex18', '26' = 'Ex19', '27' = 'Ex16', '28' = 'Ex6', '29' = 'Ex20', 
    '30' = 'Ex21', '31' = 'Ex22', '32' = 'Ex7', '33' = 'Ex23', '34' = 'Ex10',
    '35' = 'Ex15', '36' = 'Ex24', 
    '37' = 'Ex25', '38' = 'Ex10', '39' = 'Ex3', 
    '40' = 'Ex23', '41' = 'Ex1', '42' = 'Ex24', '43' = 'Ex7', '44' = 'Ex17', 
    '45' = 'Ex15', '46' = 'Ex26', '47' = 'Ex22', '48' = 'Ex6', '49' = 'Ex27', 
    '50' = 'Ex24', '51' = 'Ex21', '52' = 'Ex16', '53' = 'Ex28', '54' = 'Ex8', 
    '55' = 'Ex29','56' = 'Ex8'
)
                   

#names(main_type_anno) 
Tel_Glu <- RenameIdents(Tel_Glu, main_type_anno)
Tel_Glu$Glu_level_celltype_total <- Idents(Tel_Glu)

# delete Ex5 
Tel_Glu <- subset(Tel_Glu, subset = Glu_level_celltype_total != "Ex5")
table(Tel_Glu$Glu_level_celltype_total)


#Glu_level_subcelltype_annotation
Idents(Tel_Glu) <- Tel_Glu$`Glu_level_celltype_total`

main_type_anno <- c( 
               'Ex1' = 'DACH2_CHD7_Ex', 'Ex2' = 'SATB2_FOXP2_Ex', 'Ex3' = 'DACH2_SV2C_Ex', 'Ex4' = 'DACH2_NDST4_Ex',
               'Ex6' = 'DACH2_LYPD1_Ex', 'Ex7' = 'DACH2_GRIK1_Ex', 'Ex8' = 'DACH2_CCN2_Ex', 'Ex9' = 'SATB2_BCL6_Ex',
               'Ex10' = 'DACH2_MEIS2_Ex', 'Ex11' = 'CACNA1H_SLIT1_Ex', 'Ex12' = 'SATB2_CCK_Ex', 'Ex13' = 'DACH2_BDNF_Ex',
               'Ex14' = 'SATB2_TFAP2D_Ex', 'Ex15' = 'SATB2_KIAA1217_Ex', 'Ex16' = 'SATB2_SCUBE1_Ex', 'Ex17' = 'SATB2_ZNF385B_Ex',
               'Ex18' = 'SATB2_RARB_Ex', 'Ex19' = 'SATB2_TAC1_Ex', 'Ex20' = 'CACNA1H_NRP1_Ex', 'Ex21' = 'DACH2_ZNF381_Ex',
               'Ex22' = 'SATB2_GLRA3_Ex', 'Ex23' = 'CACNA1H_SLIT1_Ex', 'Ex24' = 'SATB2_ZNF385B_Ex', 'Ex25' = 'DACH2_MEIS2_Ex',
               'Ex26' = 'CACNA1H_ZBTB20_Ex', 'Ex27' = 'SATB2_KIAA1217_Ex', 'Ex28' = 'SATB2_KIAA1217_Ex', 'Ex29' = 'DACH2_TAC1_Ex'
)
                   

#names(main_type_anno) 
Tel_Glu <- RenameIdents(Tel_Glu, main_type_anno)
Tel_Glu$Glu_level_subcelltype_annotation <- Idents(Tel_Glu)

cell_colors <- c(
'DACH2_CHD7_Ex' = '#807087', #FABDBB
'SATB2_FOXP2_Ex' = '#E45A5F',
'DACH2_SV2C_Ex' = '#e1a4c6',   
'DACH2_NDST4_Ex' = '#B2E5FB',
'DACH2_LYPD1_Ex' = '#c7deef',
'DACH2_GRIK1_Ex' = '#4A9D47',
'DACH2_CCN2_Ex' = '#d2e0ac',    
'SATB2_BCL6_Ex' = '#74a893',
'DACH2_MEIS2_Ex' = '#7587b1',   
'CACNA1H_SLIT1_Ex' = '#edeaa4',  
'SATB2_CCK_Ex' = '#e3d1db',
'DACH2_BDNF_Ex' = '#F19294',
'SATB2_TFAP2D_Ex' = '#C0937E',
'SATB2_KIAA1217_Ex' = '#684797',
'SATB2_SCUBE1_Ex' = '#BDA7CB',
 'SATB2_ZNF385B_Ex' = "#5EC0AD",   
'SATB2_RARB_Ex' = '#e0bc58',
'SATB2_TAC1_Ex' = '#dc8e97', 
"CACNA1H_NRP1_Ex" = "#82D4FB",
"DACH2_ZNF381_Ex" = "#FF33A1",
"SATB2_GLRA3_Ex" = "#FEE082",
"CACNA1H_ZBTB20_Ex" = "#C7B889",
"DACH2_TAC1_Ex" = "#8C33FF"
)
DimPlot(Tel_Glu, label = F, reduction = "umap",raster=FALSE,group.by ="Glu_level_subcelltype_annotation", alpha = 0.6,cols = cell_colors)



#Glu_level_subclass_annotation
Idents(Tel_Glu) <- Tel_Glu$`Glu_level_celltype_total`
main_type_anno <- c( 
               'Ex1' = 'DACH2_GRIA4_Ex', 'Ex2' = 'SATB2_SOX6_Ex', 'Ex3' = 'DACH2_GRIA4_Ex', 'Ex4' = 'DACH2_GRIA4_Ex',
               'Ex6' = 'DACH2_GRIA4_Ex', 'Ex7' = 'DACH2_GRIK1_Ex', 'Ex8' = 'DACH2_CCN2_Ex', 'Ex9' = 'SATB2_KIAA1217_Ex',
               'Ex10' = 'DACH2_MEIS2_Ex', 'Ex11' = 'CACNA1H_CCN2_Ex', 'Ex12' = 'SATB2_SOX6_Ex', 'Ex13' = 'DACH2_GRIA4_Ex',
               'Ex14' = 'SATB2_SOX6_Ex', 'Ex15' = 'SATB2_KIAA1217_Ex', 'Ex16' = 'SATB2_SCUBE1_Ex', 'Ex17' = 'SATB2_SATB1_Ex',
               'Ex18' = 'SATB2_RARB_Ex', 'Ex19' = 'SATB2_SATB1_Ex', 'Ex20' = 'CACNA1H_CCN2_Ex', 'Ex21' = 'DACH2_GRIA4_Ex',
               'Ex22' = 'SATB2_RARB_Ex', 'Ex23' = 'CACNA1H_CCN2_Ex', 'Ex24' = 'SATB2_SATB1_Ex', 'Ex25' = 'DACH2_MEIS2_Ex',
               'Ex26' = 'CACNA1H_CCN2_Ex', 'Ex27' = 'SATB2_KIAA1217_Ex', 'Ex28' = 'SATB2_KIAA1217_Ex', 'Ex29' = 'DACH2_GRIK1_Ex'
)
                   

#names(main_type_anno) 
Tel_Glu <- RenameIdents(Tel_Glu, main_type_anno)
Tel_Glu$Glu_level_subclass_annotation <- Idents(Tel_Glu)




#annotation_cluster
Idents(Tel_Glu) <- Tel_Glu$`Glu_integrated_snn_res.4`
main_type_anno <- c( 
'0' = 'Tel_Glu_1', '1' = 'Tel_Glu_2', '2' = 'Tel_Glu_3', '3' = 'Tel_Glu_4', '4' = 'Tel_Glu_5', 
'6' = 'Tel_Glu_6', '7' = 'Tel_Glu_7', '8' = 'Tel_Glu_8', '9' = 'Tel_Glu_9', '10' = 'Tel_Glu_10',
  '12' = 'Tel_Glu_11', '13' = 'Tel_Glu_12', '14' = 'Tel_Glu_13', '15' = 'Tel_Glu_14', '16' = 'Tel_Glu_15', 
    '17' = 'Tel_Glu_16', '18' = 'Tel_Glu_17', '19' = 'Tel_Glu_18', '20' = 'Tel_Glu_19', '21' = 'Tel_Glu_20',
    '22' = 'Tel_Glu_21', '23' = 'Tel_Glu_22', '24' = 'Tel_Glu_23', '25' = 'Tel_Glu_24', '26' = 'Tel_Glu_25',
    '27' = 'Tel_Glu_26', '28' = 'Tel_Glu_27', '29' = 'Tel_Glu_28', '30' = 'Tel_Glu_29', '31' = 'Tel_Glu_30', 
    '32' = 'Tel_Glu_31', '33' = 'Tel_Glu_32', '34' = 'Tel_Glu_33', '35' = 'Tel_Glu_34', '36' = 'Tel_Glu_35',
    '37' = 'Tel_Glu_36', '38' = 'Tel_Glu_37', '39' = 'Tel_Glu_38', '40' = 'Tel_Glu_39', '41' = 'Tel_Glu_40',
    '42' = 'Tel_Glu_41', '43' = 'Tel_Glu_42', '44' = 'Tel_Glu_43', '45' = 'Tel_Glu_44', '46' = 'Tel_Glu_45',
    '47' = 'Tel_Glu_46', '48' = 'Tel_Glu_47', '49' = 'Tel_Glu_48', '50' = 'Tel_Glu_49', '51' = 'Tel_Glu_50',
    '52' = 'Tel_Glu_51', '53' = 'Tel_Glu_52', '54' = 'Tel_Glu_53', '55' = 'Tel_Glu_54', '56' = 'Tel_Glu_55'
)
                   

#names(main_type_anno) 
Tel_Glu <- RenameIdents(Tel_Glu, main_type_anno)
Tel_Glu$Glu_level_cluster <- Idents(Tel_Glu)


mainmarkers <- c('MEIS2','NDST4','LYPD1','BDNF','ZNF831','CHD7','SV2C','CCN2','GRIK1','TAC1',
                 'ZNF385B','GLRA3','RARB','TFAP2D','CCK',
                 'FOXP2', 'SCUBE1','BCL6','KIAA1217',
                 'NRP1','ZBTB20','SLIT1','CACNA1H','SATB2',
                 'DACH2','NTS','SATB1',
                  'GRIA4', 'SOX6','ETV1'
                 )
#plot
p = DotPlot(Tel_Glu, features = mainmarkers, group.by = "Glu_level_subcelltype_annotation", scale = TRUE)
dat <- p$data
ggplot(dat, aes(features.plot, id,size=pct.exp, fill=avg.exp.scaled)) + 
  geom_point(shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill=NA))) + 
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=12),
    axis.text.x = element_text(color='black',size=12, angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_gradientn(colours = c('#7CACCE', '#FFFFFF', '#EA85A8'))
ggsave("/data/work/5month_Tel/Glu/ggplot_Tel_Glu_marker.pdf", width = 10, height = 6)

saveRDS(Tel_Glu, file = '/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_Glu_annotated_!Ex5.rds')