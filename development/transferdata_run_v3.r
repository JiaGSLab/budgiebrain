library("Seurat")
budgerigar_dev_Tel_Glu <- readRDS('')
Tel_Glu <- readRDS('10_budgerigar_brain_5months_Tel_Glu_annotated_notEx5.rds')
Idents(budgerigar_dev_Tel_Glu) <- budgerigar_dev_Tel_Glu$`harmony_snn_res.0.9`
main_type_anno <- c( 
'0' = 'Progenitor_Ex', '1' = 'Nido/Hyper_GRIK1_Ex', '2' = 'Nido/Hyper_GRIA4_Ex', '3' = 'Pallium_Ipc', '4' = 'Meso_GAS2_Ipc', 
    '5' = 'Meso_SOX6_Ex', '6' = 'Arco_Ex', '7' = 'Meso_SCUBE1_Ex', '8' = 'RFX6_ITGA4_Ex', '9' = 'Meso_SOX6_Ex', 
    '10' = 'Meso_KIAA1217_Ex', '11' = 'Meso_RARB_Ex', '12' = 'Meso_KIAA1217_Ex', '13' = 'Meso_RARB_Ex', '14' = 'Hp_Ex', 
    '15' = 'Ento_Ex', '16' = 'NLc_Ex', '17' = 'Nido/Hyper_GRIK1_Ex', '18' = 'Ento_Ex', '19' = 'Meso_SOX6_Immature_Ex' 
)                  
#names(main_type_anno) 
budgerigar_dev_Tel_Glu <- RenameIdents(budgerigar_dev_Tel_Glu, main_type_anno)
budgerigar_dev_Tel_Glu$dev_celltype <- Idents(budgerigar_dev_Tel_Glu)
DefaultAssay(budgerigar_dev_Tel_Glu) <- "RNA"
features <- rownames(Tel_Glu)
# find anchors
anchors <- FindTransferAnchors(reference = budgerigar_dev_Tel_Glu, query = Tel_Glu,features = features)
# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = budgerigar_dev_Tel_Glu$dev_celltype)
Tel_Glu <- AddMetaData(object = Tel_Glu, metadata = predictions)
saveRDS(Tel_Glu, file = 'Tel_dev_Glu_subcellype_annotation_v3.rds')
write.csv(predictions,"metadata_v3.csv")
