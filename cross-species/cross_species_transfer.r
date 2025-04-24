#transfer chicken brain(Zaremba et al., 2025)
library("Seurat")
library("dplyr")
library('ggplot2')
library("tidyverse")
chicken <- readRDS('Gg_adult_snRNA_seg_srt')
Tel_Glu <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_Glu_annotated_!Ex5.rds')
common_genes <- rownames(chicken)[rownames(Tel_Glu) %in% rownames(chicken)]
length(common_genes)
chicken <- chicken[rownames(chicken) %in% common_genes,]
Tel_Glu <- Tel_Glu[rownames(Tel_Glu) %in% common_genes,]
DefaultAssay(chicken) <- "RNA"
features <- rownames(chicken)
# find anchors
anchors <- FindTransferAnchors(reference = Tel_Glu, query = chicken,features = features)
# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = Tel_Glu$subcelltype_annotation)
chicken <- AddMetaData(object = chicken, metadata = predictions)
saveRDS(chicken, file = '')
write.csv(predictions,"metadata_v5.csv")


