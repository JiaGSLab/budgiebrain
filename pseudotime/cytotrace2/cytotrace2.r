library(CytoTRACE2)
library(Seurat)
Tel_Glu <- readRDS('/data/work/Develop/12_budgerigar_dev_Tel_Glu_harmony.rds')
Idents(Tel_Glu) <- Tel_Glu$`dev_celltype`
result_sce <- cytotrace2(Tel_Glu,
                    is_seurat = TRUE,
                    slot_type = "counts",
                    species = "human",
                    seed = 1234)
annotation <- data.frame(phenotype=result_sce@meta.data$dev_celltype) %>% set_rownames(., colnames(result_sce))
CytoTRACE2_Score <- as.data.frame(result_sce$CytoTRACE2_Score)
write.csv(CytoTRACE2_Score, file = "/data/work/Develop/cytotrace/Tel_Glu_dev_cytotrace2_result.csv")
CytoTRACE2_Potency <- as.data.frame(result_sce$CytoTRACE2_Potency)
CytoTRACE2_Relative <- as.data.frame(result_sce$CytoTRACE2_Relative)
preKNN_CytoTRACE2_Score <- as.data.frame(result_sce$preKNN_CytoTRACE2_Score)
preKNN_CytoTRACE2_Potency <- as.data.frame(result_sce$preKNN_CytoTRACE2_Potency)
write.csv(CytoTRACE2_Potency, file = "/data/work/Develop/cytotrace/Tel_Glu_dev_CytoTRACE2_Potency.csv")
write.csv(CytoTRACE2_Relative, file = "/data/work/Develop/cytotrace/Tel_Glu_dev_CytoTRACE2_Relative.csv")
write.csv(preKNN_CytoTRACE2_Score, file = "/data/work/Develop/cytotrace/Tel_Glu_dev_preKNN_CytoTRACE2_Score.csv")
write.csv(preKNN_CytoTRACE2_Potency, file = "/data/work/Develop/cytotrace/Tel_Glu_dev_preKNN_CytoTRACE2_Potency.csv")