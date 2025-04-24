#edu+/-:single-nucleus; edu2+/-:single-cell
#EdU
EdU_smart <- readRDS('/data/work/neuronal_ipc/EDU/edu_gene.rds')
EdU_smart <- subset(EdU_smart, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA >2000 & nCount_RNA<20000)
table(EdU_smart$Sample)
#adult_data
budgerigar_5months_Tel <- readRDS('/data/work/seurat_integrated/09_budgerigar_brain_5months_Tel_annotated.rds')
downsampled_cells <- c()
annotations <- unique(budgerigar_5months_Tel@meta.data$subcelltype_annotation)
Idents(budgerigar_5months_Tel) <- budgerigar_5months_Tel$`subcelltype_annotation`
for (annotation in annotations) {
    cells_in_cluster <- WhichCells(budgerigar_5months_Tel, ident = annotation)
    if (length(cells_in_cluster) > 100) {
        set.seed(42) 
        sampled_cells <- sample(cells_in_cluster, 100)
    } else {
        sampled_cells <- cells_in_cluster
    }
    
    downsampled_cells <- c(downsampled_cells, sampled_cells)
}

seurat_downsampled <- subset(budgerigar_5months_Tel, cells = downsampled_cells)
DefaultAssay(seurat_downsampled) <- "RNA"
features <- rownames(EdU_smart)
# find anchors
anchors <- FindTransferAnchors(reference = seurat_downsampled, query = EdU_smart,features = features)
# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = seurat_downsampled$subcelltype_annotation)
EdU_smart <- AddMetaData(object = EdU_smart, metadata = predictions)
# standard analysis
EdU_smart <- NormalizeData(EdU_smart, normalization.method = "LogNormalize", scale.factor = 10000)
EdU_smart <- FindVariableFeatures(EdU_smart, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(EdU_smart), 10)
all.genes <- rownames(EdU_smart)
EdU_smart <- ScaleData(EdU_smart, features = all.genes)
EdU_smart <- RunPCA(EdU_smart, features = VariableFeatures(object = EdU_smart))
EdU_smart <- FindNeighbors(EdU_smart, dims = 1:30)
EdU_smart <- FindClusters(EdU_smart, resolution = 0.4)
head(Idents(EdU_smart), 5)
EdU_smart <- RunUMAP(EdU_smart, dims = 1:30)
VlnPlot(EdU_smart, features = c("nCount_RNA","nFeature_RNA"), group.by = "Sample",ncol=2)

cell_colors <- c(
'DACH2_CHD7_Ex' = '#807087', 
'SATB2_FOXP2_Ex' = '#E45A5F',
'DACH2_SV2C_Ex' = '#e1a4c6',   
'DACH2_CCN2_Ex' = '#d2e0ac',    
'SATB2_BCL6_Ex' = '#74a893',
'DACH2_MEIS2_Ex' = '#7587b1',   
'CACNA1H_SLIT1_Ex' = '#edeaa4',  
'SATB2_CCK_Ex' = '#e3d1db',
'SATB2_RARB_Ex' = '#e0bc58',
"SATB2_GLRA3_Ex" = "#FEE082",
"DACH2_TAC1_Ex" = "#8C33FF",
 'LGE_FOXP2_SLIT2_Inh' = '#FF6B95',   
'SNCG_NRP1_Inh' = '#B2E5FB',
'LGE_PCP4_EBF1_Inh' = '#c7deef',
'MGE_SST_NPY_Inh' = '#EBEBF3',   
'LGE_PCP4_IL17REL_Inh' = "#5EC0AD",   
"LGE_PCP4_SLC1A3_Inh" = "#FEE082",   
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF"
)
DimPlot(EdU_smart,group.by = 'predicted.id',label = F,repel = F,cols = cell_colors)
ggsave("/data/work/neuronal_ipc/EDU/edu_predicted.id.pdf",bg = "transparent", width =9.5, height = 6)

#EdU positive excitatory neurons
EdU_smart <- subset(EdU_smart, subset = Sample != "SHC-EDU--index8")
edu_P <- subset(EdU_smart, subset = Sample != "SHC-EDU2--index7")
edu_P <- subset(edu_P, subset = predicted.id %in% c("DACH2_CCN2_Ex", "DACH2_MEIS2_Ex",  
 'DACH2_CHD7_Ex', 
'CACNA1H_SLIT1_Ex',  
 'DACH2_SV2C_Ex',     
'SATB2_CCK_Ex',
'SATB2_FOXP2_Ex'))
DACH2_MEIS2_Ex_genes <- c( "VCAN",
                           "NR2F2", "SYNE2",
                           "EPHA5",
                            "ZBTB18", 
                            "DCX" 
                          )
#AddModuleScore
edu_P <- AddModuleScore(edu_P, features = list(DACH2_MEIS2_Ex_genes), name = "DACH2_MEIS2_Ex_Score", ctrl = 100)
# plot
p1 <- ggplot()
vln_data <- FetchData(edu_P, vars = c("predicted.id", "DACH2_MEIS2_Ex_Score1"))
vln_data$group <- ifelse(vln_data$predicted.id == "DACH2_MEIS2_Ex", "DACH2_MEIS2_Ex", "Other")
t_test_result <- wilcox.test(DACH2_MEIS2_Ex_Score1 ~ group, data = vln_data)
print(t_test_result)
library(ggplot2)
library(ggpubr)
p1 <- ggplot(vln_data, aes(x = group, y = DACH2_MEIS2_Ex_Score1, fill = group)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("DACH2_MEIS2_Ex" = "#7587b1", "Other" = "gray")) +
  labs(x = "Group", y = "DACH2_MEIS2_Ex_Score1") +
  theme_classic() +
  theme(legend.position = "none")

p1 + stat_compare_means(
  method = "wilcox.test", 
  comparisons = list(c("DACH2_MEIS2_Ex", "Other")),
  label = "p.signif", 
  label.y = max(vln_data$DACH2_MEIS2_Ex_Score1) * 1.1 
)
ggsave("/data/work/neuronal_ipc/EDU/DACH2_MEIS2_Ex_Score1_ggplot2.pdf", width = 5, height = 5)