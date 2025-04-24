library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")
budgerigar_brain_obj <- readRDS('/data/work/seurat_integrated/07_budgerigar_brain_obj_unannotated.rds')
Idents(budgerigar_brain_obj) <- budgerigar_brain_obj$`integrated_snn_res.0.4`

main_type_anno <- c( '0' = 'Glutamatergic',
                           '1' = 'GABAergic',
                           '2' = 'Glutamatergic',
                           '3' = 'Astrocyte',
                           '4' = 'Glutamatergic',
                           '5' = 'Oligodendrocyte',
                           '6' = 'MSN',
                           '7' = 'Glutamatergic',
                           '8' = 'Glutamatergic',
                           '9' = 'GABAergic',
                           '10' = 'Glutamatergic',
                           '11' = 'Glutamatergic',
                           '12' = 'GABAergic',
                           '13' = 'OPC',
                           '14' = 'MSN',
                           '15' = 'Glutamatergic',
                           '16' = 'Glutamatergic',
                           '17' = 'GABAergic',
                           '18' = 'Microglia',
                           '19' = 'Glutamatergic',
                           '20' = 'Vascular cell',
                           '21' = 'APC',
                           '22' = 'OPC',
                           '23' = 'Mural cell',
                           '24' = 'GABAergic',
                           '25' = 'Blood',
                           '26' = 'MSN',
                           '27' = 'Purkinje',
                           '28' = 'Astroependymal',
                           '29' = 'Glutamatergic',
                           '30' = 'Astroependymal',
                           '31' = 'Glutamatergic',
                           '32' = 'Astrocyte',
                           '33' = 'MSN',
                           '34' = 'GABAergic')
                   

#names(main_type_anno) 
budgerigar_brain_obj <- RenameIdents(budgerigar_brain_obj, main_type_anno)
budgerigar_brain_obj$main_cell_type_res_0.4 <- Idents(budgerigar_brain_obj)

region_colors <- c(
   'Diencephalon'= '#C2C1E0',
    'Cerebellum' ='#ffd92f',
    'Midbrain' = '#D3E2B7',
    'Telencephalon' = '#FFC0CB'
)
DimPlot(budgerigar_brain_obj, label = FALSE, reduction = "umap",raster=TRUE,group.by ="Region_1", alpha = 0.2,cols = region_colors)
ggsave("/data/work/seurat_integrated/region_umap.pdf", width = 30, height = 10)

age_colors <- c(
   'E14'= '#D55D4C',
    'P01' ='#ffd92f',
    '2months' = '#D3E2B7',
    '5months' = '#B9E8EA'
)
DimPlot(budgerigar_brain_obj, label = FALSE, reduction = "umap",raster=TRUE,group.by ="Age", alpha = 0.2,cols = age_colors)
ggsave("/data/work/seurat_integrated/age_umap.pdf")

cell_colors <- c(
  "MSN" = "#d2e0ac",
  "APC" = "#c7deef",
  "Purkinje" = "#F3BAA5",
'Glutamatergic' = '#E84e40', 
"GABAergic" = "#A4D38E",
"Astrocyte" = "#82D4FB",
"Oligodendrocyte" = "#FF33A1",
"OPC" = "#FEE082",
"Microglia" = "#7D33FF",
"Vascular cell" = "#EF859B",
"Mural cell" = "#FFD133",
"Blood" = "#8C33FF",
"Purkinje" = "#F3BAA5",
"Astroependymal" = "#33A1FF"
)
DimPlot(budgerigar_brain_obj, label = FALSE, reduction = "umap",raster=TRUE,group.by ="main_cell_type_res_0.4", alpha = 0.2,cols = cell_colors)