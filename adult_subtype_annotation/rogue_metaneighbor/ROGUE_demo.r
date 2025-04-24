library("tidyverse")
library(ROGUE)
library(ggplot2)
Tel_Glu <- readRDS('/data/work/5month_Tel/10_budgerigar_brain_5months_Tel_Glu_annotated.rds')
DefaultAssay(Tel_Glu) <- "RNA"
# downsample
total_cells <- ncol(Tel_Glu)
print(paste("Total cells in dataset:", total_cells))
if (total_cells > 10000) {
  set.seed(123) # seed
  sampled_cells <- sample(colnames(Tel_Glu), 10000)
  Tel_Glu_subset <- subset(Tel_Glu, cells = sampled_cells)
} else {
  Tel_Glu_subset <- Tel_Glu
}

resolutions <- c(1, 2, 3, 4, 5, 6)
columns <- paste0("Glu_integrated_snn_res.", resolutions)
average_rogue_values <- numeric(length(columns))
for (i in seq_along(columns)) {
  column <- columns[i]
  if (column %in% colnames(Tel_Glu_subset@meta.data)) {
    expr <- GetAssayData(Tel_Glu_subset, slot = 'counts') %>% as.matrix()
    meta <- Tel_Glu_subset@meta.data
    
    # check NA
    if (any(is.na(meta[[column]]))) {
      warning(paste("Column", column, "contains NA values. Skipping this column."))
      average_rogue_values[i] <- NA
    } else {
      span_value <- 0.8  
      
      rogue_result <- rogue(expr = expr, 
                            labels = meta[[column]], 
                            samples = meta$Sample_ID, 
                            platform = 'UMI',
                            span = span_value)
      
      column_means <- apply(rogue_result, 2, mean, na.rm = TRUE)
      average_rogue_values[i] <- mean(column_means, na.rm = TRUE)
    }
  } else {
    average_rogue_values[i] <- NA
  }
}
#plot
rogue_df <- data.frame(
  Clusters = columns,
  Average_ROGUE = average_rogue_values
)
ggplot(rogue_df, aes(x = Clusters, y = Average_ROGUE, group = 1)) +
  geom_line(color = "#ff1691") +
  geom_point(color = "#ff1691") +
  theme_bw() +
  labs(
    title = "Average ROGUE for Different Clustering Resolutions",
    x = "Clustering Resolution",
    y = "Average ROGUE"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/data/work/5month_Tel/Glu_rogue/TeL_Glu_1-6.pdf",bg = "transparent", width = 6, height = 5.5)

write.csv(rogue_df, file = "/data/work/5month_Tel/Glu_rogue/Tel_Glu_1-6.csv")